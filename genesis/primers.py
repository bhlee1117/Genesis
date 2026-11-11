"""Stage 9: junctions, primer design and primer validation.

Primers are built as ``[overlap tail] + [annealing region]``. Only the
annealing region determines the extension temperature, so Tm is always
reported for that part alone - a very common source of failed PCRs is
quoting the Tm of the whole oligo including its Gibson tail.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict

from .config import Config
from .seqtools import (clean_seq, gc_fraction, hairpin_score, max_homopolymer,
                       revcomp, self_dimer_3prime, slice_region, tm, translate)


@dataclass
class Primer:
    name: str
    sequence: str
    anneal: str                 # template-binding 3' portion
    tail: str = ""              # overlap / adapter 5' portion
    direction: str = "F"
    template: str = ""
    target_start: int = 0
    target_end: int = 0
    warnings: list = field(default_factory=list)
    #: For a reverse primer, how many bases its anchor was pulled back
    #: from the junction so the pair could share an annealing step. The
    #: next fragment's forward tail has to be this much longer, or those
    #: bases appear in no amplicon and the product comes out short.
    retreat: int = 0
    #: Enough of the template to work out what this primer actually binds,
    #: kept because the tail is rewritten after the primer is built - to
    #: carry a folded element, a stretch of recipient, or a neighbour's
    #: retreat - and anything precomputed against the old tail would be
    #: quietly wrong afterwards. Dropped from `to_dict`, which is what
    #: gets stored and sent.
    tpl_seq: str = ""
    tpl_circular: bool = False
    anneal_at: int = 0

    @property
    def tail_binds(self) -> int:
        """How many bases at the 3' end of `tail` also pair with the template.

        An overlap arm is the neighbouring fragment's sequence, so in
        principle none of it binds here - but junctions fall at Kozaks
        and start codons, where both sides of the join look alike, and
        the last bases of the arm often match anyway. They pair in the
        first cycle like any others, so they belong in the Tm: on the
        iGluSnFR4f forward primer two such bases are worth 4.0 C, which
        was the difference between a pair reported 4.4 C apart and one
        that is really 0.4 C apart.

        Computed on demand rather than stored, because the tail changes
        under it.
        """
        if not self.tpl_seq or not self.tail:
            return 0
        return tail_match(self.tail, self.tpl_seq, self.anneal_at,
                          self.direction, self.tpl_circular)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def binding(self) -> str:
        """Everything at this primer's 3' end that pairs with the template."""
        return (self.tail[len(self.tail) - self.tail_binds:] + self.anneal
                if self.tail_binds else self.anneal)

    @property
    def arm(self) -> str:
        """The part of the tail that does not pair with the template."""
        return (self.tail[:len(self.tail) - self.tail_binds]
                if self.tail_binds else self.tail)

    def to_dict(self, cfg: Config | None = None) -> dict:
        cfg = cfg or Config()
        d = asdict(self)
        d.pop("tpl_seq", None)            # kilobases, and only needed here
        d["length"] = self.length
        d["tail_binds"] = self.tail_binds
        d["binding"] = self.binding
        d["arm"] = self.arm
        d["anneal_tm"] = tm(self.binding, cfg.primer_conc_nm,
                            cfg.na_mm, cfg.mg_mm)
        d["full_tm"] = tm(self.sequence, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
        d["gc"] = round(gc_fraction(self.binding), 3)
        return d


@dataclass
class Junction:
    """A boundary between two modules in the final construct."""
    index: int
    left_name: str
    right_name: str
    left_type: str
    right_type: str
    method: str                  # gibson | hifi | ligation | none
    sequence: str = ""           # the junction as it appears in the product
    overlap_length: int = 0
    frame_shift: int = 0
    added_residues: str = ""
    #: Whether two amplicons actually meet here, as opposed to this
    #: being a boundary between two elements inside one of them.
    is_join: bool = True
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------------------------------------
# Annealing-region selection
# --------------------------------------------------------------------------

def anneal_penalty(region: str, cfg: Config,
                   unique_in: str = "", circular: bool = False) -> float:
    """How poor a primer this stretch of template would make. Lower is better.

    The criteria are NEB's for Q5: a melting temperature near the target,
    40-60% GC, a G or C at the 3' end but not four of them in the last
    five, no long single-base run, and no 3' self-complementarity. None of
    them is a veto on its own - a site that fails one and passes the rest
    beats a site that fails three - so they are weighed rather than
    filtered.
    """
    if not region:
        return 1e6
    t = tm(region, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
    if t != t:                                  # ambiguous bases
        return 1e6
    score = abs(t - cfg.primer_tm_target)
    if abs(t - cfg.primer_tm_target) > cfg.primer_tm_tolerance:
        score += 6.0

    gc = gc_fraction(region)
    if gc < cfg.primer_ideal_min_gc:
        score += (cfg.primer_ideal_min_gc - gc) * 40
    elif gc > cfg.primer_ideal_max_gc:
        score += (gc - cfg.primer_ideal_max_gc) * 40

    # A G or C at the 3' end clamps the extending end; four of them in the
    # last five is the opposite problem and misprimes.
    if region[-1] in "GC":
        score -= 1.5
    if region[-5:].count("G") + region[-5:].count("C") >= 4:
        score += 4.0

    if max_homopolymer(region) >= 5:
        score += 5.0
    if self_dimer_3prime(region) >= 4:
        score += 5.0
    hair = hairpin_score(region)
    if hair and hair["stem"] >= 6:
        score += 4.0

    if unique_in:
        from .seqtools import find_all
        n = len(find_all(unique_in, region, circular)) + \
            len(find_all(unique_in, revcomp(region), circular))
        if n > 1:
            score += 25.0 * (n - 1)
    return score


#: What a degree of mismatch between a pair costs, per degree past the
#: limit. Large enough that no combination of the single-primer criteria
#: outweighs it, because a pair that cannot share an annealing
#: temperature is not a pair.
PAIR_PENALTY = 8.0

#: What one base of retreat costs. Small, because a dozen extra bases on
#: one oligo is a far better trade than a pair that cannot share an
#: annealing step - but not zero, or the search wanders off the junction
#: for a tenth of a degree.
RETREAT_PENALTY = 0.15


def rank_anneal_sites(template: str, start: int, direction: str, cfg: Config,
                      circular: bool = False, unique_in: str = "",
                      keep: int = 0, tail: str = "") -> list[dict]:
    """Every site worth considering, best first, with the Tm it would have.

    `pick_anneal` answers "which is the best site for this primer", which
    is the wrong question when the answer has to agree with the other
    primer's.  This returns the shortlist so the two can be chosen
    together.
    """
    out: list[dict] = []
    for length in range(cfg.primer_min_len, cfg.primer_max_len + 1):
        if direction == "F":
            region = slice_region(template, start, start + length, circular)
        else:
            lo = start - length
            if lo < 0 and not circular:
                break
            region = revcomp(slice_region(template, lo, start, circular))
        if len(region) < length:
            break
        # The arm this primer will carry may go on matching the template
        # past the junction, and those bases pair too. Rank on what the
        # finished oligo will do, not on the region in isolation: that
        # mismatch put eleven of forty-five pairs outside the limit while
        # every one of them had been chosen to be inside it.
        ext = tail_match(tail, template, start, direction, circular) if tail else 0
        binding = tail[len(tail) - ext:] + region if ext else region
        t = tm(binding, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
        if t != t:
            continue
        out.append({"region": region, "tm": t,
                    "penalty": anneal_penalty(binding, cfg, unique_in,
                                              circular)})
    # Sorted for the single-primer caller; NOT truncated by default,
    # because the pair search needs the cool candidates a ranking by
    # single-primer penalty would cut off.  The pool is one length
    # window, so keeping it whole costs nothing and a truncated one
    # silently hid the only sites that could have matched the other
    # primer.
    out.sort(key=lambda c: c["penalty"])
    return out[:keep] if keep else out


def pair_anneal_sites(template: str, start: int, end: int, cfg: Config,
                      circular: bool = False, unique_in: str = "",
                      max_retreat: int = 0, f_tail: str = "",
                      r_tail: str = "") -> tuple:
    """Choose both sites of one reaction together.

    A forward and a reverse primer go into the same tube and anneal at
    the same temperature, so the pair is what has to be right. Choosing
    each one's best site on its own optimises a quantity nobody runs:
    it put 12 of 13 pairs on this inventory outside the 5 C that NEB's
    calculator allows, one of them by 19 C, because the two sites were
    picked from different stretches of sequence with no reference to
    each other.

    `max_retreat` lets the reverse anchor move back into the fragment.
    The forward primer cannot move - the bases it would skip have nowhere
    to go, because reverse primers here carry no tail - but the reverse
    one can, because whatever it stops short of is written onto the NEXT
    fragment's forward tail instead. That is the only free variable at a
    junction, and it is usually enough: the backbone junction that sat
    6.7 C apart sits 0.5 C apart twelve bases earlier.

    The junction only moves to fix a mismatch, never to shave a tenth of
    a degree off a pair that already matches - a retreat costs the
    neighbouring oligo a base for every base of it, and the first version
    of this spent eight of them on a 0.1 C gap. So the pinned junction is
    tried first and kept if it works, and otherwise the nearest position
    that brings the pair inside the limit wins.

    Returns (f_region, r_region, gap, retreat).
    """
    limit = getattr(cfg, "primer_pair_max_diff", 5.0)
    fwd = rank_anneal_sites(template, start, "F", cfg, circular, unique_in,
                            tail=f_tail)
    fallback = (pick_anneal(template, start, "F", cfg, circular),
                pick_anneal(template, end, "R", cfg, circular),
                float("nan"), 0)
    if not fwd:
        return fallback

    def best_against(rev: list):
        """The best pairing of this reverse shortlist with the forward one."""
        chosen, cost = None, float("inf")
        for a in fwd:
            for b in rev:
                gap = abs(a["tm"] - b["tm"])
                c = (a["penalty"] + b["penalty"]
                     + PAIR_PENALTY * max(0.0, gap - limit))
                if c < cost:
                    cost, chosen = c, (a, b, gap)
        return chosen

    here = best_against(
        rank_anneal_sites(template, end, "R", cfg, circular, unique_in,
                          tail=r_tail))
    if here is None:
        return fallback
    if here[2] <= limit or max_retreat <= 0:
        return here[0]["region"], here[1]["region"], here[2], 0

    # The junction has to move. Take the first position that works, and
    # failing that the position that comes closest.
    nearest = (here[2], here, 0)
    for back in range(1, max_retreat + 1):
        rev = rank_anneal_sites(template, end - back, "R", cfg, circular,
                                unique_in, tail=r_tail)
        if not rev:
            continue
        cand = best_against(rev)
        if cand is None:
            continue
        if cand[2] <= limit:
            return cand[0]["region"], cand[1]["region"], cand[2], back
        if cand[2] < nearest[0]:
            nearest = (cand[2], cand, back)
    _gap, best, back = nearest
    a, b, gap = best
    return a["region"], b["region"], gap, back


def pick_anneal(template: str, start: int, direction: str, cfg: Config,
                circular: bool = False) -> str:
    """Choose the region of the primer that hybridises to the template.

    Length is held inside the window Addgene specifies for PCR cloning
    (18-21 bp by default), and within that window the length whose melting
    temperature sits closest to the target is taken. Only this part of the
    oligo is measured: the Tm that matters is "the Tm of the portion of the
    primer that hybridizes to the sequence to be amplified", not the Tm of
    the whole thing with its tail.

    `direction` 'F' extends rightwards from `start`; 'R' extends leftwards
    and the returned sequence is already reverse-complemented.
    """
    best = ""
    best_delta = float("inf")
    for length in range(cfg.primer_min_len, cfg.primer_max_len + 1):
        if direction == "F":
            region = slice_region(template, start, start + length, circular)
        else:
            left = start - length
            if left < 0 and not circular:
                break
            region = revcomp(slice_region(template, left, start, circular))
        if len(region) < length:
            break
        t = tm(region, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
        if t != t:  # NaN from ambiguous bases
            continue
        delta = abs(t - cfg.primer_tm_target)
        # Prefer a G or C at the 3' end: it clamps the extending end.
        if region[-1] in "GC":
            delta -= 0.4
        if delta < best_delta:
            best_delta, best = delta, region
        if t > cfg.primer_tm_target + cfg.primer_tm_tolerance:
            break
    return best or slice_region(template, start, start + cfg.primer_min_len, circular)


def validate_primer(primer: Primer, template: str, cfg: Config,
                    circular: bool = False) -> list[str]:
    """Run the Stage 9 primer checks and return human-readable warnings."""
    w: list[str] = []
    anneal = primer.anneal
    full = primer.sequence
    # The Tm check asks what the oligo does on the template, so it runs
    # on everything that pairs; the length checks below ask what the
    # designer chose, so they stay on the region that was designed.
    binding = primer.binding

    t = tm(binding, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
    if t != t:
        n_ambig = sum(1 for ch in binding if ch not in "ACGT")
        w.append(f"the source sequence has {n_ambig} ambiguous base(s) at this "
                 f"primer site, so no reliable primer can be designed here; "
                 f"re-sequence this region or choose another source plasmid")
    else:
        if t < cfg.primer_tm_target - cfg.primer_tm_tolerance:
            w.append(f"annealing Tm {t:.1f} C is below the {cfg.primer_tm_target:.0f} C target")
        elif t > cfg.primer_tm_target + cfg.primer_tm_tolerance + 4:
            w.append(f"annealing Tm {t:.1f} C is well above the {cfg.primer_tm_target:.0f} C target")

    if not (cfg.primer_min_len <= len(anneal) <= cfg.primer_max_len):
        w.append(f"the hybridising region is {len(anneal)} nt, outside the "
                 f"{cfg.primer_min_len}-{cfg.primer_max_len} nt window this "
                 f"protocol uses")
    if len(full) > 60:
        w.append(f"oligo is {len(full)} nt; synthesis cost and error rate rise "
                 f"sharply past ~60 nt")
    if len(anneal) < cfg.primer_min_len:
        w.append(f"annealing region is only {len(anneal)} nt")

    gc = gc_fraction(anneal)
    if gc < cfg.primer_min_gc:
        w.append(f"GC content {gc:.0%} is low")
    elif gc > cfg.primer_max_gc:
        w.append(f"GC content {gc:.0%} is high")

    hp = max_homopolymer(full)
    if hp >= 5:
        w.append(f"{hp} bp homopolymer run")

    hair = hairpin_score(full)
    if hair and hair["stem"] >= 6:
        w.append(f"possible hairpin with a {hair['stem']} bp stem")

    dimer = self_dimer_3prime(full)
    if dimer >= 4:
        w.append(f"{dimer} nt self-complementarity at the 3' end; primer-dimer risk")

    if anneal[-5:].count("G") + anneal[-5:].count("C") >= 4:
        w.append("GC-rich 3' end (4+ of the last 5 bases); may misprime")

    # Mispriming: does the annealing region occur more than once?
    if template:
        from .seqtools import find_all
        n_fwd = len(find_all(template, anneal, circular))
        n_rev = len(find_all(template, revcomp(anneal), circular))
        if n_fwd + n_rev > 1:
            w.append(f"annealing region occurs {n_fwd + n_rev} times in the "
                     f"template; expect mispriming")
        # A weaker check on just the 3' half, which drives priming.
        half = anneal[len(anneal) // 2:]
        if len(half) >= 12:
            n_half = len(find_all(template, half, circular)) + \
                len(find_all(template, revcomp(half), circular))
            if n_half > 3:
                w.append(f"3' half matches {n_half} sites in the template")
    return w


# --------------------------------------------------------------------------
# Gibson / HiFi primer pairs
# --------------------------------------------------------------------------

def design_gibson_pair(fragment_index: int, fragment_name: str,
                       template: str, start: int, end: int,
                       upstream_seq: str, downstream_seq: str,
                       cfg: Config, circular_template: bool = False,
                       carry_downstream: bool = True,
                       unique_in: str = "",
                       max_retreat: int = 0, extra_arm: int = 0) -> tuple[Primer, Primer]:
    """Design the primer pair that amplifies one fragment with HiFi overlaps.

    The forward primer carries the 3' end of the *upstream* fragment as a
    tail.  `carry_downstream` adds the mirror image to the reverse primer.

    A junction needs its homology written once.  When every fragment in the
    assembly is amplified, the next fragment's forward primer already
    carries this one's end, and adding it here too means each amplicon
    finishes with its neighbour's start instead of its own: the two ends
    that are supposed to anneal then share nothing, and the assembly does
    not close.  Pass False in that case.  Keep it True when the neighbour is
    not amplified - a digested vector cannot be given a tail, so the insert
    has to carry both of its junctions.
    """
    # Chosen together. The two primers of one reaction anneal in the same
    # tube at the same temperature, so the pair is the unit that has to be
    # right; picking each one's best site alone optimises a quantity
    # nobody runs.
    # Tails first. They are the neighbours' sequence and do not depend on
    # where this fragment's primers land - but where the primers land has
    # to depend on them, because an arm that goes on matching the template
    # is part of the duplex.
    f_tail = pick_overlap(upstream_seq, cfg, from_end=True,
                          unique_in=unique_in)
    # `extra_arm` is the previous fragment's retreat: the bases it stopped
    # short of, which exist in no amplicon unless this tail writes them.
    # They sit at the END of the upstream sequence, which is where the arm
    # is taken from, so reaching further back picks them up and leaves the
    # homology that actually anneals the same length as before.
    if extra_arm and upstream_seq:
        want = len(f_tail) + extra_arm
        f_tail = upstream_seq[-want:] if want <= len(upstream_seq) \
            else upstream_seq
    r_tail = (revcomp(pick_overlap(downstream_seq, cfg, from_end=False,
                                   unique_in=unique_in))
              if carry_downstream else "")

    f_anneal, r_anneal, tm_gap, retreat = pair_anneal_sites(
        template, start, end, cfg, circular_template, unique_in=template,
        max_retreat=max_retreat, f_tail=f_tail, r_tail=r_tail)

    fwd = Primer(
        name=f"{fragment_name}-F", sequence=f_tail + f_anneal, anneal=f_anneal,
        tail=f_tail, direction="F", template=fragment_name,
        target_start=start, target_end=end)
    rev = Primer(
        name=f"{fragment_name}-R", sequence=r_tail + r_anneal, anneal=r_anneal,
        tail=r_tail, direction="R", template=fragment_name,
        target_start=start, target_end=end - retreat, retreat=retreat)

    # Where each primer sits, so `binding` can be worked out later even
    # though the tails get rewritten between here and the order sheet.
    for pr, at in ((fwd, start), (rev, end - retreat)):
        pr.tpl_seq, pr.tpl_circular, pr.anneal_at = (
            template, circular_template, at)

    fwd.warnings = validate_primer(fwd, template, cfg, circular_template)
    rev.warnings = validate_primer(rev, template, cfg, circular_template)

    # The gap the chooser worked with was computed before the tails
    # existed, so it knew nothing about this. Re-measure.
    _a = tm(fwd.binding, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
    _b = tm(rev.binding, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
    if _a == _a and _b == _b:
        tm_gap = abs(_a - _b)

    limit = getattr(cfg, "primer_pair_max_diff", 5.0)
    if tm_gap == tm_gap and tm_gap > limit:
        for pr in (fwd, rev):
            pr.warnings.append(
                f"this pair's Tms differ by {tm_gap:.1f} C, more than the "
                f"{limit:.0f} C two primers sharing one annealing step "
                f"should; no combination of lengths"
                + (f" or of junction positions up to {max_retreat} bp"
                   if max_retreat else "")
                + " does better, so run a gradient over the lower Tm")

    # The tail is homology and nothing else now, so all of it is the arm.
    for p, tail in ((fwd, f_tail), (rev, r_tail)):
        if tail:
            t_ov = tm(tail, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
            if t_ov == t_ov and t_ov < cfg.overlap_tm_min:
                p.warnings.append(
                    f"overlap tail Tm {t_ov:.1f} C is below the "
                    f"{cfg.overlap_tm_min:.0f} C needed for efficient HiFi joining")
    return fwd, rev


#: Alternatives tried when the default leader would create a site of its
#: own. All are 6 bp with no run longer than two, so they do not trip the
#: homopolymer check either.
_LEADERS = ("TTAATT", "ATCTAT", "TACTAC", "AGTAGT", "CATCAT", "GTAGTA")


def choose_leader(site: str, cfg: Config, stocked: tuple = ()) -> str:
    """Pick the extra bases that sit 5' of a restriction site.

    Addgene calls for 3-6 bp so the enzyme has something to grip. The choice
    is not arbitrary: a leader that happens to complete another recognition
    sequence where it abuts the site would put an extra cut in the amplicon,
    so each candidate is checked against the site it precedes and against
    every enzyme the lab stocks.
    """
    from .enzymes import CATALOG
    from .seqtools import find_all

    candidates = [cfg.primer_leader] + [x for x in _LEADERS
                                        if x != cfg.primer_leader]
    names = stocked or tuple(CATALOG)
    for leader in candidates:
        junction = leader + site
        clash = False
        for name in names:
            enzyme = CATALOG.get(name)
            if enzyme is None:
                continue
            # Count distinct cut positions. A palindromic site is found by
            # both the forward and the reverse-complement search at the same
            # place, and counting that as two would reject every leader.
            positions = set(find_all(junction, enzyme.site))
            positions |= set(find_all(junction, revcomp(enzyme.site)))
            expected = 1 if enzyme.site in (site, revcomp(site)) else 0
            if len(positions) > expected:
                clash = True
                break
        if not clash:
            return leader
    return cfg.primer_leader


def design_re_pair(fragment_name: str, template: str, start: int, end: int,
                   left_enzyme: str, right_enzyme: str, cfg: Config,
                   circular_template: bool = False) -> tuple[Primer, Primer]:
    """Primers that add restriction sites, with a leader so they can be cut.

    Structure follows Addgene's PCR cloning protocol: a 3-6 bp leader, then
    the 6-8 bp recognition site, then the 18-21 bp hybridising region.
    """
    from .enzymes import CATALOG
    # Same reasoning as the assembly pair: both go in one tube.
    f_anneal, r_anneal, tm_gap, _retreat = pair_anneal_sites(
        template, start, end, cfg, circular_template)
    f_tail = r_tail = ""
    if left_enzyme:
        site = CATALOG[left_enzyme].site
        f_tail = choose_leader(site, cfg, cfg.stocked_enzymes) + site
    if right_enzyme:
        site = revcomp(CATALOG[right_enzyme].site)
        r_tail = choose_leader(site, cfg, cfg.stocked_enzymes) + site

    fwd = Primer(f"{fragment_name}-F", f_tail + f_anneal, f_anneal, f_tail, "F",
                 fragment_name, start, end)
    rev = Primer(f"{fragment_name}-R", r_tail + r_anneal, r_anneal, r_tail, "R",
                 fragment_name, start, end)
    for pr, at in ((fwd, start), (rev, end)):
        pr.tpl_seq, pr.tpl_circular, pr.anneal_at = (
            template, circular_template, at)
    fwd.warnings = validate_primer(fwd, template, cfg, circular_template)
    rev.warnings = validate_primer(rev, template, cfg, circular_template)
    limit = getattr(cfg, "primer_pair_max_diff", 5.0)
    if tm_gap == tm_gap and tm_gap > limit:
        for pr in (fwd, rev):
            pr.warnings.append(
                f"this pair's Tms differ by {tm_gap:.1f} C, more than the "
                f"{limit:.0f} C two primers sharing one annealing step "
                f"should; no combination of lengths here does better, so "
                f"run a gradient over the lower Tm")
    return fwd, rev


# --------------------------------------------------------------------------
# Junction analysis
# --------------------------------------------------------------------------

TRANSLATED = {"CDS", "trafficking_motif", "linker", "p2A", "tag"}

#: Homology each kit needs, as (minimum, maximum) base pairs.
#:
#: NEBuilder HiFi's higher-fidelity enzyme mix works with 15-25 bp, where
#: standard Gibson typically wants 20-40 bp. Shorter homology means shorter
#: oligos, so HiFi is the easier of the two on primer design, not the
#: harder. Both want the arm's Tm at or above 48 C, since the reaction runs
#: at 50 C.
#: https://www.neb.com/en-us/tools-and-resources/usage-guidelines/optimization-tips-for-nebuilder-hifi-dna-assembly-and-neb-gibson-assembly
ASSEMBLY_HOMOLOGY = {
    "hifi": (15, 25),
    "gibson": (20, 40),
}


def homology_range(cfg: Config) -> tuple:
    kit = (cfg.preferred_assembly_kit or "").lower()
    key = "gibson" if "gibson" in kit and "hifi" not in kit else "hifi"
    low, high = ASSEMBLY_HOMOLOGY[key]
    # The kit decides the ceiling. `overlap_length` is only a floor-raiser
    # for someone who wants longer arms than the kit strictly needs.
    return max(low, cfg.overlap_min_length), max(high, cfg.overlap_length)


#: How far an arm may be grown to get past a repeat. Beyond this the oligo
#: stops being orderable as a stock primer.
UNIQUE_ARM_MAX = 60


def _occurrences(arm: str, host: str) -> int:
    """How many times an arm appears in a circular molecule, either strand."""
    if not arm or not host:
        return 1
    arm = arm.upper()
    host = host.upper()
    wrapped = host + host[:len(arm) - 1]
    return wrapped.count(arm) + wrapped.count(revcomp(arm))


def tail_match(tail: str, template: str, at: int, direction: str,
               circular: bool = False) -> int:
    """How far the tail goes on matching the template past the junction.

    `at` is where the annealing region begins on the template for a
    forward primer, or where it ends for a reverse one. Walking outward
    from there, the arm is compared with the template on the primer's own
    strand and the run of agreement counted.

    This is not the annealing site having been moved. Nobody put these
    bases there: the arm is the neighbouring fragment's end and it simply
    happens to start the way the template does - which is common, because
    a junction usually falls on a Kozak or a start codon and both sides
    of it carry GCCACC or ATG. On the iGluSnFR4f forward primer it is two
    bases, and they are worth 4.0 C.
    """
    if not tail or not template:
        return 0
    k = len(tail)
    if direction == "F":
        before = slice_region(template, at - k, at, circular) \
            if (circular or at - k >= 0) else template[:at]
    else:
        after = slice_region(template, at, at + k, circular) \
            if (circular or at + k <= len(template)) else template[at:]
        before = revcomp(after)
    n = 0
    while (n < len(tail) and n < len(before)
           and tail[len(tail) - 1 - n].upper()
           == before[len(before) - 1 - n].upper()):
        n += 1
    return n


def pick_overlap(neighbour: str, cfg: Config, from_end: bool = True,
                 unique_in: str = "") -> str:
    """Take the shortest homology arm that still anneals at the assembly
    temperature, and that occurs only once.

    Using the maximum every time makes every oligo longer than it needs to
    be. The arm is grown from the kit's minimum until its Tm reaches the
    threshold, and only then does it stop.

    Length is not the only requirement. An arm has to say *where* a fragment
    goes, which it cannot do if the same bases occur twice. A Cre-lox
    construct carries each lox site twice by design, so an arm taken from
    one of them matches both, and the payload assembles either way round -
    in the simulation, and in the tube. Where that happens the arm keeps
    growing until it reaches past the repeat into sequence that occurs once.
    """
    low, high = homology_range(cfg)
    if not neighbour:
        return ""

    def arm_of(n: int) -> str:
        return neighbour[-n:] if from_end else neighbour[:n]

    best = ""
    limit = min(high, len(neighbour))
    for length in range(low, limit + 1):
        arm = arm_of(length)
        best = arm
        t = tm(arm, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
        if t == t and t >= cfg.overlap_tm_min:
            if _occurrences(arm, unique_in) <= 1:
                return arm
            break

    if unique_in and best and _occurrences(best, unique_in) > 1:
        for length in range(limit + 1, min(UNIQUE_ARM_MAX,
                                           len(neighbour)) + 1):
            arm = arm_of(length)
            if _occurrences(arm, unique_in) <= 1:
                return arm
    return best


def analyse_junctions(layout: list[dict], seq: str, method: str,
                      cfg: Config, joins: set | None = None) -> list[Junction]:
    """Describe every module boundary, and which of them are real joins.

    Two different things used to be called a junction. A boundary
    between two annotated elements matters for the reading frame - a
    CDS that is not a multiple of three puts the next one out of frame
    whether or not anything is joined there. A *join* is where two
    amplicons meet, and only there is there an overlap arm whose GC and
    Tm decide whether the assembly works.

    Reporting all of them as HiFi joins with a 25 bp overlap gave
    seventeen rows for a three-fragment assembly, and - worse - warned
    that "the overlap arm is only 24% GC" at the boundary between the
    AmpR promoter and AmpR, which sit in the middle of one PCR product
    with no overlap arm between them at all. That is a redesign of a
    primer that does not exist.

    ``joins`` is the set of construct positions where one piece ends and
    the next begins. Without it every boundary is treated as a join,
    which is the old behaviour and is right for a caller that has not
    worked out its pieces yet.
    """
    junctions: list[Junction] = []
    for i in range(len(layout) - 1):
        left, right = layout[i], layout[i + 1]
        pos = left["end"]
        window = seq[max(pos - 15, 0): pos + 15]
        is_join = joins is None or pos in joins
        j = Junction(
            index=i + 1, left_name=left["name"], right_name=right["name"],
            left_type=left["type"], right_type=right["type"],
            method=method if is_join else "",
            sequence=window, is_join=is_join,
            overlap_length=(cfg.overlap_length
                            if is_join and method in ("gibson", "hifi")
                            else 0),
        )


        both_coding = left["type"] in TRANSLATED and right["type"] in TRANSLATED
        if both_coding:
            left_len = left["end"] - left["start"]
            if left_len % 3 != 0:
                j.frame_shift = left_len % 3
                j.warnings.append(
                    f"{left['name']} is {left_len} bp, so {right['name']} starts "
                    f"{j.frame_shift} base(s) out of frame")
            codon = seq[pos - 3:pos]
            if translate(codon) == "*":
                j.warnings.append(
                    f"{left['name']} ends in a stop codon, which will terminate "
                    f"translation before {right['name']}")

        # GC-rich junctions are the ones that fail in HiFi assembly -
        # and only where there is an arm to be GC-rich.
        if is_join and method in ("gibson", "hifi"):
            arm = seq[max(pos - cfg.overlap_length, 0):pos]
            g = gc_fraction(arm)
            if g > 0.75:
                j.warnings.append(f"overlap arm is {g:.0%} GC; HiFi joining may be poor")
            elif g < 0.25:
                j.warnings.append(f"overlap arm is only {g:.0%} GC; overlap Tm will be low")
            t_arm = tm(arm, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
            if t_arm == t_arm and t_arm < cfg.overlap_tm_min:
                j.warnings.append(f"overlap arm Tm {t_arm:.1f} C is low")

        junctions.append(j)
    return junctions


def check_overlap_uniqueness(layout: list[dict], seq: str, cfg: Config) -> list[str]:
    """Overlaps must be unique, or fragments will assemble in the wrong order."""
    warnings: list[str] = []
    arms = []
    for i in range(len(layout) - 1):
        pos = layout[i]["end"]
        arms.append((i, seq[max(pos - cfg.overlap_length, 0):pos]))
    for i, (idx_a, a) in enumerate(arms):
        for idx_b, b in arms[i + 1:]:
            if not a or not b:
                continue
            if a == b:
                warnings.append(
                    f"junctions {idx_a + 1} and {idx_b + 1} have identical "
                    f"overlap arms; fragments could assemble in the wrong order")
            elif len(a) >= 15 and (a[-15:] == b[-15:]):
                warnings.append(
                    f"junctions {idx_a + 1} and {idx_b + 1} share their last "
                    f"15 bp; consider shifting one boundary")
    return warnings
