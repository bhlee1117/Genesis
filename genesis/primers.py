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
    #: The part of `tail` that is template the annealing site stepped over,
    #: sitting at its 3' end. Anything spliced into the tail later has to go
    #: in FRONT of these, or the amplicon reads them out of order. Declared
    #: last on purpose: other call sites build a Primer positionally, and
    #: inserting a field mid-list silently shifted "F" into this one and
    #: every argument after it.
    shift_tail: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_dict(self, cfg: Config | None = None) -> dict:
        cfg = cfg or Config()
        d = asdict(self)
        d["length"] = self.length
        d["anneal_tm"] = tm(self.anneal, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
        d["full_tm"] = tm(self.sequence, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
        d["gc"] = round(gc_fraction(self.anneal), 3)
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

def anneal_penalty(region: str, cfg: Config, shift: int = 0,
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

    # Every base of shift is a base of oligo somebody pays for.
    score += 0.25 * shift

    if unique_in:
        from .seqtools import find_all
        n = len(find_all(unique_in, region, circular)) + \
            len(find_all(unique_in, revcomp(region), circular))
        if n > 1:
            score += 25.0 * (n - 1)
    return score


def pick_anneal_site(template: str, start: int, direction: str, cfg: Config,
                     circular: bool = False, max_shift: int = 0,
                     unique_in: str = "") -> tuple[str, int]:
    """Choose where the primer binds, and how far in from the junction.

    `pick_anneal` can only change how long the primer is, not where it
    sits, so a junction falling in an AT-rich stretch produced a primer
    that was reported as too cold and used regardless. The site may now
    move *into* the fragment; the bases it steps over are added to the
    primer's tail, which leaves the product identical and makes the oligo
    longer by exactly the shift.

    It never moves outward, because that would amplify template bases the
    construct does not call for.

    Returns the annealing region and the shift it took.
    """
    best, best_score = ("", 0), float("inf")
    for shift in range(0, max(0, max_shift) + 1):
        for length in range(cfg.primer_min_len, cfg.primer_max_len + 1):
            if direction == "F":
                lo = start + shift
                region = slice_region(template, lo, lo + length, circular)
            else:
                hi = start - shift
                lo = hi - length
                if lo < 0 and not circular:
                    break
                region = revcomp(slice_region(template, lo, hi, circular))
            if len(region) < length:
                break
            score = anneal_penalty(region, cfg, shift, unique_in, circular)
            if score < best_score:
                best_score, best = score, (region, shift)
    if best[0]:
        return best
    return pick_anneal(template, start, direction, cfg, circular), 0


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

    t = tm(anneal, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
    if t != t:
        n_ambig = sum(1 for ch in anneal if ch not in "ACGT")
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
                       unique_in: str = "") -> tuple[Primer, Primer]:
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
    shift = max(0, getattr(cfg, "primer_max_shift", 0))
    f_anneal, f_shift = pick_anneal_site(template, start, "F", cfg,
                                         circular_template, shift,
                                         unique_in=template)
    r_anneal, r_shift = pick_anneal_site(template, end, "R", cfg,
                                         circular_template, shift,
                                         unique_in=template)

    f_tail = pick_overlap(upstream_seq, cfg, from_end=True,
                          unique_in=unique_in)
    r_tail = (revcomp(pick_overlap(downstream_seq, cfg, from_end=False,
                                   unique_in=unique_in))
              if carry_downstream else "")

    # Whatever the annealing site stepped over has to be put back, or the
    # amplicon is short by exactly that much and the junction it was
    # designed for no longer exists. It goes on the tail, between the
    # overlap and the binding region, which is where the polymerase will
    # copy it from.
    f_overlap, r_overlap = len(f_tail), len(r_tail)
    if f_shift:
        f_tail = f_tail + slice_region(template, start, start + f_shift,
                                       circular_template)
    if r_shift:
        r_tail = r_tail + revcomp(slice_region(template, end - r_shift, end,
                                               circular_template))

    fwd = Primer(
        name=f"{fragment_name}-F", sequence=f_tail + f_anneal, anneal=f_anneal,
        tail=f_tail, shift_tail=f_tail[f_overlap:], direction="F",
        template=fragment_name,
        target_start=start, target_end=end)
    rev = Primer(
        name=f"{fragment_name}-R", sequence=r_tail + r_anneal, anneal=r_anneal,
        tail=r_tail, shift_tail=r_tail[r_overlap:], direction="R",
        template=fragment_name,
        target_start=start, target_end=end)

    fwd.warnings = validate_primer(fwd, template, cfg, circular_template)
    rev.warnings = validate_primer(rev, template, cfg, circular_template)

    # Say when the site was moved. A primer that sits 9 bases inside the
    # fragment rather than at its edge is a deliberate choice, and the
    # person checking it against their own map should not have to work
    # out why the numbers do not line up with the junction.
    for p, moved in ((fwd, f_shift), (rev, r_shift)):
        if moved:
            p.warnings.append(
                f"annealing site moved {moved} bp into the fragment to reach "
                f"better bases; those {moved} bp are carried on the tail, so "
                f"the product is unchanged")

    # Only the homology arm is an overlap. The bases added to reach a
    # better annealing site are template, not homology, and measuring
    # their melting temperature as if they were a HiFi arm reported a
    # 3 nt tail as -118 C.
    for p, tail in ((fwd, f_tail[:f_overlap]), (rev, r_tail[:r_overlap])):
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
    f_anneal = pick_anneal(template, start, "F", cfg, circular_template)
    r_anneal = pick_anneal(template, end, "R", cfg, circular_template)
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
    fwd.warnings = validate_primer(fwd, template, cfg, circular_template)
    rev.warnings = validate_primer(rev, template, cfg, circular_template)
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
