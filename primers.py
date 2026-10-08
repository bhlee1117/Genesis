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
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------------------------------------
# Annealing-region selection
# --------------------------------------------------------------------------

def pick_anneal(template: str, start: int, direction: str, cfg: Config,
                circular: bool = False) -> str:
    """Grow an annealing region from `start` until it reaches the target Tm.

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
                       cfg: Config, circular_template: bool = False) -> tuple[Primer, Primer]:
    """Design the primer pair that amplifies one fragment with HiFi overlaps.

    The forward primer carries the 3' end of the *upstream* fragment as a
    tail; the reverse primer carries the reverse complement of the 5' start
    of the *downstream* fragment.
    """
    ov = cfg.overlap_length
    f_anneal = pick_anneal(template, start, "F", cfg, circular_template)
    r_anneal = pick_anneal(template, end, "R", cfg, circular_template)

    f_tail = upstream_seq[-ov:] if upstream_seq else ""
    r_tail = revcomp(downstream_seq[:ov]) if downstream_seq else ""

    fwd = Primer(
        name=f"{fragment_name}-F", sequence=f_tail + f_anneal, anneal=f_anneal,
        tail=f_tail, direction="F", template=fragment_name,
        target_start=start, target_end=end)
    rev = Primer(
        name=f"{fragment_name}-R", sequence=r_tail + r_anneal, anneal=r_anneal,
        tail=r_tail, direction="R", template=fragment_name,
        target_start=start, target_end=end)

    fwd.warnings = validate_primer(fwd, template, cfg, circular_template)
    rev.warnings = validate_primer(rev, template, cfg, circular_template)

    for p, tail in ((fwd, f_tail), (rev, r_tail)):
        if tail:
            t_ov = tm(tail, cfg.primer_conc_nm, cfg.na_mm, cfg.mg_mm)
            if t_ov == t_ov and t_ov < cfg.overlap_tm_min:
                p.warnings.append(
                    f"overlap tail Tm {t_ov:.1f} C is below the "
                    f"{cfg.overlap_tm_min:.0f} C needed for efficient HiFi joining")
    return fwd, rev


def design_re_pair(fragment_name: str, template: str, start: int, end: int,
                   left_enzyme: str, right_enzyme: str, cfg: Config,
                   circular_template: bool = False) -> tuple[Primer, Primer]:
    """Primers that add restriction sites plus a 5' clamp for efficient cutting."""
    from .enzymes import CATALOG
    clamp = "TTTT"
    f_anneal = pick_anneal(template, start, "F", cfg, circular_template)
    r_anneal = pick_anneal(template, end, "R", cfg, circular_template)
    f_tail = clamp + CATALOG[left_enzyme].site if left_enzyme else ""
    r_tail = clamp + revcomp(CATALOG[right_enzyme].site) if right_enzyme else ""

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


def analyse_junctions(layout: list[dict], seq: str, method: str,
                      cfg: Config) -> list[Junction]:
    """Describe every module boundary, including its effect on the frame."""
    junctions: list[Junction] = []
    for i in range(len(layout) - 1):
        left, right = layout[i], layout[i + 1]
        pos = left["end"]
        window = seq[max(pos - 15, 0): pos + 15]
        j = Junction(
            index=i + 1, left_name=left["name"], right_name=right["name"],
            left_type=left["type"], right_type=right["type"], method=method,
            sequence=window,
            overlap_length=cfg.overlap_length if method in ("gibson", "hifi") else 0,
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

        # GC-rich junctions are the ones that fail in HiFi assembly.
        if method in ("gibson", "hifi"):
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
