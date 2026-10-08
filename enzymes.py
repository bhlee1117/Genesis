"""Restriction enzyme catalogue and digest simulation (Stage 7 foundation).

Cut positions follow the usual convention: ``top`` and ``bottom`` are offsets
from the start of the recognition site on the top strand.

    EcoRI  G^AATTC   top=1 bottom=5  ->  5' overhang AATT
    KpnI   GGTAC^C   top=5 bottom=1  ->  3' overhang GTAC
    EcoRV  GAT^ATC   top=3 bottom=3  ->  blunt

``bottom - top`` is therefore positive for a 5' overhang, negative for a 3'
overhang and zero for a blunt cut. Type IIS enzymes simply have offsets that
fall outside the recognition sequence.
"""

from __future__ import annotations

from dataclasses import dataclass

from .seqtools import find_all, iupac_regex, revcomp, slice_region


@dataclass(frozen=True)
class Enzyme:
    name: str
    site: str              # recognition sequence, IUPAC, top strand
    top: int               # top-strand cut offset
    bottom: int            # bottom-strand cut offset
    dam_sensitive: bool = False
    dcm_sensitive: bool = False
    cpg_sensitive: bool = False
    star_activity: bool = False
    buffer: str = "rCutSmart"
    notes: str = ""

    @property
    def overhang_len(self) -> int:
        return self.bottom - self.top

    @property
    def end_type(self) -> str:
        if self.overhang_len > 0:
            return "5'"
        if self.overhang_len < 0:
            return "3'"
        return "blunt"

    @property
    def is_type_iis(self) -> bool:
        """Cuts outside its recognition sequence."""
        return self.top > len(self.site) or self.bottom > len(self.site)

    @property
    def palindromic(self) -> bool:
        return self.site == revcomp(self.site)


def _e(*args, **kwargs) -> Enzyme:
    return Enzyme(*args, **kwargs)


#: A working catalogue covering the lab's stocked set plus common extras.
CATALOG: dict[str, Enzyme] = {e.name: e for e in [
    # --- 5' overhang, 6-cutters ---------------------------------------
    _e("EcoRI", "GAATTC", 1, 5, star_activity=True, buffer="rCutSmart"),
    _e("BamHI", "GGATCC", 1, 5, star_activity=True),
    _e("HindIII", "AAGCTT", 1, 5),
    _e("XbaI", "TCTAGA", 1, 5, dam_sensitive=True,
       notes="blocked when overlapped by dam methylation (GATC)"),
    _e("XhoI", "CTCGAG", 1, 5),
    _e("SalI", "GTCGAC", 1, 5, star_activity=True),
    _e("NheI", "GCTAGC", 1, 5),
    _e("SpeI", "ACTAGT", 1, 5),
    _e("AgeI", "ACCGGT", 1, 5),
    _e("BglII", "AGATCT", 1, 5),
    _e("MluI", "ACGCGT", 1, 5),
    _e("BsrGI", "TGTACA", 1, 5),
    _e("NcoI", "CCATGG", 1, 5),
    _e("AflII", "CTTAAG", 1, 5),
    _e("AvrII", "CCTAGG", 1, 5),
    _e("BspEI", "TCCGGA", 1, 5, dam_sensitive=True),
    _e("EagI", "CGGCCG", 1, 5, cpg_sensitive=True),
    _e("BsiWI", "CGTACG", 1, 5),
    _e("XmaI", "CCCGGG", 1, 5),
    _e("AscI", "GGCGCGCC", 2, 6, cpg_sensitive=True),
    _e("NotI", "GCGGCCGC", 2, 6, cpg_sensitive=True),
    _e("ClaI", "ATCGAT", 2, 4, dam_sensitive=True),
    _e("NdeI", "CATATG", 2, 4),
    _e("BstBI", "TTCGAA", 2, 4),
    _e("FseI", "GGCCGGCC", 6, 2),
    # --- 3' overhang ---------------------------------------------------
    _e("KpnI", "GGTACC", 5, 1),
    _e("SacI", "GAGCTC", 5, 1),
    _e("PstI", "CTGCAG", 5, 1),
    _e("ApaI", "GGGCCC", 5, 1),
    _e("SphI", "GCATGC", 5, 1),
    _e("AatII", "GACGTC", 5, 1),
    _e("NsiI", "ATGCAT", 5, 1),
    _e("SbfI", "CCTGCAGG", 6, 2),
    _e("PacI", "TTAATTAA", 5, 3),
    _e("SacII", "CCGCGG", 4, 2, cpg_sensitive=True),
    _e("PvuI", "CGATCG", 4, 2),
    # --- blunt ---------------------------------------------------------
    _e("EcoRV", "GATATC", 3, 3),
    _e("SmaI", "CCCGGG", 3, 3),
    _e("HpaI", "GTTAAC", 3, 3),
    _e("StuI", "AGGCCT", 3, 3, dcm_sensitive=True),
    _e("NruI", "TCGCGA", 3, 3),
    _e("PvuII", "CAGCTG", 3, 3, star_activity=True),
    _e("ScaI", "AGTACT", 3, 3),
    _e("DraI", "TTTAAA", 3, 3),
    _e("SnaBI", "TACGTA", 3, 3),
    _e("PmeI", "GTTTAAAC", 4, 4),
    _e("SwaI", "ATTTAAAT", 4, 4),
    # --- Type IIS (cut outside the recognition site) --------------------
    _e("BsaI", "GGTCTC", 7, 11, notes="Golden Gate; 4 nt 5' overhang"),
    _e("BsmBI", "CGTCTC", 7, 11, notes="Golden Gate; 4 nt 5' overhang"),
    _e("BbsI", "GAAGAC", 8, 12),
    _e("SapI", "GCTCTTC", 8, 11),
    _e("BspQI", "GCTCTTC", 8, 11),
]}


# --------------------------------------------------------------------------
# Cut-site finding
# --------------------------------------------------------------------------

@dataclass
class CutSite:
    enzyme: str
    site_start: int        # 0-based start of the recognition sequence
    strand: int            # which strand the recognition site was read on
    top_cut: int           # absolute position of the top-strand nick
    bottom_cut: int        # absolute position of the bottom-strand nick

    @property
    def position(self) -> int:
        return self.top_cut


def find_sites(seq: str, enzyme: Enzyme, circular: bool = False) -> list[CutSite]:
    """Locate every cut this enzyme makes in `seq`."""
    n = len(seq)
    sites: list[CutSite] = []
    if n == 0:
        return sites

    for pos in find_all(seq, enzyme.site, circular=circular):
        top = pos + enzyme.top
        bottom = pos + enzyme.bottom
        if not circular and (top > n or bottom > n):
            continue
        sites.append(CutSite(enzyme.name, pos, 1, top % n if circular else top,
                             bottom % n if circular else bottom))

    if not enzyme.palindromic:
        rc_site = revcomp(enzyme.site)
        for pos in find_all(seq, rc_site, circular=circular):
            # Mirror the offsets onto the bottom strand.
            end = pos + len(enzyme.site)
            top = end - enzyme.bottom
            bottom = end - enzyme.top
            if not circular and (top < 0 or bottom < 0):
                continue
            sites.append(CutSite(enzyme.name, pos, -1,
                                 top % n if circular else top,
                                 bottom % n if circular else bottom))

    # A palindromic site read on both strands yields the same cut twice.
    unique: dict[tuple, CutSite] = {}
    for s in sites:
        unique.setdefault((s.top_cut, s.bottom_cut, s.enzyme), s)
    return sorted(unique.values(), key=lambda s: s.top_cut)


def count_sites(seq: str, enzyme_name: str, circular: bool = False) -> int:
    enz = CATALOG.get(enzyme_name)
    if enz is None:
        return 0
    return len(find_sites(seq, enz, circular))


# --------------------------------------------------------------------------
# Digest simulation
# --------------------------------------------------------------------------

@dataclass
class Fragment:
    seq: str               # top-strand sequence spanning the full fragment
    start: int
    end: int
    left_enzyme: str | None
    right_enzyme: str | None
    left_overhang: str = ""     # single-stranded extension, 5'->3'
    right_overhang: str = ""
    left_end_type: str = "blunt"
    right_end_type: str = "blunt"

    @property
    def length(self) -> int:
        return len(self.seq)


def digest(seq: str, enzyme_names: list[str], circular: bool = False) -> list[Fragment]:
    """Simulate a single or multiple digest and return the fragments.

    Fragment boundaries are taken at the top-strand nick; the resulting
    single-stranded overhang is recorded so that end compatibility can be
    checked later.
    """
    n = len(seq)
    cuts: list[CutSite] = []
    for name in enzyme_names:
        enz = CATALOG.get(name)
        if enz is None:
            raise KeyError(f"unknown enzyme: {name}")
        cuts.extend(find_sites(seq, enz, circular))

    cuts.sort(key=lambda c: c.top_cut)
    if not cuts:
        return [Fragment(seq, 0, n, None, None)]

    def overhang_for(cut: CutSite) -> tuple[str, str]:
        enz = CATALOG[cut.enzyme]
        olen = enz.overhang_len
        if olen == 0:
            return "", "blunt"
        lo, hi = sorted((cut.top_cut, cut.bottom_cut))
        oh = slice_region(seq, lo, hi, circular=circular)
        return oh, enz.end_type

    frags: list[Fragment] = []
    if circular:
        for i, cut in enumerate(cuts):
            nxt = cuts[(i + 1) % len(cuts)]
            start, end = cut.top_cut, nxt.top_cut
            sub = slice_region(seq, start, end if end != start else start + n,
                               circular=True)
            loh, lt = overhang_for(cut)
            roh, rt = overhang_for(nxt)
            frags.append(Fragment(sub, start, end, cut.enzyme, nxt.enzyme,
                                  loh, roh, lt, rt))
    else:
        boundaries = [0] + [c.top_cut for c in cuts] + [n]
        for i in range(len(boundaries) - 1):
            start, end = boundaries[i], boundaries[i + 1]
            if end <= start:
                continue
            left_cut = cuts[i - 1] if i > 0 else None
            right_cut = cuts[i] if i < len(cuts) else None
            loh, lt = overhang_for(left_cut) if left_cut else ("", "blunt")
            roh, rt = overhang_for(right_cut) if right_cut else ("", "blunt")
            frags.append(Fragment(
                seq[start:end], start, end,
                left_cut.enzyme if left_cut else None,
                right_cut.enzyme if right_cut else None,
                loh, roh, lt, rt))
    return frags


def ends_compatible(a_overhang: str, a_type: str,
                    b_overhang: str, b_type: str) -> bool:
    """Can these two ends be ligated together?

    Two 5' overhangs anneal when one is the reverse complement of the other.
    Blunt ends always join (inefficiently). A 5' and a 3' overhang never do.
    """
    if a_type == "blunt" and b_type == "blunt":
        return True
    if a_type != b_type:
        return False
    if not a_overhang or not b_overhang:
        return False
    if len(a_overhang) != len(b_overhang):
        return False
    return a_overhang.upper() == revcomp(b_overhang).upper()


def methylation_warnings(seq: str, enzyme: Enzyme, circular: bool = False) -> list[str]:
    """Flag dam/dcm overlap, which silently blocks digestion of plasmid DNA
    grown in standard E. coli."""
    warnings: list[str] = []
    if enzyme.dam_sensitive:
        for site in find_sites(seq, enzyme, circular):
            context = slice_region(seq, max(site.site_start - 4, 0),
                                   site.site_start + len(enzyme.site) + 4,
                                   circular=circular)
            if "GATC" in context.upper():
                warnings.append(
                    f"{enzyme.name} site at {site.site_start} overlaps a dam "
                    f"site (GATC); digestion will be blocked in dam+ strains")
    if enzyme.dcm_sensitive:
        rx = iupac_regex("CCWGG")
        for site in find_sites(seq, enzyme, circular):
            context = slice_region(seq, max(site.site_start - 4, 0),
                                   site.site_start + len(enzyme.site) + 4,
                                   circular=circular)
            if rx.search(context):
                warnings.append(
                    f"{enzyme.name} site at {site.site_start} overlaps a dcm "
                    f"site (CCWGG); use a dcm- strain")
    return warnings
