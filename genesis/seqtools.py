"""Pure-sequence utilities: complement, translation, Tm, searching, repeats.

Everything here is stdlib-only and side-effect free so it can be unit tested
without a database or any lab files.
"""

from __future__ import annotations

import math
import re
from typing import Iterator

# --------------------------------------------------------------------------
# Alphabet handling
# --------------------------------------------------------------------------

_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def clean_seq(seq: str) -> str:
    """Strip whitespace, digits and GenBank line numbering; uppercase."""
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def iupac_regex(pattern: str) -> re.Pattern:
    """Compile an IUPAC pattern into a regex with character classes."""
    out = []
    for ch in pattern.upper():
        expand = IUPAC.get(ch, ch)
        out.append(expand if len(expand) == 1 else f"[{expand}]")
    return re.compile("".join(out), re.IGNORECASE)


# --------------------------------------------------------------------------
# Searching (linear and circular)
# --------------------------------------------------------------------------

_PLAIN_DNA = re.compile(r"^[ACGTacgt]+$")


def find_all(haystack: str, needle: str, circular: bool = False) -> list[int]:
    """All 0-based start offsets of `needle`, allowing IUPAC codes.

    Overlapping matches are reported. For circular molecules, matches that
    wrap the origin are reported at their start offset in the original
    coordinate frame.
    """
    if not needle or not haystack:
        return []
    n = len(haystack)
    target = haystack + (haystack[: len(needle) - 1] if circular else "")
    hits: list[int] = []

    if _PLAIN_DNA.match(needle):
        # Fast path: plain C-level substring search, stepping by 1 to keep
        # overlapping hits.
        needle_u = needle.upper()
        target_u = target.upper()
        i = target_u.find(needle_u)
        while i != -1:
            if i >= n:
                break
            hits.append(i)
            i = target_u.find(needle_u, i + 1)
        return hits

    # IUPAC path: a lookahead keeps finditer overlapping.
    rx = re.compile(f"(?=({iupac_regex(needle).pattern}))", re.IGNORECASE)
    for m in rx.finditer(target):
        if m.start() >= n:
            break
        hits.append(m.start())
    return hits


def _slice_circular(seq: str, start: int, end: int) -> str:
    """Slice allowing end > len(seq) to wrap around the origin."""
    n = len(seq)
    if n == 0:
        return ""
    if end <= n:
        return seq[start:end]
    return seq[start:] + seq[: end - n]


def slice_region(seq: str, start: int, end: int, circular: bool = False) -> str:
    """Extract [start, end) honouring wrap-around on circular molecules.

    Callers write a wrapped span either way round: with `end` past the last
    base, which is how the occurrence table records one, or reduced so that
    `end` lands at or before `start`. Both mean the same arc.

    How long that arc is must be worked out before either coordinate is
    reduced. Reducing `start` alone and leaving `end` where it was turned
    an 18 base window at the origin into a whole turn of the plasmid plus
    18, which is how primer design came to anneal a 1,174 nt oligo to a
    1,156 bp template.
    """
    if circular:
        n = len(seq)
        if n == 0:
            return ""
        span = end - start
        if span <= 0:
            span += n
        start %= n
        return _slice_circular(seq, start, start + span)
    return seq[start:end]


# --------------------------------------------------------------------------
# Translation
# --------------------------------------------------------------------------

CODON_TABLE = {}
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _aa in enumerate(_AAS):
    CODON_TABLE["".join((_BASES[_i // 16], _BASES[(_i // 4) % 4], _BASES[_i % 4]))] = _aa

START_CODONS = {"ATG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate(seq: str, frame: int = 0, to_stop: bool = False) -> str:
    """Translate DNA to protein. Ambiguous codons become 'X'."""
    seq = clean_seq(seq)[frame:]
    out = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        if aa == "*" and to_stop:
            break
        out.append(aa)
    return "".join(out)


def find_orfs(seq: str, min_aa: int = 50, circular: bool = False) -> list[dict]:
    """Find ORFs (ATG..stop) on both strands. Coordinates are 0-based."""
    orfs = []
    for strand, s in ((1, seq), (-1, revcomp(seq))):
        search = s + (s[:300] if circular else "")
        for m in re.finditer(r"ATG", search):
            start = m.start()
            if start >= len(s):
                continue
            for i in range(start, len(search) - 2, 3):
                if search[i:i + 3] in STOP_CODONS:
                    length_aa = (i - start) // 3
                    if length_aa >= min_aa:
                        orfs.append({
                            "strand": strand,
                            "start": start if strand == 1 else len(s) - (i + 3),
                            "end": (i + 3) if strand == 1 else len(s) - start,
                            "aa_length": length_aa,
                            "protein": translate(search[start:i]),
                        })
                    break
    return sorted(orfs, key=lambda o: -o["aa_length"])


# --------------------------------------------------------------------------
# Thermodynamics - SantaLucia (1998) unified nearest-neighbour parameters
# --------------------------------------------------------------------------

# (dH kcal/mol, dS cal/mol/K)
_NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
_R = 1.987  # cal/(mol*K)


def tm(seq: str, primer_nm: float = 500.0, na_mm: float = 50.0,
       mg_mm: float = 2.0, dntp_mm: float = 0.2) -> float:
    """Nearest-neighbour melting temperature in degrees C.

    Uses the SantaLucia unified parameters with a von Ahsen divalent-cation
    correction, which is the regime most high-fidelity PCR kits run in.
    """
    seq = clean_seq(seq)
    if len(seq) < 2 or any(c not in "ACGT" for c in seq):
        return float("nan")

    dh, ds = 0.0, 0.0
    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        if pair not in _NN:
            return float("nan")
        d1, d2 = _NN[pair]
        dh += d1
        ds += d2

    # Helix initiation penalties, by terminal base pair identity.
    for end in (seq[0], seq[-1]):
        if end in "GC":
            dh += 0.1
            ds += -2.8
        else:
            dh += 2.3
            ds += 4.1

    # Effective monovalent cation concentration (von Ahsen et al. 2001):
    # free Mg2+ dominates, and dNTPs chelate it.
    free_mg = max(mg_mm - dntp_mm, 0.0)
    na_eq = (na_mm + 120.0 * math.sqrt(free_mg)) / 1000.0
    if na_eq <= 0:
        na_eq = 1e-6
    ds += 0.368 * (len(seq) - 1) * math.log(na_eq)

    ct = primer_nm * 1e-9
    try:
        t = (dh * 1000.0) / (ds + _R * math.log(ct / 4.0)) - 273.15
    except (ValueError, ZeroDivisionError):
        return float("nan")
    return round(t, 1)


def hairpin_score(seq: str, min_stem: int = 4, min_loop: int = 3) -> dict | None:
    """Longest self-complementary stem that could form a hairpin."""
    seq = clean_seq(seq)
    best = None
    n = len(seq)
    for i in range(n):
        for j in range(i + min_stem + min_loop, n + 1):
            for stem in range(min_stem, (j - i - min_loop) // 2 + 1):
                left = seq[i:i + stem]
                right = seq[j - stem:j]
                if left == revcomp(right):
                    if best is None or stem > best["stem"]:
                        best = {"stem": stem, "start": i, "end": j,
                                "loop": (j - stem) - (i + stem)}
    return best


def self_dimer_3prime(seq: str, window: int = 6) -> int:
    """Length of 3'-end self-complementarity (the dimer type that matters)."""
    seq = clean_seq(seq)
    tail = seq[-window:]
    rc = revcomp(seq)
    best = 0
    for k in range(3, len(tail) + 1):
        if tail[-k:] in rc:
            best = k
    return best


def max_homopolymer(seq: str) -> int:
    runs = re.findall(r"(A+|C+|G+|T+)", clean_seq(seq))
    return max((len(r) for r in runs), default=0)


def longest_repeat(seq: str, min_len: int = 20, circular: bool = False) -> list[dict]:
    """Maximal direct repeats >= min_len, which jeopardise Gibson/HiFi specificity.

    A sliding k-mer scan reports one hit per offset, so a single 200 bp
    repeat would otherwise surface as ~180 near-identical warnings. Hits
    sharing an offset difference and lying at consecutive positions belong
    to the same repeat, so they are merged and reported once at full length.
    """
    seq = clean_seq(seq)
    seen: dict[str, int] = {}
    raw: list[tuple[int, int]] = []
    for i in range(len(seq) - min_len + 1):
        kmer = seq[i:i + min_len]
        if kmer in seen:
            raw.append((seen[kmer], i))
        else:
            seen[kmer] = i

    by_offset: dict[int, list[tuple[int, int]]] = {}
    for first, second in raw:
        by_offset.setdefault(second - first, []).append((first, second))

    merged: list[dict] = []
    for delta, pairs in by_offset.items():
        pairs.sort()
        run_start = pairs[0]
        prev = pairs[0]
        for pair in pairs[1:] + [(None, None)]:
            if pair[0] is not None and pair[0] == prev[0] + 1:
                prev = pair
                continue
            length = (prev[0] - run_start[0]) + min_len
            merged.append({
                "seq": seq[run_start[0]:run_start[0] + length],
                "first": run_start[0], "second": run_start[1], "length": length,
            })
            if pair[0] is not None:
                run_start = prev = pair
    return sorted(merged, key=lambda r: -r["length"])


def hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def iter_kmers(seq: str, k: int, circular: bool = False) -> Iterator[tuple[int, str]]:
    n = len(seq)
    limit = n if circular else n - k + 1
    for i in range(max(limit, 0)):
        yield i, _slice_circular(seq, i, i + k) if circular else seq[i:i + k]
