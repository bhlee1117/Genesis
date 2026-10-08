"""Stage 3: controlled feature vocabulary, plus a motif-based annotator.

Two complementary strategies are used, because they fail in different ways:

* **Label rules** map whatever the depositor typed ("Amp resistance",
  "CMV promotor") onto a controlled ``feature_type``. Cheap, but only works
  on files that were annotated in the first place.
* **Motif search** finds elements by sequence. Recombinase sites, ITRs and
  polyA signals are matched on DNA because they are near-invariant; peptide
  motifs (2A, ER export, trafficking) are matched on *protein* because codon
  usage differs between constructs while the peptide does not.

The original imported label is always preserved alongside the standard one.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .seqtools import revcomp, translate, find_all

# --------------------------------------------------------------------------
# Controlled vocabulary
# --------------------------------------------------------------------------

FEATURE_TYPES = [
    "backbone", "promoter", "enhancer", "CDS", "trafficking_motif", "linker",
    "p2A", "IRES", "Cre-lox", "lox_site", "AAV_ITR", "resistance", "origin",
    "tag", "polyA", "WPRE", "restriction_site", "primer_binding_site",
    "terminator", "intron", "unknown",
]

#: "DIO" was the old name for the family. Specs and parts saved with it keep
#: working; it now means Cre-lox in its DIO arrangement.
FEATURE_TYPE_ALIASES = {"DIO": "Cre-lox", "FLEX": "Cre-lox", "DO": "Cre-lox",
                        "FAS": "Cre-lox", "Cre_lox": "Cre-lox",
                        "crelox": "Cre-lox"}


def canonical_feature_type(value: str) -> str:
    """Map an old or informal type name onto the current vocabulary."""
    if not value:
        return ""
    if value in FEATURE_TYPES:
        return value
    return FEATURE_TYPE_ALIASES.get(value, FEATURE_TYPE_ALIASES.get(
        value.replace("-", "_"), value))

#: Feature types that must not be disrupted by a cut or a PCR boundary.
PROTECTED_TYPES = {"AAV_ITR", "lox_site", "Cre-lox", "origin", "resistance"}

#: Shortest length at which an annotation of each type is a believable part.
#: Lab files contain stray annotations such as an 8 bp feature labelled with
#: a fluorophore name; promoting those to reusable parts makes them look like
#: excellent candidates, because a short sequence matches many plasmids.
MIN_PART_LENGTH = {
    "CDS": 150, "promoter": 50, "enhancer": 50, "IRES": 300,
    "resistance": 300, "origin": 250, "AAV_ITR": 100, "DIO": 100,
    "WPRE": 100, "polyA": 30, "lox_site": 30, "p2A": 45,
    "trafficking_motif": 12, "linker": 9, "tag": 12, "backbone": 1000,
}
DEFAULT_MIN_PART_LENGTH = 10


def min_part_length(feature_type: str) -> int:
    return MIN_PART_LENGTH.get(feature_type, DEFAULT_MIN_PART_LENGTH)

# --------------------------------------------------------------------------
# Label-based rules.  Order matters: first match wins.
# --------------------------------------------------------------------------

_LABEL_RULES: list[tuple[str, str]] = [
    (r"\bITR\b|inverted terminal repeat", "AAV_ITR"),
    # loxFAS is a lox variant, and this lab writes it "FAS SITE".  It belongs
    # with the other sites, because that is what the recombination logic
    # looks for.  Bare "FAS" is not enough: Fas is also a gene name.
    (r"\blox(P|2272|N|511|66|71|FAS)?\b|\bFRT\b|\brox\b"
     r"|\bFAS[- _]?site\b", "lox_site"),
    (r"\bDIO\b|\bDO\b|double[- ]?floxed|\bFLEX\b|\bCreOn\b|\bCreOff\b",
     "Cre-lox"),
    (r"\bWPRE\b|woodchuck", "WPRE"),
    (r"poly[- ]?A|polyadenylation|\bbGH\b|\bSV40 ?pA\b", "polyA"),
    (r"\bIRES\b", "IRES"),
    (r"\b[PTEF]2A\b|\b2A (peptide|self)|self[- ]?cleav", "p2A"),
    (r"promot|\bCMV\b|\bEF1|\bhSyn\b|\bCaMKII|\bCAG\b|\bU6\b|\bTRE\b|\bPGK\b|\bSV40 ?prom", "promoter"),
    (r"enhancer", "enhancer"),
    (r"resistance|\bAmpR?\b|\bKanR?\b|\bbla\b|\bCmR\b|\bpuro|\bblast|\bneo\b|\bhygro", "resistance"),
    (r"\borigin\b|\bori\b|\bpUC\b|\bColE1\b|\bp15A\b|\bf1\b|\bSV40 ori", "origin"),
    (r"terminat", "terminator"),
    (r"\bintron\b", "intron"),
    (r"traffick|export|\bER2\b|\bKGC\b|\bTS\b signal|golgi|dendritic|\bCAAX\b|membrane target", "trafficking_motif"),
    (r"\blinker\b|\bXTEN\b|\bGS linker\b|flexible", "linker"),
    (r"\bHA tag\b|\bFLAG\b|\bMyc\b|\b6?x?His\b|\bV5\b|\bStrep[- ]?tag", "tag"),
    (r"primer|sequencing site|\bM13\b|\bT7\b|\bSP6\b", "primer_binding_site"),
    (r"backbone|vector", "backbone"),
]

#: GenBank/SnapGene feature keys that already carry meaning.
_KEY_RULES = {
    "cds": "CDS",
    "promoter": "promoter",
    "enhancer": "enhancer",
    "terminator": "terminator",
    "polya_signal": "polyA",
    "polya_site": "polyA",
    "rep_origin": "origin",
    "protein_bind": "lox_site",
    "primer_bind": "primer_binding_site",
    "intron": "intron",
    "ltr": "promoter",
    "sig_peptide": "trafficking_motif",
    "repeat_region": "unknown",
}


#: Keys whose own meaning beats the label rules.
#:
#: Normally the label is the better evidence: a feature keyed ``CDS`` and
#: labelled "AmpR" really is a resistance gene, and that refinement is the
#: point of the label rules.  These two keys are the exceptions.
#:
#: ``None`` means the key wins outright; a set means it wins only against
#: those verdicts.
_KEY_OVERRIDES: dict[str, set[str] | None] = {
    # A primer is named after what it anneals to, never after itself, so
    # SnapGene's stock "Amp-R" and "SV40pA-R" are 20 bp primer sites and not
    # a resistance gene or a polyA signal.  Nothing a primer is called can
    # make it the element it sits next to.
    "primer_bind": None,
    # SnapGene files the WPRE under enhancer, so the label still decides
    # there.  But a depositor who keyed a feature as an enhancer has already
    # distinguished it from a promoter, and "CMV enhancer" should not be
    # overruled into a promoter just because it says CMV.
    "enhancer": {"promoter"},
}


def classify_label(label: str, key: str = "") -> str:
    """Map a free-text label plus GenBank key onto the controlled vocabulary."""
    text = f"{label} {key}".strip()
    k = (key or "").lower()

    verdict = ""
    for pattern, ftype in _LABEL_RULES:
        if re.search(pattern, text, re.IGNORECASE):
            verdict = ftype
            break

    if k in _KEY_OVERRIDES and k in _KEY_RULES:
        beaten = _KEY_OVERRIDES[k]
        if not verdict or beaten is None or verdict in beaten:
            return _KEY_RULES[k]

    if verdict:
        return verdict
    if k in _KEY_RULES:
        mapped = _KEY_RULES[k]
        # repeat_region is only an ITR if the label says so
        if mapped == "unknown" and re.search(r"ITR", label, re.I):
            return "AAV_ITR"
        return mapped
    return "unknown"


# --------------------------------------------------------------------------
# DNA motif library - near-invariant elements
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Motif:
    name: str
    feature_type: str
    seq: str
    level: str = "dna"        # "dna" or "protein"
    max_mismatch: int = 0
    note: str = ""


DNA_MOTIFS: list[Motif] = [
    Motif("loxP", "lox_site", "ATAACTTCGTATAATGTATGCTATACGAAGTTAT"),
    Motif("lox2272", "lox_site", "ATAACTTCGTATAAAGTATCCTATACGAAGTTAT"),
    Motif("loxN", "lox_site", "ATAACTTCGTATAATGTATACTATACGAAGTTAT"),
    Motif("lox511", "lox_site", "ATAACTTCGTATAATGTATACTATACGAAGTTAT"),
    # loxFAS carries mutated inverted repeats (ACAACTTCGTATA / TATACGAAGTTGT)
    # so it recombines with itself but not with loxP. Lab files label it
    # "FAS SITE".
    Motif("loxFAS", "lox_site", "ACAACTTCGTATATACCTTTCTATACGAAGTTGT"),
    Motif("F5", "lox_site", "GAAGTTCCTATTCTTCAAAAGGTATAGGAACTTC"),
    Motif("FRT", "lox_site", "GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC"),
    Motif("Kozak", "primer_binding_site", "GCCACCATGG", note="strong Kozak context"),
    # AAV2 ITR: match the invariant A-A' stem rather than the whole hairpin,
    # which is GC-rich, palindromic and frequently truncated in lab files.
    Motif("AAV2 ITR (A element)", "AAV_ITR",
          "TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGC", max_mismatch=3),
    Motif("bGH polyA", "polyA",
          "CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCG", max_mismatch=3),
    Motif("SV40 polyA", "polyA",
          "AACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAA", max_mismatch=3),
    Motif("WPRE", "WPRE",
          "AATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCC",
          max_mismatch=3),
    Motif("EMCV IRES", "IRES",
          "GCCCCTCTCCCTCCCCCCCCCCTAACGTTACTGGCCGAAGCCGCTTGGAATAAGGCC", max_mismatch=3),
]

# --------------------------------------------------------------------------
# Protein motif library - matched on translation, codon-usage independent
# --------------------------------------------------------------------------

PROTEIN_MOTIFS: list[Motif] = [
    Motif("P2A", "p2A", "ATNFSLLKQAGDVEENPGP", level="protein"),
    Motif("T2A", "p2A", "EGRGSLLTCGDVEENPGP", level="protein"),
    Motif("E2A", "p2A", "QCTNYALLKLAGDVESNPGP", level="protein"),
    Motif("F2A", "p2A", "VKQTLNFDLLKLAGDVESNPGP", level="protein"),
    # Trafficking / export motifs used heavily in opsin and GEVI constructs.
    Motif("ER2 export motif", "trafficking_motif", "FCYENEV", level="protein"),
    Motif("KGC / Kir2.1 TS", "trafficking_motif", "KSRITSEGEYIPLDQIDINV", level="protein"),
    Motif("Kir2.1 ER export", "trafficking_motif", "FCYENEV", level="protein"),
    Motif("CAAX box", "trafficking_motif", "KMSKDGKKKKKKSKTKCVIM", level="protein"),
    Motif("Golgi export (TS)", "trafficking_motif", "KSRITSEGEYIPLDQIDINV", level="protein"),
    Motif("HA tag", "tag", "YPYDVPDYA", level="protein"),
    Motif("FLAG tag", "tag", "DYKDDDDK", level="protein"),
    Motif("Myc tag", "tag", "EQKLISEEDL", level="protein"),
    Motif("V5 tag", "tag", "GKPIPNPLLGLDST", level="protein"),
    Motif("6xHis tag", "tag", "HHHHHH", level="protein"),
    Motif("XTEN linker", "linker", "SGSETPGTSESATPES", level="protein"),
]


def _fuzzy_find(haystack: str, needle: str, max_mismatch: int,
                circular: bool = False) -> list[int]:
    """Find `needle` allowing up to `max_mismatch` substitutions."""
    if max_mismatch == 0:
        return find_all(haystack, needle, circular=circular)
    n, m = len(haystack), len(needle)
    if m == 0 or n == 0:
        return []
    target = haystack + (haystack[:m - 1] if circular else "")

    # Seed-and-extend: split the motif into (max_mismatch + 1) blocks. Any
    # match with <= max_mismatch substitutions must contain at least one
    # block intact, so exact search on blocks gives a cheap candidate list
    # instead of comparing every offset character by character.
    nblocks = max_mismatch + 1
    block = m // nblocks
    hits: set[int] = set()
    if block < 4:
        candidates = range(min(len(target) - m + 1, n))
    else:
        cand: set[int] = set()
        for b in range(nblocks):
            off = b * block
            seed = needle[off:off + block]
            start = target.find(seed)
            while start != -1:
                origin = start - off
                if 0 <= origin < n:
                    cand.add(origin)
                start = target.find(seed, start + 1)
        candidates = sorted(cand)

    for i in candidates:
        window = target[i:i + m]
        if len(window) < m:
            continue
        mism = 0
        for a, b in zip(window, needle):
            if a != b:
                mism += 1
                if mism > max_mismatch:
                    break
        else:
            hits.add(i)
    return sorted(hits)


def find_motifs(seq: str, circular: bool = False,
                include_protein: bool = True) -> list[dict]:
    """Scan a sequence for every known motif, both strands.

    Returns dicts with 0-based ``start``/``end``, ``strand``, ``name``,
    ``feature_type`` and ``level``.
    """
    results: list[dict] = []
    seq = seq.upper()

    for motif in DNA_MOTIFS:
        for strand, target in ((1, seq), (-1, revcomp(seq))):
            for pos in _fuzzy_find(target, motif.seq, motif.max_mismatch, circular):
                start = pos if strand == 1 else len(seq) - pos - len(motif.seq)
                results.append({
                    "name": motif.name, "feature_type": motif.feature_type,
                    "start": start % max(len(seq), 1),
                    "end": (start + len(motif.seq)) % max(len(seq), 1) or len(seq),
                    "strand": strand, "level": "dna", "note": motif.note,
                })

    if include_protein:
        for strand, target in ((1, seq), (-1, revcomp(seq))):
            for frame in range(3):
                prot = translate(target, frame=frame)
                for motif in PROTEIN_MOTIFS:
                    idx = 0
                    while True:
                        idx = prot.find(motif.seq, idx)
                        if idx < 0:
                            break
                        dna_start = frame + idx * 3
                        dna_end = dna_start + len(motif.seq) * 3
                        if strand == 1:
                            s, e = dna_start, dna_end
                        else:
                            s, e = len(seq) - dna_end, len(seq) - dna_start
                        results.append({
                            "name": motif.name, "feature_type": motif.feature_type,
                            "start": s, "end": e, "strand": strand,
                            "level": "protein", "note": motif.note,
                        })
                        idx += 1

    # Deduplicate identical spans, preferring protein-level evidence.
    seen: dict[tuple, dict] = {}
    for r in results:
        key = (r["start"], r["end"], r["strand"], r["name"])
        if key not in seen or r["level"] == "protein":
            seen[key] = r
    return sorted(seen.values(), key=lambda r: (r["start"], r["end"]))


def detect_dio(features: list[dict]) -> list[dict]:
    """Detect DIO / FLEX architecture: two pairs of inverted lox sites.

    A DIO cassette is a CDS flanked by heterotypic lox pairs (e.g. loxP and
    lox2272) in opposing orientations. Finding this matters because the
    expressed ORF only exists *after* Cre inversion.
    """
    lox = [f for f in features if f.get("feature_type") == "lox_site"]
    findings = []
    by_name: dict[str, list[dict]] = {}
    for site in lox:
        by_name.setdefault(site["name"], []).append(site)

    for name, sites in by_name.items():
        if len(sites) < 2:
            continue
        for i in range(len(sites)):
            for j in range(i + 1, len(sites)):
                a, b = sites[i], sites[j]
                if a["strand"] != b["strand"]:
                    findings.append({
                        "type": "inverted_lox_pair", "site": name,
                        "left": min(a["start"], b["start"]),
                        "right": max(a["end"], b["end"]),
                        "detail": f"inverted {name} pair flanking "
                                  f"{abs(b['start'] - a['end'])} bp",
                    })
    return findings


# --------------------------------------------------------------------------
# Composite parts
# --------------------------------------------------------------------------

#: Element types that are commonly fused into a CDS annotation rather than
#: being annotated separately.
EMBEDDABLE_TYPES = ("p2A", "trafficking_motif", "tag", "linker")


def find_embedded_elements(seq: str, feature_type: str = "CDS") -> list[dict]:
    """Find elements fused inside a part's own sequence.

    Lab files frequently annotate a whole fusion as one CDS - "Voltron2"
    whose sequence actually ends in a P2A, or a CDS that already carries its
    trafficking motif. Treating that as a plain CDS and then adding an
    explicit P2A module would build a construct with two tandem 2A peptides,
    so the embedded elements have to be found and reported.

    Matching is done on the translation in all three frames, because the
    peptide is conserved where the codons are not.
    """
    if feature_type not in ("CDS", "Cre-lox", "DIO") or len(seq) < 30:
        return []

    found: list[dict] = []
    for frame in range(3):
        prot = translate(seq, frame=frame)
        if not prot:
            continue
        for motif in PROTEIN_MOTIFS:
            if motif.feature_type not in EMBEDDABLE_TYPES:
                continue
            idx = 0
            while True:
                idx = prot.find(motif.seq, idx)
                if idx < 0:
                    break
                start = frame + idx * 3
                end = start + len(motif.seq) * 3
                found.append({
                    "name": motif.name, "feature_type": motif.feature_type,
                    "start": start, "end": end, "frame": frame,
                    "at_end": end >= len(seq) - 6,
                    "at_start": start <= 6,
                })
                idx += 1

    # Collapse hits that describe the same span.
    unique: dict[tuple, dict] = {}
    for f in found:
        unique.setdefault((f["start"], f["end"], f["feature_type"]), f)
    return sorted(unique.values(), key=lambda f: f["start"])


def trim_embedded(seq: str, elements: list[dict], feature_type: str) -> tuple[str, list[dict]]:
    """Remove an embedded element that sits at one end of the sequence.

    Only terminal elements are trimmed: cutting something out of the middle
    of a coding sequence would change the protein, which is never what the
    user meant.
    """
    targets = [e for e in elements if e["feature_type"] == feature_type
               and (e["at_end"] or e["at_start"])]
    if not targets:
        return seq, []
    trimmed = seq
    removed = []
    for e in sorted(targets, key=lambda e: -e["start"]):
        if e["at_end"]:
            trimmed = trimmed[:e["start"]]
            removed.append(e)
        elif e["at_start"]:
            trimmed = trimmed[e["end"]:]
            removed.append(e)
    return trimmed, removed


#: Elements that separate one protein from the next inside a fused CDS.
#: A 2A peptide genuinely ends one protein and starts another; a linker or
#: trafficking signal between two coding blocks is doing the same job in
#: practice, which is why lab annotations run them all together.
SPLIT_TYPES = ("p2A", "linker", "trafficking_motif")

#: Below this a segment is a fragment, not a protein worth cataloguing.
MIN_SEGMENT_BP = 150


def split_at_elements(seq: str, elements: list[dict],
                      split_types: tuple = SPLIT_TYPES,
                      min_segment: int = MIN_SEGMENT_BP) -> list[dict]:
    """Cut a fused coding sequence into its separate proteins.

    ``LR-Voltron2-LR-CheRiff_eYFP`` is annotated as one CDS but is really
    two proteins joined by a 2A peptide. Cutting at that boundary recovers
    both, while the uncut original stays in the catalogue because sometimes
    the whole fusion is exactly what you want.

    Separators that sit next to each other are treated as one boundary.
    Voltron2 ends in KGC, then ER2, then the 2A: those are its own
    C-terminal trafficking tags followed by the real protein boundary, not
    three boundaries with slivers of protein between them. A boundary is
    only taken where a full-length coding block lies on each side.

    Returns pieces in order, each ``{kind, feature_type, name, start, end,
    seq}`` with ``kind`` of ``"segment"`` or ``"separator"``. Returns
    nothing unless at least two real coding segments come out, since a
    single tag on the end of one protein is not a fusion.
    """
    candidates = sorted(
        (e for e in elements if e["feature_type"] in split_types),
        key=lambda e: e["start"])
    if not candidates:
        return []

    # Group separators that are adjacent, or separated by less coding
    # sequence than could be a protein, into a single boundary.
    clusters: list[list[dict]] = []
    for element in candidates:
        if clusters and element["start"] - clusters[-1][-1]["end"] < min_segment:
            clusters[-1].append(element)
        else:
            clusters.append([element])

    pieces: list[dict] = []
    cursor = 0
    for cluster in clusters:
        # Where a cluster contains a 2A or an IRES, that element is the real
        # protein boundary and anything beside it belongs to the protein it
        # is attached to: Voltron2's KGC and ER2 are its own C-terminal
        # trafficking tags, not part of the junction to the next protein.
        hard = [e for e in cluster if e["feature_type"] in ("p2A", "IRES")]
        if hard:
            lo, hi = hard[0]["start"], hard[-1]["end"]
            cluster = hard
        else:
            lo = cluster[0]["start"]
            hi = cluster[-1]["end"]
        before = lo - cursor
        after = len(seq) - hi
        # Only a boundary if a whole protein sits on each side of it.
        if before < min_segment or after < min_segment:
            continue
        pieces.append({"kind": "segment", "feature_type": "CDS", "name": "",
                       "start": cursor, "end": lo, "seq": seq[cursor:lo]})
        pieces.append({
            "kind": "separator",
            "feature_type": cluster[0]["feature_type"],
            "name": " + ".join(dict.fromkeys(e["name"] for e in cluster)),
            "start": lo, "end": hi, "seq": seq[lo:hi],
            "elements": [e["name"] for e in cluster],
        })
        cursor = hi

    if not pieces:
        return []
    if cursor < len(seq):
        pieces.append({"kind": "segment", "feature_type": "CDS", "name": "",
                       "start": cursor, "end": len(seq), "seq": seq[cursor:]})

    coding = [p for p in pieces if p["kind"] == "segment"
              and len(p["seq"]) >= min_segment]
    if len(coding) < 2:
        return []
    return pieces
