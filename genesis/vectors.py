"""Vector families: what a finished AAV or lentiviral plasmid must contain.

A recipient plasmid is not just a backbone with a promoter and a CDS. An
AAV transfer plasmid carries two ITRs flanking the payload, a WPRE and a
polyA inside them, and outside them the bacterial machinery that lets the
plasmid be grown and selected at all. A lentiviral transfer plasmid carries
LTRs, the packaging signal, RRE and cPPT instead.

Getting this wrong is quiet and expensive: a construct missing its WPRE
still clones fine and still expresses, just worse, and a construct missing
an ITR is simply not packageable. So the family is detected from the
features a plasmid actually has, and a design that claims to be one is
checked against the list.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass(frozen=True)
class VectorFamily:
    key: str
    label: str
    #: Features that must be present for the vector to work at all.
    required: tuple
    #: Features that are normally present and whose absence is worth saying.
    expected: tuple
    #: Label patterns that identify this family when the feature types alone
    #: are ambiguous.
    #: (pattern, plain-English description) pairs. The description is what
    #: a person reads; nobody should be shown a regex as evidence.
    signatures: tuple
    note: str = ""


#: Names as they appear in lab annotations, matched case-insensitively.
AAV = VectorFamily(
    key="AAV",
    label="AAV transfer plasmid",
    required=("AAV_ITR", "promoter", "polyA", "origin", "resistance"),
    expected=("WPRE",),
    signatures=((r"\bITR\b", "the name mentions an ITR"),
                (r"inverted terminal repeat", "the name mentions an inverted "
                                              "terminal repeat"),
                (r"\bpAAV", "the name begins pAAV"),
                (r"\bAAV\d?\b", "the name mentions AAV")),
    note="payload between the ITRs must stay under about 4.7 kb",
)

LENTI = VectorFamily(
    key="lentiviral",
    label="lentiviral transfer plasmid",
    required=("promoter", "polyA", "origin", "resistance"),
    expected=("WPRE",),
    signatures=((r"lenti", "the name says lenti"),
                (r"\bHIV-?1\b", "the name mentions HIV-1"),
                (r"\bFUGW\b|\bFUW\b", "the name is a FUGW-family vector"),
                (r"\bFCK\b", "the name is an FCK-family vector"),
                (r"pLJM|pHR", "the name is a known lentiviral backbone"),
                (r"\bLTR\b", "the name mentions an LTR"),
                (r"\bRRE\b|\bcPPT\b|\bpsi\b", "the name mentions a "
                                                  "packaging element")),
    note="LTRs, the packaging signal, RRE and cPPT are what make it package",
)

PLAIN = VectorFamily(
    key="plasmid",
    label="plain plasmid",
    required=("origin", "resistance"),
    expected=(),
    signatures=(),
)

FAMILIES = {f.key: f for f in (AAV, LENTI, PLAIN)}

#: Lentiviral elements that are annotated as free text rather than as a
#: controlled feature type, so they have to be matched by name.
LENTI_ELEMENTS = {
    "5' LTR": (r"5'? ?LTR", r"RSV.*LTR", r"LTR.*hybrid"),
    "3' LTR": (r"3'? ?LTR", r"\bSIN\b.*LTR"),
    "psi packaging signal": (r"\bpsi\b", r"packaging signal", r"\bHIV-1 psi\b"),
    "RRE": (r"\bRRE\b", r"rev response"),
    "cPPT": (r"\bcPPT\b", r"central polypurine"),
}

AAV_ELEMENTS = {
    "left ITR": (r"5'? ?ITR", r"left ITR", r"\bL-?ITR\b"),
    "right ITR": (r"3'? ?ITR", r"right ITR", r"\bR-?ITR\b"),
}


def _matches(patterns, text: str) -> bool:
    return any(re.search(p, text, re.IGNORECASE) for p in patterns)


def classify_vector(annotations: list[dict], name: str = "") -> dict:
    """Work out which family a plasmid belongs to, and why.

    Both families are scored rather than short-circuiting on the first
    signal, because the two sets of evidence can genuinely conflict: a
    plasmid named pLenti that also carries ITR annotations is either a
    hybrid or mis-annotated, and either way the person should be told
    rather than handed a confident wrong answer.
    """
    labels = " | ".join(a.get("name") or "" for a in annotations)
    haystack = f"{name} | {labels}"

    aav_evidence: list[str] = []
    lenti_evidence: list[str] = []

    n_itr = sum(1 for a in annotations if a.get("feature_type") == "AAV_ITR")
    if n_itr >= 2:
        aav_evidence.append(f"{n_itr} ITR annotations")
    elif n_itr == 1:
        aav_evidence.append("one ITR annotation, which is not packageable alone")

    lenti_hits = [k for k, pats in LENTI_ELEMENTS.items()
                  if _matches(pats, haystack)]
    if lenti_hits:
        lenti_evidence.append("lentiviral elements: " + ", ".join(lenti_hits))

    for pattern, description in AAV.signatures:
        if re.search(pattern, name, re.IGNORECASE):
            aav_evidence.append(description)
            break
    for pattern, description in LENTI.signatures:
        if re.search(pattern, name, re.IGNORECASE):
            lenti_evidence.append(description)
            break

    # Weight structural evidence above names: an annotation describes the
    # molecule, a filename describes whatever someone typed.
    aav_name_hit = any(e.startswith("the name") for e in aav_evidence)
    lenti_name_hit = any(e.startswith("the name") for e in lenti_evidence)
    aav_score = (2 if n_itr >= 2 else 0) + (1 if n_itr == 1 else 0) \
        + (1 if aav_name_hit else 0)
    lenti_score = 2 * (len(lenti_hits) >= 2) + (1 if len(lenti_hits) == 1 else 0) \
        + (1 if lenti_name_hit else 0)

    conflict = aav_score >= 2 and lenti_score >= 2
    if aav_score > lenti_score:
        family, evidence = AAV, aav_evidence
    elif lenti_score > aav_score:
        family, evidence = LENTI, lenti_evidence
    elif aav_score == lenti_score and aav_score > 0:
        family = AAV if n_itr >= 2 else LENTI
        evidence = aav_evidence + lenti_evidence
        conflict = True
    else:
        family, evidence = PLAIN, []

    result = {"family": family.key, "label": family.label,
              "evidence": evidence, "note": family.note,
              "conflict": conflict}
    if conflict:
        result["warning"] = (
            "this plasmid shows both AAV and lentiviral evidence "
            f"(AAV: {'; '.join(aav_evidence) or 'none'}; "
            f"lentiviral: {'; '.join(lenti_evidence) or 'none'}); "
            f"treating it as {family.label}, but check the annotation")
    return result


def missing_elements(annotations: list[dict], family_key: str,
                     name: str = "") -> dict:
    """Which parts of a vector of this family are absent."""
    family = FAMILIES.get(family_key, PLAIN)
    types = {a.get("feature_type") for a in annotations}
    haystack = f"{name} | " + " | ".join(a.get("name") or ""
                                         for a in annotations)

    missing_required = [t for t in family.required if t not in types]
    missing_expected = [t for t in family.expected if t not in types]

    # An AAV plasmid needs two ITRs, not just the feature type present once.
    if family is AAV:
        n_itr = sum(1 for a in annotations
                    if a.get("feature_type") == "AAV_ITR")
        if 0 < n_itr < 2 and "AAV_ITR" not in missing_required:
            missing_required.append("second AAV_ITR")

    named_missing = []
    catalogue = AAV_ELEMENTS if family is AAV else (
        LENTI_ELEMENTS if family is LENTI else {})
    # Two ITRs annotated by feature type already satisfy "left and right";
    # most lab files just call both of them "AAV2 ITR" with no side given.
    if family is AAV and sum(
            1 for a in annotations
            if a.get("feature_type") == "AAV_ITR") >= 2:
        catalogue = {}
    for element, patterns in catalogue.items():
        if not _matches(patterns, haystack):
            named_missing.append(element)

    return {"family": family.key, "label": family.label,
            "missing_required": missing_required,
            "missing_expected": missing_expected,
            "missing_named": named_missing,
            "complete": not missing_required and not named_missing}


def describe(family_key: str) -> dict:
    family = FAMILIES.get(family_key, PLAIN)
    return {"key": family.key, "label": family.label,
            "required": list(family.required),
            "expected": list(family.expected), "note": family.note}
