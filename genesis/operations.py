"""Stage 6: the laboratory operation library.

Every step in a generated plan carries an explicit operation code, so the
output reads as bench instructions rather than as an abstract graph.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict


@dataclass(frozen=True)
class Operation:
    """One laboratory operation.

    Two separate time estimates, because conflating them is misleading:
    ``default_minutes`` is elapsed turnaround (a sequencing run is a day,
    a synthesis order is a week) while ``hands_on_minutes`` is the active
    bench time a person actually spends. Planning a week of work needs the
    first; deciding what to do this afternoon needs the second.
    """
    code: str
    meaning: str
    inputs: tuple
    output: str
    mandatory: bool = False
    #: Worth doing and not required. A checkpoint gel tells you whether
    #: the last reaction worked; skipping it does not stop the cloning,
    #: it only moves the moment you find out. Marked so the protocol can
    #: say so, and so an optional step's bench time is not added to the
    #: hours the protocol actually demands.
    optional: bool = False
    changes_sequence: bool = False
    combinable_with: tuple = ()
    default_minutes: int = 60      # elapsed, including incubations and waits
    hands_on_minutes: int = 30     # active bench time
    notes: str = ""


OPERATIONS: dict[str, Operation] = {o.code: o for o in [
    Operation("PCR", "Amplify a DNA fragment and/or add overlap sequences",
              ("template plasmid", "primer pair"), "amplicon",
              changes_sequence=True, default_minutes=90, hands_on_minutes=20,
              combinable_with=("DPNI_DIGEST",),
              notes="High-fidelity polymerase; anneal ~3 C below the lower primer Tm."),
    Operation("RE_DIGEST", "Restriction-enzyme digestion",
              ("DNA", "enzyme(s)", "buffer"), "digested DNA",
              default_minutes=60, hands_on_minutes=10, combinable_with=("DNA_CLEANUP",),
              notes="Check buffer compatibility before setting up a double digest."),
    Operation("GEL_EXTRACTION", "Separate and recover the desired DNA band",
              ("digested DNA or PCR product",), "size-selected DNA",
              default_minutes=75, hands_on_minutes=45,
              notes="Use when several bands are expected or a parental band must be removed."),
    Operation("DNA_CLEANUP", "Column or bead cleanup",
              ("PCR or digest reaction",), "clean DNA", default_minutes=20, hands_on_minutes=15),
    Operation("GIBSON_ASSEMBLY", "Isothermal assembly of overlapping fragments",
              ("backbone", "inserts"), "assembled plasmid",
              changes_sequence=True, default_minutes=60, hands_on_minutes=15,
              notes="50 C, 15-60 min; 2-3 fold molar excess of insert over vector."),
    Operation("HI_FI_ASSEMBLY", "HiFi-style homologous assembly",
              ("backbone", "inserts"), "assembled plasmid",
              changes_sequence=True, default_minutes=60, hands_on_minutes=15,
              notes="NEBuilder HiFi; 50 C, 15 min for 2-3 fragments."),
    Operation("LIGATION", "Enzyme-mediated joining of compatible ends",
              ("digested vector", "insert"), "ligated plasmid",
              changes_sequence=True, default_minutes=60, hands_on_minutes=10,
              notes="T4 DNA ligase; 3:1 insert:vector molar ratio."),
    Operation("DEPHOSPHORYLATION", "Remove 5' phosphates to suppress self-ligation",
              ("digested vector",), "dephosphorylated vector",
              default_minutes=30, hands_on_minutes=10, combinable_with=("RE_DIGEST",),
              notes="Recommended when the vector is cut with a single enzyme."),
    Operation("DPNI_DIGEST", "Remove methylated parental template after PCR",
              ("PCR reaction",), "template-free amplicon", default_minutes=30, hands_on_minutes=5,
              combinable_with=("PCR",),
              notes="Only works on template grown in a dam+ strain."),
    Operation("TRANSFORMATION", "Introduce assembled DNA into bacteria",
              ("assembly or ligation product",), "bacterial colonies",
              mandatory=True, default_minutes=1020, hands_on_minutes=30),
    Operation("COLONY_SCREEN", "Screen colonies by PCR or restriction digest",
              ("colony",), "candidate clone", mandatory=True, default_minutes=960, hands_on_minutes=60),
    Operation("SEQUENCE_VERIFY", "Confirm the final construct by sequencing",
              ("candidate clone",), "verified plasmid", mandatory=True,
              default_minutes=1440, hands_on_minutes=20),
    Operation("QC_GEL", "Check the last reaction worked before going on",
              ("an aliquot of the reaction",), "a go / no-go answer",
              default_minutes=45, hands_on_minutes=20, optional=True,
              notes="Optional, and the cheapest step in the protocol and "
                    "the one that saves the most: a failed PCR costs an "
                    "hour here and two days if it is found on the plate "
                    "instead."),
    Operation("SYNTHESIS", "Order a synthetic DNA fragment",
              ("sequence",), "synthetic fragment", changes_sequence=True,
              default_minutes=10080, hands_on_minutes=10, notes="Use only when no lab template exists."),
]}


@dataclass
class Step:
    """One concrete bench step in a generated plan."""
    index: int
    operation: str
    description: str
    inputs: list = field(default_factory=list)
    output: str = ""
    details: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        op = OPERATIONS.get(self.operation)
        d["meaning"] = op.meaning if op else ""
        d["notes"] = op.notes if op else ""
        d["minutes"] = op.default_minutes if op else 0       # elapsed
        d["hands_on"] = op.hands_on_minutes if op else 0     # active bench time
        d["optional"] = bool(op.optional) if op else False
        return d


def renumber(steps: list[Step]) -> list[Step]:
    for i, s in enumerate(steps, start=1):
        s.index = i
    return steps
