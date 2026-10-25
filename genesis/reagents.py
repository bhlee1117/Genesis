"""The kits this lab actually buys, and how each one is run.

Genesis used to print a generic recipe per operation - "2x assembly
master mix, 10 uL" - which is right for every kit and exact for none.
The purchasing record says which products the bench really uses, and
those have published conditions that differ from the generic ones in
ways that matter: NEB Stable cells recover at 30 C in their own medium
rather than at 37 C in SOC, Phusion extends plasmid template at 15 s/kb
rather than 30, and an HF enzyme pair in rCutSmart is Time-Saver
qualified and cuts in 15 minutes rather than an hour.

What is here is the reaction setup: volumes, temperatures and times.
Those are facts about a reaction, not the vendor's text, and each entry
carries the catalogue number and the URL of the page it came from so
the original - including the vendor's own figures, which are theirs and
are not reproduced here - is one click away.

Checked against the vendor protocols on 2026-10-25.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Product:
    """One catalogue item and the reaction it is used in."""
    vendor: str
    catalogue: str
    name: str
    url: str
    #: (amount, component) - the order a bench protocol lists them.
    setup: tuple = ()
    title: str = ""
    #: (stage, condition) for a thermocycler or incubation programme.
    program: tuple = ()
    notes: tuple = ()
    #: Operations this product is the default for.
    operations: tuple = ()
    #: Seconds of 72 C extension per kb on a plasmid template. The rate
    #: belongs to the polymerase - Q5 wants 20-30 s/kb, Phusion 15 on
    #: the same template - so the figure Genesis prints for a step is
    #: taken from here rather than written out a second time.
    extension_s_per_kb: int = 0

    @property
    def label(self) -> str:
        return f"{self.vendor} #{self.catalogue}  {self.name}"


#: Everything below is keyed by what the purchasing record shows. Where
#: the lab owns two products for one job the more recent and more often
#: reordered one is the default and the other is listed as an
#: alternative, so the protocol names something that is on the shelf.
CATALOGUE: tuple = (
    Product(
        vendor="NEB", catalogue="M0492", name="Q5 High-Fidelity 2X Master Mix",
        url="https://www.neb.com/en-us/protocols/"
            "protocol-for-q5-high-fidelity-2x-master-mix-m0492",
        title="Q5 High-Fidelity 2X Master Mix (25 uL)",
        operations=("PCR",), extension_s_per_kb=30,
        setup=(("12.5 uL", "Q5 High-Fidelity 2X Master Mix"),
               ("1.25 uL", "forward primer, 10 uM (0.5 uM final)"),
               ("1.25 uL", "reverse primer, 10 uM (0.5 uM final)"),
               ("1 pg - 1 ng", "plasmid template"),
               ("to 25 uL", "nuclease-free water")),
        program=(("initial denaturation", "98 C, 30 s"),
                 ("25-35 cycles", "98 C 5-10 s / anneal 10-30 s / "
                                  "72 C 20-30 s per kb"),
                 ("final extension", "72 C, 2 min"),
                 ("hold", "4-10 C")),
        notes=("Assemble on ice and move to a block already at 98 C.",
               "Anneal 3 C above the lower primer Tm, from the NEB Tm "
               "calculator with Q5 selected - Q5 anneals hotter than "
               "Phusion for the same primers.",
               "Above 6 kb, extend at 40-50 s per kb."),
    ),
    Product(
        vendor="NEB", catalogue="M0531",
        name="Phusion High-Fidelity PCR Master Mix with HF Buffer",
        url="https://www.neb.com/en/protocols/"
            "protocol-phusion-high-fidelity-pcr-master-mix-with-hf-buffer-m0531",
        title="Phusion High-Fidelity PCR Master Mix, HF Buffer (50 uL)",
        operations=(), extension_s_per_kb=15,
        setup=(("25 uL", "Phusion High-Fidelity PCR Master Mix, 2X"),
               ("2.5 uL", "forward primer, 10 uM (0.5 uM final)"),
               ("2.5 uL", "reverse primer, 10 uM (0.5 uM final)"),
               ("1 pg - 10 ng", "plasmid template"),
               ("to 50 uL", "nuclease-free water")),
        program=(("initial denaturation", "98 C, 30 s"),
                 ("25-35 cycles", "98 C 5-10 s / anneal 10-30 s / "
                                  "72 C 15 s per kb"),
                 ("final extension", "72 C, 5-10 min"),
                 ("hold", "4-10 C")),
        notes=("15 s per kb is for low-complexity template - plasmid, "
               "lambda, BAC. Genomic DNA wants 30 s per kb.",
               "Anneal at the lower primer Tm + 3 C for primers over "
               "20 nt, or at the lower Tm itself for shorter ones."),
    ),
    Product(
        vendor="NEB", catalogue="E2621",
        name="NEBuilder HiFi DNA Assembly Master Mix",
        url="https://www.neb.com/en-us/products/"
            "e2621-nebuilder-hifi-dna-assembly-master-mix",
        title="NEBuilder HiFi DNA Assembly (20 uL)",
        operations=("HI_FI_ASSEMBLY", "GIBSON_ASSEMBLY"),
        setup=(("10 uL", "NEBuilder HiFi DNA Assembly Master Mix, 2X"),
               ("X uL", "linearised vector and inserts (see note)"),
               ("to 20 uL", "nuclease-free water")),
        program=(("2-3 fragments", "50 C, 15 min"),
                 ("4-6 fragments", "50 C, 60 min"),
                 ("then", "on ice, or -20 C until transformation")),
        notes=("2-3 fragments: 1:2 vector:insert, 0.03-0.2 pmol of "
               "fragments in total. 4-6 fragments: 1:1, 0.2-0.5 pmol.",
               "Under 15 min gains nothing and loses yield.",
               "Transform 2 uL of the finished assembly."),
    ),
    Product(
        vendor="NEB", catalogue="M2200", name="Quick Ligation Kit",
        url="https://www.neb.com/en-us/protocols/quick-ligation-protocol",
        title="Quick Ligation (20 uL)",
        operations=("LIGATION",),
        setup=(("10 uL", "Quick Ligase Reaction Buffer, 2X"),
               ("50 ng", "digested vector"),
               ("3-fold molar excess", "insert over vector"),
               ("1 uL", "Quick Ligase (added last)"),
               ("to 20 uL", "nuclease-free water")),
        program=(("ligation", "25 C, 5 min"),
                 ("then", "chill on ice and transform 1-5 uL")),
        notes=("Thaw and resuspend the 2X buffer at room temperature; "
               "set the reaction up on ice and add the ligase last.",
               "Do not heat inactivate - it collapses the "
               "transformation efficiency.",
               "Past 5 min at 25 C the reaction has finished.",
               "Keep vector + insert between 1 and 10 ng/uL; below "
               "1 ng/uL the fragments circularise on themselves.",
               "Run a vector-only ligation alongside, as the background."),
    ),
    Product(
        vendor="NEB", catalogue="R0176", name="DpnI",
        url="https://www.neb.com/en-us/products/r0176-dpni",
        title="DpnI treatment",
        operations=("DPNI_DIGEST",),
        setup=(("all", "completed PCR reaction"),
               ("1 uL", "DpnI (20 U)")),
        program=(("digest", "37 C, 15-60 min"),),
        notes=("Added straight to the PCR - no clean-up first.",
               "DpnI reads methylation rather than a rare sequence: it "
               "cuts GATC only where the adenine is dam-methylated, so "
               "it destroys a plasmid grown in a standard strain and "
               "leaves the PCR product, built from unmethylated dNTPs, "
               "untouched. Cutting the template at many sites is how "
               "the step works."),
    ),
    Product(
        vendor="NEB", catalogue="rCutSmart",
        name="HF restriction enzymes in rCutSmart Buffer",
        url="https://www.neb.com/en/tools-and-resources/usage-guidelines/"
            "double-digests",
        title="Restriction digest (50 uL)",
        operations=("RE_DIGEST",),
        setup=(("1 ug", "DNA"),
               ("5 uL", "rCutSmart Buffer, 10X"),
               ("1 uL each", "enzyme (10-20 U)"),
               ("to 50 uL", "nuclease-free water")),
        program=(("Time-Saver pair", "37 C, 15 min"),
                 ("otherwise", "37 C, 1 h")),
        notes=("Both HF enzymes together in one rCutSmart reaction; "
               "5-10 units per ug of DNA.",
               "Keep the total enzyme under 5 uL in 50 uL - above 5 % "
               "glycerol the enzymes start to show star activity.",
               "A long digest is the other way to get star activity; "
               "a Time-Saver pair does not need one.",
               "Check the band on a gel before going on."),
    ),
    Product(
        vendor="NEB", catalogue="C3040",
        name="NEB Stable Competent E. coli (High Efficiency)",
        url="https://www.neb.com/en/protocols/2013/10/30/"
            "high-efficiency-transformation-protocol-c3040h",
        title="Transformation, NEB Stable",
        operations=("TRANSFORMATION",),
        setup=(("50 uL", "NEB Stable Competent E. coli, thawed on ice 10 min"),
               ("1-2 uL", "DNA (100 pg - 100 ng)"),
               ("950 uL", "NEB 10-beta / Stable Outgrowth Medium")),
        program=(("on ice", "30 min"),
                 ("heat shock", "42 C, exactly 30 s"),
                 ("on ice", "5 min"),
                 ("outgrowth", "30 C, 60 min, shaking at 250 rpm"),
                 ("plate", "50-100 uL; 24 h at 30 C")),
        notes=("30 C, not 37 C. This strain is bought for repeats and "
               "inverted terminal repeats - an AAV transfer plasmid is "
               "exactly that - and the lower temperature is what keeps "
               "them from recombining out.",
               "Use the outgrowth medium supplied with the cells. NEB "
               "says not to substitute SOC for this strain.",
               "Mix by flicking, never by vortexing.",
               "Miniprep from a fresh transformant, off a plate no more "
               "than three days old, and keep an unstable construct as "
               "DNA rather than as a glycerol stock."),
    ),
)

#: Other things on the shelf for the same jobs, named so that a protocol
#: can be run with what is actually in the freezer.
ALTERNATIVES = {
    "PCR": ("NEB #M0531 Phusion High-Fidelity PCR Master Mix with HF Buffer",
            "Thermo #F548S Phusion Flash High-Fidelity PCR Master Mix"),
    "TRANSFORMATION": ("NEB #C2987H NEB 5-alpha Competent E. coli (High "
                       "Efficiency) - for a construct with no repeats, "
                       "recovered at 37 C in SOC",
                       "NEB #C2988J NEB 5-alpha Subcloning Efficiency"),
    "GEL_EXTRACTION": ("Thermo #K0691 GeneJET Gel Extraction Kit",),
    "DNA_CLEANUP": ("Thermo #K0691 GeneJET Gel Extraction Kit (column "
                    "clean-up without a gel)",),
}

_BY_OPERATION = {}
for _p in CATALOGUE:
    for _op in _p.operations:
        _BY_OPERATION.setdefault(_op, _p)


def for_operation(operation: str) -> Product | None:
    """The product this lab runs that operation with, if one is known."""
    return _BY_OPERATION.get(operation)


def alternatives(operation: str) -> tuple:
    return ALTERNATIVES.get(operation, ())


def used_by(operations) -> list:
    """The catalogue entries a plan touches, in protocol order."""
    seen, out = set(), []
    for op in operations:
        product = _BY_OPERATION.get(op)
        if product is not None and product.catalogue not in seen:
            seen.add(product.catalogue)
            out.append(product)
    return out
