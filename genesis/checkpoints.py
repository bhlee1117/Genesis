"""Stop and look, before the next step spends a day hiding the answer.

A cloning protocol fails most often at the PCR, and the failure is
silent: the amplicon that is not there goes into the assembly anyway,
the assembly goes into cells anyway, and the first hard evidence
arrives two days later as a plate of colonies that are all parental
vector. The information needed to catch it existed at minute ninety.

So each step that makes a molecule gets a checkpoint after it, and the
checkpoint is written against *this* plan rather than in general. It
knows the band to expect because the planner computed it, the template
band that would mean the digest did not go, the annealing temperature
to drop if nothing comes up, and the no-insert control to compare the
plate against. "Run a gel" is advice; "one band at 5,791 bp, and a band
at 6,539 bp is uncut template" is a test.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Checkpoint:
    """What to do after a step to find out whether it worked."""
    operation: str                 # the operation this follows
    title: str
    how: tuple = ()                # what to set up
    expect: tuple = ()             # what a good result looks like
    #: (what you see, what it means, what to do next)
    troubles: tuple = ()
    minutes: int = 45
    hands_on: int = 20

    def to_dict(self) -> dict:
        return {"operation": self.operation, "title": self.title,
                "how": list(self.how), "expect": list(self.expect),
                "troubles": [list(t) for t in self.troubles],
                "minutes": self.minutes, "hands_on": self.hands_on}


def _gel_percent(length: int) -> str:
    """The gel that resolves a band of this size from its neighbours."""
    if not length:
        return "1 %"
    if length < 500:
        return "2 %"
    if length < 1000:
        return "1.5 %"
    if length < 4000:
        return "1 %"
    return "0.8 %"


def after_pcrs(steps: list, route: dict, templates: dict) -> Checkpoint:
    """One gel for a run of PCRs, with a lane and a size for each.

    Two amplifications set up back to back are run on one gel - nobody
    pours a second gel to look at the second tube - so the checkpoint
    covers the run and names every band it should show.
    """
    if len(steps) == 1:
        return after_pcr(steps[0], route, templates)

    rows, longest = [], 0
    for step in steps:
        details = step.get("details") or {}
        length = details.get("length") or 0
        longest = max(longest, length)
        name = (step.get("output") or "amplicon").replace(" amplicon", "")
        rows.append((name, length, details.get("template") or ""))

    expect = [f"{name}: one band at {length:,} bp" if length
              else f"{name}: one clean band" for name, length, _t in rows]
    for _name, _length, template in rows:
        tpl = templates.get(template)
        if tpl:
            expect.append(f"nothing at {tpl:,} bp - that size is {template} "
                          f"itself, carried through from the template")
    base = after_pcr(steps[0], route, templates)
    return Checkpoint(
        operation="PCR",
        title=("Check both PCRs on one gel before going on"
               if len(steps) == 2 else
               f"Check all {len(steps)} PCRs on one gel before going on"),
        how=(f"Run 3 uL of each reaction in its own lane on a "
             f"{_gel_percent(longest)} agarose gel, against a ladder that "
             f"spans {longest:,} bp.",
             "Keep the rest of each reaction; this is an analytical gel, "
             "not a preparative one."),
        expect=tuple(dict.fromkeys(expect)),
        troubles=base.troubles,
        minutes=45, hands_on=25,
    )


def after_pcr(step: dict, route: dict, templates: dict) -> Checkpoint:
    """The gel that says whether the amplicon exists, and what it is.

    Everything here is a number the planner already has: the product
    length, the template length - which is the band that means the
    polymerase copied the plasmid rather than the region - and the
    annealing temperature, which is the first thing to change.
    """
    details = step.get("details") or {}
    length = details.get("length") or 0
    name = (step.get("output") or "the amplicon").replace(" amplicon", "")
    template = details.get("template") or ""
    tpl_bp = templates.get(template)
    gc = details.get("gc")

    expect = [f"one band at {length:,} bp" if length else "one clean band"]
    if tpl_bp:
        expect.append(
            f"nothing at {tpl_bp:,} bp - that size is {template} itself, "
            f"carried through from the template")
    expect.append("no smear, and nothing below ~100 bp beyond faint primer")

    troubles = [
        ("no band at all",
         "the primers did not prime, or there was no template",
         "re-run with the annealing temperature 5 C lower, or as a "
         "gradient across +/- 5 C; check the template on the same gel"),
        ("a smear rather than a band",
         "too much template, too many cycles, or annealing too low",
         "drop the template tenfold, cut to 25 cycles, raise the "
         "annealing temperature 3-5 C"),
        ("several bands",
         "the primers also prime somewhere else",
         "raise the annealing temperature, or gel purify the band at "
         + (f"{length:,} bp" if length else "the right size")
         + " and carry that forward"),
        (f"one band, but not at {length:,} bp" if length
         else "a band at the wrong size",
         "the primers found a different pair of sites",
         "check each primer against the template before changing the "
         "cycling; a wrong-size band is a sequence problem, not a "
         "temperature one"),
    ]
    if gc and gc > 0.65:
        troubles.append(
            (f"weak product from a {gc:.0%} GC template",
             "GC-rich template needs help melting",
             "add the high-GC enhancer, or raise the denaturation to 98 C "
             "for 15 s"))

    return Checkpoint(
        operation="PCR",
        title=f"Check the {name} PCR on a gel before going on",
        how=(f"Run 3 uL of the reaction on a {_gel_percent(length)} agarose "
             f"gel against a ladder that spans {length:,} bp."
             if length else
             "Run 3 uL of the reaction on an agarose gel against a ladder.",
             "Keep the rest of the reaction; this is an analytical gel, "
             "not a preparative one."),
        expect=tuple(expect),
        troubles=tuple(troubles),
        minutes=45, hands_on=20,
    )


def after_digest(step: dict, route: dict, bands: list) -> Checkpoint:
    """A digest either went or it did not, and the gel says which."""
    sizes = ", ".join(f"{b:,} bp" for b in bands[:4]) if bands else ""
    expect = []
    if sizes:
        expect.append(f"bands at {sizes}")
    expect.append("no band at the uncut plasmid's size, and no supercoiled "
                  "smear: either means the digest is partial")
    return Checkpoint(
        operation="RE_DIGEST",
        title="Check the digest went to completion",
        how=("Run 2 uL on a 1 % agarose gel beside 100 ng of the uncut "
             "plasmid, so partial digestion is obvious by comparison.",),
        expect=tuple(expect),
        troubles=(
            ("the uncut band is still there",
             "a partial digest, which is the usual cause of a plate full "
             "of re-closed empty vector",
             "digest longer, add enzyme, or check the buffer is the one "
             "both enzymes work in"),
            ("one band where two were expected",
             "only one enzyme cut, or the two sites are too close to "
             "resolve",
             "run each enzyme singly in parallel to find out which"),
            ("a faint ladder of extra bands",
             "star activity from too much enzyme, too much glycerol or "
             "too long an incubation",
             "keep the enzyme under 5 uL in 50, and use the Time-Saver "
             "incubation rather than an overnight"),
        ),
        minutes=40, hands_on=15,
    )


def after_assembly(step: dict, route: dict) -> Checkpoint:
    """The assembly cannot be seen, so the control is the measurement."""
    return Checkpoint(
        operation=step.get("operation", "HI_FI_ASSEMBLY"),
        title="Plate a no-insert control alongside the assembly",
        how=("Set up a second reaction with the vector and water in place "
             "of the insert, and transform and plate it beside the real "
             "one.",
             "An assembly reaction cannot be read on a gel; the plate is "
             "the readout, and a plate on its own says nothing without "
             "the control to compare it against."),
        expect=("many more colonies on the assembly than on the "
                "vector-only plate - tenfold is comfortable, threefold is "
                "workable",
                "a nearly empty vector-only plate"),
        troubles=(
            ("both plates are crowded",
             "the vector is re-closing on itself: an incomplete digest, "
             "or a single-enzyme cut that was never dephosphorylated",
             "go back to the digest gel; dephosphorylate a single-cut "
             "vector"),
            ("both plates are empty",
             "the assembly or the transformation, not the insert",
             "check the competent cells with 100 pg of an intact "
             "plasmid; check the fragments are still on the gel"),
            ("colonies on the assembly plate, but every clone is empty",
             "parental plasmid carried through the PCR",
             "DpnI the amplicon, or gel purify it away from its "
             "template"),
        ),
        minutes=20, hands_on=10,
    )


#: The operations worth stopping after, and what to ask at each.
CHECKED = ("PCR", "RE_DIGEST", "HI_FI_ASSEMBLY", "GIBSON_ASSEMBLY",
           "LIGATION")


def for_route(route: dict, templates: dict | None = None,
              bands: dict | None = None) -> list:
    """A checkpoint for each step of this route that makes a molecule.

    Returns (index of the step it follows, Checkpoint), so a caller can
    splice them in without having to know how the route was built.
    """
    templates = templates or {}
    bands = bands or {}
    steps = route.get("steps") or []
    out, seen_assembly, i = [], False, 0
    while i < len(steps):
        op = steps[i].get("operation")
        if op == "PCR":
            run = [steps[i]]
            while (i + 1 < len(steps)
                   and steps[i + 1].get("operation") == "PCR"):
                i += 1
                run.append(steps[i])
            out.append((i, after_pcrs(run, route, templates)))
        elif op == "RE_DIGEST":
            out.append((i, after_digest(steps[i], route, bands.get(i) or [])))
        elif op in ("HI_FI_ASSEMBLY", "GIBSON_ASSEMBLY", "LIGATION"):
            if not seen_assembly:
                seen_assembly = True
                out.append((i, after_assembly(steps[i], route)))
        i += 1
    return out
