# How Genesis plans a cloning

A specification names an architecture — `pAAV-hSyn-ASAP7y-P2A-EGFP` — and a
protocol comes out the other end. This is what happens in between, in the
order it happens, and why each step is there.

Most of the "why" is a fault that was found by building the thing and
looking at what came out. Those are written down here because the reason a
step exists is the only guard against someone removing it.

The pipeline splits in two. Everything up to step 4 is cheap and decides
**what to build from**; everything after is expensive and depends on that
choice, which is why the Design tab stops in between (`Scan the source`,
then `Generate protocol`).

---

## Phase A — Read the request  (`scan_sources`)

### 1. Resolve each module against the inventory
`construct.resolve`

A module says `CDS: ASAP7y`. This turns it into a concrete sequence plus
the list of plasmids that carry it.

Matching is by **protein as well as DNA**. Two genes can encode the same
protein and share almost no 40-mers, so a codon-optimised copy is invisible
to a DNA search. Protein matching finds 11 GFP parts that DNA matching
misses.

### 2. Read the modules as an ordered architecture
`swap.architecture_candidates`

`hSyn - ASAP7y - P2A - EGFP` is an ordered thing, not a bag of parts. The
useful question is how much of it already exists in one tube.

A recombinase cassette is a **position in the architecture**, not a part:
four lox sites of two kinds. It has no sequence of its own, so it counts
towards the run positionally — the plasmid either has a cassette or does
not — and what follows it starts after the sites rather than straight
after the promoter. Leaving it out made `pAAV-hSyn-DIO-XX` and
`pAAV-hSyn-XX` the same request, so a plasmid that already had the cassette
scored no better than one that would need it built.

### 3. Rank the candidate recipients
same

The architecture is tried **whole first and then shortened from the
right**: `pAAV-PromA-CDS1-P2A-CDS2`, then `pAAV-PromA-CDS1-P2A`, then
`pAAV-PromA-CDS1`, then `pAAV-PromA`. The first length anything matches is
the answer, and `matched` on the result says which length that was.

**Only exact matches count towards the run**, and they have to be in order
and close enough together to be one arrangement — a promoter may be a
Kozak, a 5' UTR or a whole intron away from the coding sequence it drives
(1500 bp), but inside the reading frame there is only ever a linker or a 2A
between one module and the next (150 bp), and anything larger is a
different gene.

Nothing looser counts, because a coding sequence that *resembles* the one
asked for is almost always an **earlier version of the same indicator** —
ASAP4c where ASAP7y was asked for — and handing that to someone is the
specific mistake this is here to prevent. Those plasmids are kept, marked,
and listed separately as `related`; what a plasmid really has (`carries`)
and what it only resembles (`resembles`) are never run together.

Whether a plasmid has a part is asked of its **DNA**, not of the catalogue.
The catalogue links parts through variants and protein matches, which is
what finds a candidate at all and is too loose to decide what must be
built: it called a plasmid holding ASAP4c a plasmid holding ASAP7y, and
simultaneously failed to credit one that really does contain hSyn under
another part row's name.

### 4. Cover the rest with as few donors as possible
same

Longest run first: one donor if some plasmid holds the whole remainder, two
if it only exists in pieces. A build needing `YYY` and `P2A-ZZZ` from
different plasmids is a two-donor build, and saying so beats dropping one.

Donor options are ranked by **carrying the part**, then by covering more of
the architecture, then by being a convenient size — in that order, and the
order matters twice over. Sorting by size alone put six ASAP relatives
above the two plasmids that actually hold ASAP7y. Taking the longest run
first and only then preferring exact matches dropped the exact ones before
they were ever compared: the two plasmids that really hold ASAP7y cover one
module, a plasmid holding an ASAP relative plus a P2A covers two, and the
relative won again.

Relatives are still offered, after the exact ones and labelled with what
they really carry — worth seeing, worth not reaching for by accident.

> **The scan stops here.** Nothing has been designed, and it takes under a
> second. The parents are worth settling before a protocol is written
> around them, because every primer below depends on them.

---

## Phase B — Commit to sources  (`design`)

### 5. Honour an explicit choice
`spec.source_plasmid_ids`

Naming the plasmid a part comes from outranks every ranking above.
Somebody saying "that tube" has a reason the inventory cannot see.

### 6. Take a 2A fusion from a plasmid that already carries it whole
`design._prefer_fused_source`

Asking for `X - P2A - Y` names three modules, and planning them separately
means three reactions and two junctions to design, check and screen. These
constructs are built from each other: the lab has plasmids where `X-P2A-Y`
is already in one piece.

Matching all three sequences end to end is the rare case. Usually the
plasmid joins the same two proteins through **its own** 2A, which is not
what a generic `P2A` resolves to — P2A, T2A and E2A all do the job, and
plain linkers appear too. So the pair is located and whatever sits between
them is read off that plasmid, after checking it could be one: in frame,
45–120 bp, no stop codon inside. The module is rebound to those bases,
which is what makes the three one continuous stretch rather than three that
happen to be adjacent.

### 7. Take anything else the recipient already has from the recipient
`design._prefer_backbone_as_source`

Modules are resolved against the whole catalogue before a backbone is
chosen, so each binds to whichever part the search ranked first and carries
that part's home plasmid with it. Nothing afterwards asked whether the
recipient had the same bases. A pAAV-hSyn vector is chosen *because* it has
hSyn, and the plan amplified hSyn from an unrelated, unverified plasmid
anyway.

Only an exact match counts. A promoter that merely resembles the
recipient's is a different promoter, and replacing it is the point of
naming one.

---

## Phase C — Build the intended molecule

### 8. Expand the recombinase cassette
`construct.expand_dio`

`Cre-lox` denotes an architecture, not a part: it encloses the coding
sequence after it in two heterotypic, antiparallel site pairs. DIO stores
the payload **inverted** so nothing is expressed until Cre flips it.

`expand_dio` hands back modules whose sequence is already the other strand
and merely labels them reverse. Flipping again later counted the inversion
twice and cancelled it, so the payload was amplified the right way up and
assembled into a construct that wanted it upside down — same length, every
part present, different molecule.

If the recipient already carries a cassette, its sites are reused rather
than adding four more; Cre would otherwise have several ways to cut and
none of them the intended one.

### 9. Splice the payload into the recipient
`construct.assemble_into_backbone`

The result is a real plasmid carrying its ITRs, WPRE, polyA, origin and
resistance — not a payload bolted onto the end of one. The backbone is
rotated first so the payload region is contiguous, because a real pAAV
cassette usually wraps the coordinate origin.

Four decisions live here:

- **Two slots, not one.** The payload replaces the recipient's coding
  sequence; the promoter is separate and replaced only when the
  specification names one, or asking for CAG gives you hSyn as well.
- **A lox site never appears in a design that did not ask for one.** Not a
  tidiness rule: Cre turns up in the animal, not in the plasmid, so a
  payload left between somebody else's lox sites is silent until it does,
  and an inverted one stays inverted. Three things enforce it — a
  cassette-bearing plasmid is never chosen automatically for a plain
  design, the sites inside the replaced region go with it, and any that
  survive anyway are a **validation error**, whichever recipient was
  pinned and whichever path put them there.
- **A cassette is only reused when one was asked for.** Dropping a plain
  payload between somebody else's lox sites makes the construct
  Cre-dependent, so nothing is expressed until Cre arrives — a different
  experiment from the one specified. Otherwise the lox sites are removed
  with the payload they enclosed.
- **The region replaced is where the bases actually are.** The catalogued
  hSyn is 485 bp at 144–629 of one recipient, which annotates 174–622 of
  the same promoter. Replacing the *annotated* span while supplying the
  *part* put the bases they disagree about in twice, as a 30 bp direct
  repeat nobody asked for and the read-back could not see.
- **The vector's own WPRE and polyA are reused**, when they survive the
  replacement. Naming one means the construct must have it, not that it
  needs a second copy: two near-identical stretches are a repeat no digest
  can tell apart and no assembly overlap can point at.

The arcs of the recipient that survive are recorded here
(`backbone_segments`), because the route has to amplify those rather than
the whole plasmid.

---

## Phase D — Check the molecule before planning any bench work

### 10. Classify the vector family and check it is complete
`vectors.classify_vector`, `vectors.missing_elements`

An AAV transfer plasmid needs ITRs, WPRE, polyA, origin and resistance. A
payload on its own is not a vector. Only a real plasmid's annotations are
evidence — the name someone typed for the construct they intend to build is
a wish, not an observation, and letting it vote made `LENTI-TEST` detect as
lentiviral on its own.

### 11. Validate
`validate.validate_construct`

Module order, reading frame across every junction, premature stops,
repeated elements, packaging limits.

Order is compared against the **expanded** architecture, since a DIO
expansion legitimately reorders what was typed, and only against modules
this design places — the backbone's own features were never part of the
request.

### 12. Read the finished molecule back
`verify.verify_plan`

Everything above works forward: modules are resolved, concatenated, and the
layout records where each went. That bookkeeping is exactly what a bug
would corrupt, so checking the design against it proves nothing.

This throws the layout away. It takes the finished sequence, searches it
for each part's actual bases, re-annotates it from scratch with the same
motif scanner used on imported files, and only then compares. For a Cre-lox
construct it runs twice: once on the plasmid as cloned and once on the
product of recombination, because those are two different molecules and
only the second is expressed.

---

## Phase E — Plan the bench work

### 13. Group modules into fragments
`planner.group_contiguous`

Neighbouring modules from the same plasmid that sit adjacently *there*
become one PCR. This is where most of the laboratory saving comes from.

Modules are paired with layout entries **by identity, not position**. The
layout also holds the lox sites a cassette adds and the ITRs, origin and
resistance the backbone brings, so the two lists run at different lengths;
zipping them made every fragment after the first inherit a neighbour's
orientation.

### 14. Offer the routes
`swap.analyse_swap`, `planner.plan_*_route`

Six: three two-parent swaps (PCR / digest+PCR / restriction), two
build-from-parts routes (HiFi, restriction), and one that buys the insert
instead of building it.

**Ordering the insert** (`plan_gblock_route`) keeps the longest arc of the
recipient and orders everything else as one fragment with homology arms at
both ends: one vector PCR, one thing to order, one assembly. A Cre-lox
cassette on a recipient that has none is the case it exists for — four 34
bp lox sites are too short to amplify and too short to carry an overlap on
both sides while also being the lox site, so building them means joining
five pieces to make a stretch somebody will sell as one.

Synthesis is costed by the base, because a 2.7 kb fragment is not the same
purchase as a 200 bp one; costing them alike made buying look cheaper than
any number of PCRs. So the PCR route still wins on cost and on time — a
gBlock is real money and about a week of waiting — and past
`synthesis_beats_pcr_at` reactions the plan says so out loud and prints the
sequence to send, rather than quietly reranking.

A single swap puts one piece in one place. Asked for more it used to supply
the first and note the rest, returning a plan whose architecture listed
parts its own protocol never touched — and it won on cost, because doing a
third of the job is cheaper. **A route that cannot build the construct is
not a cheaper way to build it**, so it is now marked unusable.

### 15. Amplify the recipient in one piece
`planner._merge_contiguous`, `planner._carry`

The arcs the recipient keeps and the modules the recipient supplies are the
same DNA lying next to each other. Planned separately they asked for three
PCRs and two junctions to rebuild a stretch nobody had taken apart.
Neighbouring pieces of one plasmid are merged into one reaction, across the
join at position 1 as well.

An arc too short to amplify — the 25 bp between a promoter and the coding
sequence after it — rides on the next piece's forward primer instead, after
the overlap, so the amplicon reads *overlap-arc-body*.

### 16. Design the primers
`primers.design_gibson_pair`, `primers.pick_overlap`

Addgene's PCR cloning rules: 18–21 bp annealing region, Tm computed on the
hybridising region only, 3–6 bp leader before a 6–8 bp restriction site.
NEBuilder HiFi wants 15–25 bp of homology, Gibson 20–40.

Three things that were wrong and are not obvious:

- **One overlap per junction.** Writing it onto both sides made each
  amplicon end with its neighbour's start instead of its own, so the two
  ends meant to anneal shared nothing. The exception is a digested vector,
  which cannot be given a tail, so the insert carries both of its
  junctions.
- **An arm must be unique.** A Cre-lox construct carries each lox site
  twice by design, so an arm taken from one matched both and the payload
  could assemble either way round — in the simulation and in the tube. Arms
  grow past the repeat into sequence that occurs once. Where no arm can be
  unique, the plan says which primer and why.
- **Amplify in the orientation the construct uses**, and decide that by
  which orientation of the donor slice actually equals what is needed. The
  recorded strand is sometimes wrong, and believing it inverts the
  fragment.

Overlaps are taken from the amplicon the reaction will really produce, not
from the catalogued part: the recorded boundaries of an occurrence run a
few bases wide of the part sequence often enough to matter.

### 17. Rank the routes
`planner.rank_routes`

Steps, hands-on time, elapsed time, and a cost model that counts reactions,
primers, difficult junctions, unverified templates and past failures.

---

## Phase F — Run the protocol on paper

### 18. Simulate
`simulate.simulate_route`

Verification (step 12) reads the designed molecule back. That cannot catch
a design which verifies perfectly while the steps underneath it build
something else.

So this ignores the designed sequence entirely:

1. amplify each primer pair off the template it names, trying both strands,
2. join the products on their shared ends, in either orientation,
3. compare the result with the design by **rotation-invariant checksum** —
   same length is weak evidence, the same molecule is the claim,
4. search the product for every requested part, by coverage rather than
   exact containment, since a swap keeps the recipient's own promoter and
   WPRE.

It found four primer faults on its first run, and has caught every
regression since.

### 19. Gate on it
`design`

A route whose steps lose a part is marked unusable, and the routes are
re-sorted afterwards — ranking happens before this, so a route the
simulation has just ruled out would otherwise keep the position it was
ranked into and lead the list while being unusable.

A product that differs from the design without losing a part is a warning,
not a failure: it is the two being built differently, and worth comparing
before building.

---

## Where it still bends

- A 2A fusion is detected for a pair joined by one spacer. A longer run is
  tried as its pairs.
- One recipient with unusual geometry (a 12 kb lentiviral plasmid used as
  an AAV backbone) produces a design whose promoter lands after the
  payload. The simulation reports it as not identical rather than hiding
  it; pinning a normal backbone gives a correct plan.
- Donor choice is offered per segment. There is no way yet to say "take
  this part from here *and* that one from there" in a single saved
  specification beyond `source_plasmid_ids`, which does exactly that but
  has no UI for the general case.
