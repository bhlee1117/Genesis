# Genesis

A shared cloning designer for the lab's plasmid inventory.

Genesis indexes the plasmid files on the lab share, normalises them into a
searchable inventory of plasmids and reusable parts, lets you specify a
construct as an ordered list of modules, validates it in silico, and then
generates ranked, bench-ready cloning routes where **every step carries an
explicit laboratory operation code**.

```
Backbone -> promoter -> DIO -> CDS 1 -> trafficking motif -> p2A -> CDS 2 -> polyA
```

goes in; a validated sequence, a ranked set of routes, a primer table and a
stepwise protocol come out.

## Requirements

Python 3.10 or newer. **No third-party packages.** Parsers, thermodynamics,
the web server and the UI are all standard library or vanilla JavaScript, so
this runs on a lab machine with a stock Python install and no network access.

## Quick start

```bash
cd /Users/bhlee1117/Documents/GitHub/Genesis

# 1. Write a config file (edit it afterwards if the share moves)
python3 -m genesis init

# 2. Index the share. Start with a subset to check it looks right.
python3 -m genesis scan --limit 2000

# 3. See what was found
python3 -m genesis stats
python3 -m genesis find CheRiff
python3 -m genesis parts --type promoter

# 4. Plan a construct
python3 -m genesis design examples/example_spec.yaml

# 5. Or use the web interface
python3 -m genesis serve      # http://127.0.0.1:8777
```

## Deployment

The inventory contains the lab's unpublished constructs, with owner names
attached. Genesis is therefore an **internal** tool:

- It binds to `127.0.0.1` by default. Put it on the lab network or behind the
  VPN; do not expose it to the public internet, and do not publish the
  database or a generated site to a public host.
- It **never writes to the DNA share.** Each file is recorded immutably by
  path, checksum, size and mtime; the editable normalised record lives only
  in the local SQLite database.
- Re-running `scan` skips files whose mtime is unchanged, so it is cheap to
  repeat.

## What the importer actually handles

The defaults were set by looking at what the share contains, not by guessing:

| Format | Extensions | Notes |
|---|---|---|
| GenBank | `.gb` `.gbk` `.genbank` | |
| ApE | `.ape` `.str` | GenBank underneath; ~3,600 files on the share |
| SnapGene | `.dna` | binary format read directly, including its feature XML |
| FASTA | `.fa` `.fasta` `.seq` `.txt` | |

Three details that matter in practice:

- **Mixed line endings.** Older files use CR-only or CRLF endings. Everything
  is normalised before parsing; a naive reader sees these files as a single
  line and extracts nothing.
- **Sanger reads are not plasmids.** Most `.seq` files on the share are
  `.ab1` read exports, not plasmid maps. They are indexed and kept, but
  classified as `read` and excluded from the template pool, so a design can
  never be planned off a raw trace.
- **Duplicates collapse.** Sequence checksums are invariant to rotation and
  reverse-complement, so the same plasmid saved by three people with three
  different origins becomes one plasmid with three alternate sources.

Unannotated sequences are marked `annotation_incomplete` rather than being
silently guessed at.

### Annotation

Two independent strategies, because they fail differently:

- **Label rules** map free text (`"CMV promotor"`, `"Amp resistance"`) onto a
  controlled `feature_type` vocabulary, preserving the original label.
- **Motif search** finds elements by sequence. Recombinase sites, ITRs, polyA
  signals and WPRE are matched on DNA; peptide motifs (P2A/T2A, ER2, the
  Kir2.1 trafficking signal, epitope tags) are matched on the **translated
  protein**, so they are found regardless of codon usage.

## Designing a construct

A spec is a short YAML file:

```yaml
name: hSyn-CheRiff-TS-P2A-Citrine-ER2
backbone_plasmid_id: 265
modules:
  - type: promoter
    name: hSyn
    required: true
  - type: CDS
    name: CheRiff
    translated: true
  - type: trafficking_motif
    name: KGC
  - type: p2A
    name: P2A
  - type: CDS
    name: Citrine
  - type: trafficking_motif
    name: ER2
  - type: polyA
    name: bGH polyA
```

Genesis resolves each module against the parts inventory, picks source
plasmids, assembles the intended molecule, and validates it.

### Validation

- **DNA** — module order and orientation, maximal direct repeats (which
  misprime and misassemble), ITR integrity, availability of a unique
  diagnostic enzyme.
- **Protein** — reading frame across the whole fusion and per module,
  premature stops, start/stop codons, 2A peptides flanked by coding sequence
  on both sides, trafficking motifs at the C-terminus. Start and stop codons
  are also marked on the construct map, in construct coordinates, so a
  missing ATG, an in-frame stop inside the fusion, or a fusion that never
  terminates is visible rather than buried in a warning list. Reverse-strand
  runs are mapped back correctly, so an inverted DIO payload shows its ATG at
  the higher coordinate, which is where it really is.
- **AAV packaging** — payload between ITRs against the ~4.7 kb limit.
- **Cre-lox** — an architecture, not a part. A `Cre-lox` module wraps the
  coding sequence that follows it in two heterotypic, antiparallel site
  pairs: `outer - inner - [payload] - outer' - inner'`. The pairs are
  interleaved rather than nested, which is what makes the switch one-way:
  the first pair inverts the payload, which flips the second pair into the
  same orientation so it excises, and the payload cannot flip back.

  Four arrangements, set with `mode:`

  | mode | sites | payload stored | Cre |
  |---|---|---|---|
  | `DIO` | lox2272 + loxP | inverted | switches **on** |
  | `FLEX` | lox2272 + loxP | inverted | switches **on** |
  | `FAS` | loxFAS + loxP | inverted | switches **on** |
  | `DO` | lox2272 + loxP | forward | switches **off** |

  The lab's standard AAV backbone (`PP072` by default, set in `config.py`)
  is used when an AAV is asked for and no backbone is named.

  `type: DIO` still works and means Cre-lox in its DIO arrangement. The
  loxFAS sequence was taken from the lab's own plasmids, where it is
  annotated "FAS SITE"; it carries mutated inverted repeats so it pairs only
  with itself.

  Because the payload is stored reverse-complemented, it is read off the
  reverse strand throughout: translation, reading-frame checks and codon
  markers all follow the strand the fusion is encoded on and map their
  coordinates back, which is why a DIO payload shows its ATG at a *higher*
  coordinate than its stop.
- **Composite parts** — 14% of the CDS annotations on this share have
  another element fused into them: 51 already contain a 2A peptide, and 249
  carry a trafficking motif (Voltron-2, present in 55 plasmids, carries both
  ER2 and the Kir2.1 TS). If you ask for such a CDS *and* the element it
  already contains, Genesis trims the fused copy so the construct does not
  end up with two, and says it did. An element buried inside the coding
  sequence cannot be trimmed without changing the protein, so that is
  reported instead of being removed silently. A tandem 2A pair is an error
  wherever it came from.
- **Recombination** — for DIO/FLEX designs, Cre is simulated and the ORF is
  checked in the **expressed** state. Antiparallel lox pairs invert
  (reversibly); parallel pairs excise (irreversibly, which is what drives the
  switch); only sites with identical spacers recombine. The checked product
  is the one the neuron actually expresses, not the deposited plasmid.

### Vector families

A payload is not a vector. An AAV transfer plasmid needs two ITRs around
the payload, a polyA and normally a WPRE inside them, and the bacterial
origin and resistance outside. A lentiviral one needs LTRs, the packaging
signal, RRE and cPPT instead. Set `vector:` in the specification to `AAV`,
`lentiviral` or `plasmid` and Genesis will only accept a backbone of that
family and will check the finished molecule against the full element list,
reporting anything absent. Asking for a family that the backbone does not belong to reports the
mismatch and names the evidence, rather than quietly leaving the construct
without a backbone. Left blank, the family is inferred from the features the
backbone actually carries, with ITR count outweighing whatever the file
happens to be named. Only a real plasmid's name counts as evidence: the name
someone types for the construct they intend to build is a wish, not an
observation, so a construct called `LENTI-test` is not a lentivirus. A plasmid showing both AAV and lentiviral
evidence is reported as a conflict rather than silently resolved.

### Fused coding sequences

Many CDS annotations on this share are really several proteins joined
together: `LR-Voltron2-LR-CheRiff_eYFP` is Voltron2 (with its own KGC and
ER2 tags), then a P2A, then CheRiff-eYFP. Genesis catalogues the fusion
**and** the proteins it splits into, because both are legitimate answers to
"which part is this", and marks which is which so the picker can offer the
choice.

Boundaries are taken at 2A peptides, IRES elements and at linkers or
trafficking signals that sit between two full-length coding blocks.
Separators that sit next to each other count as one boundary, and where a
cluster contains a 2A that element is the boundary — so Voltron2's
C-terminal KGC and ER2 stay with Voltron2 rather than being mistaken for
part of the junction. A tag on the end of a single protein is not a fusion
and does not split. Of the CDS parts on this share, 100 split this way, and
96 of the resulting pieces turned out to be parts the lab already had names
for.

### Building into the vector

A payload is not bolted onto the end of a backbone. Concatenating gives a
sequence of the right length that is completely wrong: in an AAV it leaves
the whole cassette outside the ITRs, where it cannot be packaged.

Genesis replaces the backbone's own expression cassette instead, keeping
everything that belongs to the vector - both ITRs, the WPRE and polyA that
follow the coding sequence, the origin, the resistance marker - so the
designed molecule is a real plasmid and the map shows all of it rather than
one opaque "backbone" block. Only the promoter-to-CDS stretch is replaced,
and only up to the first WPRE or polyA after it, unless the specification
supplies its own.

The backbone is rotated so the payload region is contiguous before any of
this is worked out, because in a real pAAV the cassette usually wraps the
origin of the coordinate system: in one of this lab's backbones the promoter
sits at 9441 and its own first CDS at 0, so the payload cannot be described
as a single `start..end` slice at all.

### Choosing a recipient and a donor

The plasmid being built into is the first candidate recipient. Proposing a
swap out of one plasmid while assembling into another describes two
different experiments, which is how a plan to build a DIO Voltron construct
came to start from a plasmid holding neither.

What a recipient already has is judged by sequence, not by catalogue
identity. The same protein is catalogued several times over at different
boundaries, so PP072 holding "Voltron-2" does not make it hold the part id
resolved for "LR-Voltron2" - and asking the sequence says it has 100% of it.
A recipient must also be a plasmid that can be grown: without an origin and
a resistance marker it is an extraction or a fragment, not something you can
transform.

Donors are ranked differently from recipients. A recipient is scored on how
much of the target it already holds; a donor only has to hold one piece, so
holding six catalogued variants of GFP does not make a plasmid a better
place to get GFP from. What matters is how little trouble it is: the same
vector family as the construct, sequence verified, and small. Because this
is a lab preference as much as a computation, the alternatives travel with
the plan.



A recipient has to share the construct's *architecture*, not just a few of
its genes. A plasmid carrying three of the right coding sequences but no lox
sites cannot be edited into a Cre-lox construct by swapping one region, and
offering it produces a plan that looks plausible and is nonsense. So a
design using Cre-lox only considers recipients that already carry lox sites,
and an AAV design only considers AAV plasmids. When nothing in the inventory
has both the parts and the architecture, Genesis says so rather than
reaching for the best-scoring wrong answer.

If the recipient already carries a recombinase cassette, the payload goes
*inside* its existing lox sites rather than bringing its own. Four more sites
would give Cre several ways to cut and none of them the intended one.

### Reading the design back

Everything else works forward: modules are resolved, concatenated, and a
layout records where each one went. That bookkeeping is exactly what a bug
would corrupt, so checking a design against it proves nothing.

So every plan is also verified independently. The layout is discarded, the
finished sequence is searched for each part's actual bases, and the molecule
is re-annotated from scratch with the same motif scanner used on imported
files. Only then is what was found compared with what was asked for. A
module that was dropped, truncated, duplicated or placed backwards shows up
here. For a Cre-lox construct the check runs twice, on the plasmid as it
would be cloned and on the product of recombination, because only the second
is expressed.

The test suite injects each of those faults deliberately and asserts the
verifier catches them; a verifier that only ever says "yes" is worth nothing.

### Designs

Saving offers three modes, because they mean different things: `version`
keeps the previous one and adds another, `overwrite` replaces a version in
place, and `save_as` stores under a new name and refuses to merge into an
existing construct's history. Deleting removes one version, or a construct
and all of its versions; the construct record goes when its last version
does.

### Routes

Four routes are generated as first-class alternatives and costed with the
same configurable weights.

**Two-parent swap** — the way most constructs are actually built. Genesis
searches the inventory for the existing plasmid closest to what you want,
works out which single region has to change, and finds a second plasmid to
take the new piece from. It reports the recipient, the exact region to cut
out, the donor, and the exact region to take, bounded by the conserved
features either side (which is what the primers anneal to, or what the
enzymes cut between).

1. `PCR` (inverse, to linearise the recipient and delete the old region) →
   `PCR` (the new piece) → `DPNI_DIGEST` → `DNA_CLEANUP` → `HI_FI_ASSEMBLY`
2. `RE_DIGEST` → `GEL_EXTRACTION` → `LIGATION`

**From scratch** — when no close parent exists, or for comparison:

3. `PCR` → `DPNI_DIGEST` → `DNA_CLEANUP` → `HI_FI_ASSEMBLY`
4. `RE_DIGEST` → `GEL_EXTRACTION` → `LIGATION`

A swap normally wins on cost, because it is two PCRs and one assembly rather
than one fragment per module. If the closest parent already contains
everything you asked for, Genesis says so and recommends checking whether
the construct already exists before cloning anything.

For a restriction swap, the two enzymes must cut **outside** the region being
replaced — one between the upstream conserved feature and the region, one
between the region and the downstream feature. Cutting inside it would leave
a stub of the old sequence behind, which is a silent way to build the wrong
plasmid.

Every step is labelled with an operation code from the library in
`operations.py` (`PCR`, `RE_DIGEST`, `GEL_EXTRACTION`, `DNA_CLEANUP`,
`GIBSON_ASSEMBLY`, `HI_FI_ASSEMBLY`, `LIGATION`, `DEPHOSPHORYLATION`,
`DPNI_DIGEST`, `TRANSFORMATION`, `COLONY_SCREEN`, `SEQUENCE_VERIFY`,
`SYNTHESIS`).

Genesis declines to propose things that fail at the bench:

- A restriction route is only offered when the insert exists as one
  contiguous piece in a single donor, or can be rebuilt in at most a
  two-piece ligation. Beyond that it is marked infeasible with the reason
  stated, rather than quietly suggesting a seven-fragment ligation.
- Amplicons outside the configured PCR window are flagged, including the
  backbone — a 14 kb vector PCR is called out and a restriction
  linearisation suggested instead.
- Fragments under ~100 bp are flagged as better carried on a neighbouring
  primer tail or ordered as annealed oligos than amplified separately.
- Gel extraction is recommended only where extra bands are actually
  expected, not after every PCR.
- Templates and backbones are chosen from maintained circular plasmids in
  preference to PCR-product records and short fragments, which are
  frequently indistinguishable by name alone.

Time is reported two ways, because conflating them misleads: **hands-on**
bench time and **elapsed** turnaround including overnight growth and
sequencing. A plan that is 5 h of work spread over 3 days reads as exactly
that.

### Primers

Every primer is computed from the template sequence. Nothing is reused from
the annotated primer-binding sites in the lab's files.

Design follows Addgene's PCR cloning protocol: the region that hybridises to
the template is **18-21 bp**, and a restriction primer is **leader (3-6 bp) +
recognition site (6-8 bp) + hybridising region**, the leader being there so
the enzyme has something to grip. The leader is not arbitrary - each
candidate is checked against the site it precedes and against every stocked
enzyme, so it cannot accidentally complete a second recognition sequence
where it abuts the site.

Primers are built as `[overlap or leader+site tail] + [annealing region]`,
and **Tm is reported for the annealing region alone** — quoting the Tm of the whole oligo
including its Gibson tail is a common cause of failed first cycles. Melting
temperatures use SantaLucia unified nearest-neighbour parameters with a
divalent-cation correction. Within the 18-21 bp window the length whose Tm
lands closest to the target is chosen; where an AT-rich template cannot
reach the target inside that window, the achievable Tm is reported rather
than the window being widened. Each primer is checked for length, GC, Tm,
homopolymer runs, hairpins, 3' self-dimers, GC-rich 3' ends and mispriming
elsewhere in the template.

### Ranking

Routes are scored with weights in `config.py` — assembly reactions, PCRs,
primers, digests, gel extractions, synthesis, long amplicons, difficult
junctions, unverified templates, non-stocked enzymes — all editable. Recorded
experimental failures feed back into the score, so the ranking can learn what
actually works in this lab.

## Exports

- `pdf` — a printable bench protocol: the cloning scheme drawn as scaled maps
  of both parent plasmids with the cut and the take marked, then reagents,
  the primer table, and numbered steps with reaction set-ups inline. Every
  remaining detail (validation, junctions, sources, route comparison, final
  sequence) is in a section at the end, so the bench pages stay short.
- `protocol` — the same protocol as Markdown
- `primers` — CSV for ordering
- `genbank` — the designed construct, annotated, which round-trips through
  Genesis's own parser

The PDF is written by a small built-in PDF writer (`pdf.py`) using the
base-14 fonts, so exporting needs no third-party package either.

## Layout

```
genesis/
  config.py      Stage 0 deployment + the pre-build checklist
  seqio.py       GenBank / ApE / SnapGene / FASTA readers
  seqtools.py    complement, translation, Tm, ORFs, repeats
  vocab.py       controlled feature vocabulary + motif annotator
  db.py          SQLite schema
  importer.py    Stages 1-2: discovery, normalisation, parts, dedup
  construct.py   Stage 4: spec, resolution, in-silico assembly
  swap.py        closest-parent search and two-parent swap analysis
  protocol.py    bench protocol, rendered to PDF or Markdown
  pdf.py         minimal PDF writer (base-14 fonts, no dependencies)
  validate.py    Stage 5: DNA / protein / Cre recombination checks
  operations.py  Stage 6: laboratory operation library
  enzymes.py     enzyme catalogue + digest simulation
  primers.py     Stage 9: junctions and primer design
  planner.py     Stages 7-8, 10: route generation and ranking
  design.py      orchestration
  export.py      Stage 12: protocol, primer CSV, GenBank
  auth.py        accounts, password hashing, sessions, roles
  vectors.py     AAV / lentiviral families and their element sets
  verify.py      independent read-back of the designed molecule
  server.py      internal web API
web/             single-page UI (vanilla JS, no build step)
tests/           189 tests, no lab share required
```

## Giving it to a labmate

`http://127.0.0.1:8777` means *this computer*. It is not an address other
people can open — on their machine it points at their machine. So the
simplest way to share Genesis is for each person to run their own copy,
which is also the safest: a loopback server is reachable by nobody else, so
it needs no login and no passwords.

Copy the prepared inventory somewhere on the lab share, then each person
runs three commands:

```bash
git clone https://github.com/bhlee1117/Genesis.git
cd Genesis
python3 -m genesis setup --from "/Volumes/cohen_lab/.../genesis.db"
python3 -m genesis serve
```

and opens **http://127.0.0.1:8777**. Copying the database takes seconds and
skips the scan entirely, and everyone sees the same inventory. They do not
need the DNA share mounted unless they want to re-scan it themselves.

To refresh later, re-copy the database, or run `python3 -m genesis scan`.

### If you would rather run one shared server

Set `host` to `0.0.0.0` in `~/.genesis/config.json` and people open
`http://<your-machine's-address>:8777` instead. That address is not
127.0.0.1. Genesis will refuse to start that way until accounts exist,
because an open server on the network exposes every unpublished construct
in the lab to anyone who can reach it. See the next section.

## Accounts and sharing with the lab

Genesis starts local-only and open. Create an account and the login turns on:

```bash
python3 -m genesis useradd bhl --admin     # prompts for a password
python3 -m genesis users                   # list accounts
python3 -m genesis passwd <name>           # set someone's password
python3 -m genesis role <name> --role viewer
python3 -m genesis role <name> --disable
```

Three roles: `viewer` reads and designs, `member` also saves designs and
records outcomes, `admin` also manages accounts and runs imports.

To let the lab reach it, set `host` to `0.0.0.0` in
`~/.genesis/config.json` and give people your machine's address.

### What the login does and does not protect

It stops someone who can reach the server from reading the inventory without
an account. Passwords go through PBKDF2-HMAC-SHA256 with a per-user salt and
600,000 iterations, and only the derived key is stored. Session tokens are
256-bit and held hashed, so a stolen copy of the database yields neither
passwords nor live sessions. Repeated failures are throttled, a wrong
password and an unknown username give the same message, and changing or
disabling an account ends its sessions immediately.

**It does not encrypt anything in transit.** Over plain HTTP on the lab
network, a password and a session cookie are both readable by anything that
can see the traffic. A login is not a substitute for TLS. Either point
`tls_cert` and `tls_key` at a certificate, or put Genesis behind a reverse
proxy that terminates TLS. The sign-in page says so when it is served
unencrypted over a non-loopback address, and the server prints a warning on
startup.

Two configurations are refused outright, because they would put every
unpublished construct in the lab on the network with nothing in front of
them: binding a reachable address with `require_auth` off, and binding a
reachable address before any account exists.

## Permissions

`.claude/settings.json` allows this project's own server command so an agent
working in the repo can restart it without a prompt:

```json
{ "permissions": { "allow": ["Bash(python3 -m genesis serve:*)"] } }
```

It is committed, so it applies for everyone who clones the repo. Move it to
`.claude/settings.local.json` to keep it to one machine.

## Tests

```bash
python3 -m unittest discover -s tests -v
```

The suite uses synthetic sequences only, so it passes with the share
unmounted. It pins down the parts that are easy to get subtly wrong: overhang
polarity and end compatibility, digest geometry on circular molecules,
reading frame across fusions, the Cre invert-then-excise logic, fuzzy motif
search against brute force, and checksum invariance.

## Status and next steps

Implemented: Stages 0-12 of the specification at the "minimal useful first
version" level, plus AAV packaging checks and the SnapGene reader.

Worth adding next:

- Circular plasmid maps in the scheme, as well as the linear ones
- Golden Gate / Type IIS routes (the enzyme catalogue already models them)
- A richer plasmid map viewer with zoom
- Lab-wide authentication and per-user permissions
- Recording outcomes from the UI so route ranking starts learning
