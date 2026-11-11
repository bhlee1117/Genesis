"""Stage 0 deployment configuration and the Section 16 checklist.

Everything the specification flags as "decide before building" lives here as
an explicit, overridable setting rather than being hard-coded. The shipped
defaults are derived from what the lab share actually contains.
"""

from __future__ import annotations

import json
from .seqtools import CALIBRATED_MG_MM, CALIBRATED_NA_MM
from dataclasses import dataclass, asdict, field
from pathlib import Path

DEFAULT_ROOT = ("/Volumes/cohen_lab/Lab/Data and Analysis/"
                "Electrochromic Protein/E coli stocks")


@dataclass
class Config:
    # --- Stage 0: deployment ---------------------------------------------
    dna_root: str = DEFAULT_ROOT
    db_path: str = "~/.genesis/genesis.db"
    backup_dir: str = "~/.genesis/backups"
    host: str = "127.0.0.1"
    port: int = 8777
    #: Require a signed-in account. Binding anywhere other than loopback is
    #: refused while this is off.
    require_auth: bool = True
    session_hours: float = 12.0
    #: Paths to a TLS certificate and key. Without them the server speaks
    #: plain HTTP, which means passwords cross the network in the clear.
    tls_cert: str = ""
    tls_key: str = ""
    copy_originals: bool = False          # index in place; never write to share
    allow_annotation_edits: bool = True
    design_visibility: str = "lab"        # private | lab | read-only

    # --- Stage 1: discovery ----------------------------------------------
    #: Folder names that are not people.
    non_person_dirs: tuple = (
        "E coli protocols", "Genes from Boyden", "Addgene", "Old", "Archive",
        "Plasmid maps", "Sequencing", "Primers",
    )
    max_depth: int = 6
    skip_dir_patterns: tuple = (".git", "__MACOSX", ".Trash", "node_modules")
    max_file_bytes: int = 20_000_000
    #: Treat short ambiguous FASTA as Sanger reads rather than templates.
    separate_sequencing_reads: bool = True

    # --- Stage 5/7: biological guard rails -------------------------------
    protect_itrs: bool = True             # never PCR across or cut inside ITRs
    max_pcr_length: int = 8000
    min_pcr_length: int = 100
    #: Longest stretch worth ordering as one synthetic fragment. A gBlock
    #: tops out around here; an eBlock goes further but costs more.
    max_synthetic_fragment: int = 3000
    #: Amplifying this many fragments is the point at which ordering the
    #: insert outright is usually the better trade.
    synthesis_beats_pcr_at: int = 3
    allow_unverified_templates: bool = False
    gel_extraction_routine: bool = True

    # --- Stage 6/7: reagents --------------------------------------------
    stocked_enzymes: tuple = (
        "EcoRI", "BamHI", "HindIII", "XbaI", "XhoI", "SalI", "NotI", "AgeI",
        "NheI", "KpnI", "SacI", "PstI", "SpeI", "BglII", "EcoRV", "AscI",
        "BsrGI", "MluI", "NcoI", "SmaI", "ApaI", "ClaI", "BsaI", "BsmBI",
    )
    #: "NEBuilder HiFi" or "Gibson". The two differ in how much homology
    #: they need, which is the main thing it changes about primer design.
    preferred_assembly_kit: str = "NEBuilder HiFi"
    #: Backbone to reach for when an AAV is requested and none is named.
    #: Vector family assumed when a specification does not mention one.
    #: Set to "" to infer it from the backbone instead.
    #: Free token from https://developers.addgene.org/. Only ever used to
    #: search the catalogue for parts the lab has not got; nothing found
    #: there is used as a donor or a PCR template. The environment
    #: variable GENESIS_ADDGENE_TOKEN takes precedence.
    addgene_token: str = ""

    default_vector: str = "AAV"
    default_aav_backbone: str = "PP072"
    default_lenti_backbone: str = ""
    include_transformation_steps: bool = True
    include_sequencing_steps: bool = True

    # --- Stage 9: primer design -----------------------------------------
    primer_tm_target: float = 60.0
    primer_tm_tolerance: float = 3.0
    #: How far the two primers of one reaction may differ in Tm. They
    #: anneal in the same tube at the same temperature, so the pair is
    #: the unit that has to be right, not each primer on its own. NEB's
    #: own calculator refuses a pair past 5 C; this lab runs tighter,
    #: because 5 C is enough for the cooler primer to misprime at a Ta
    #: the hotter one still needs. The junction can usually move far
    #: enough to pay for it.
    primer_pair_max_diff: float = 3.0
    #: How far a junction may move back into its fragment so that the two
    #: primers of one reaction can share an annealing step. The bases it
    #: stops short of are written onto the next fragment's forward tail,
    #: so the product is unchanged and the cost is that one oligo is this
    #: much longer. Set to 0 to pin every junction where the construct
    #: puts it and take the Tm mismatch instead.
    primer_max_retreat: int = 25
    #: Put a unique restriction site at each end of the cassette - in
    #: front of the Kozak and past the stop codon - so the insert can be
    #: swapped later without redesigning the whole plasmid. The sites are
    #: chosen against the finished sequence, so they cut once or not at
    #: all.
    #:
    #: Off until the planner stops tripping over the extra module. The
    #: handles themselves are right - unique, popular enzymes, no start
    #: codon introduced in front of the Kozak, product simulates
    #: identical - but a twelve-base module is below the length worth
    #: amplifying, so it is folded onto the primer tails that meet there,
    #: and a folded junction cannot move. That costs the retreat that
    #: keeps a primer pair inside its Tm limit, and three planner tests
    #: fail outright.
    junction_sites: bool = False
    #: Addgene's PCR cloning protocol: the region that hybridises to the
    #: template is 18-21 bp, which is "usually sufficient to give specificity
    #: and to also be compatible with a standard PCR reaction".
    #: https://www.addgene.org/protocols/pcr-cloning/
    primer_min_len: int = 18
    primer_max_len: int = 21
    #: Extra bases 5' of a restriction site, so the enzyme has something to
    #: hold on to. Addgene: "usually 3-6bp".
    primer_leader: str = "TTAATT"
    primer_max_gc: float = 0.65
    primer_min_gc: float = 0.35
    #: The band NEB asks for in the Q5 protocol - "ideally have a GC
    #: content of 40-60%". The wider pair above is where Genesis starts
    #: warning; this narrower one is what it aims at while choosing where
    #: to put the primer.
    primer_ideal_min_gc: float = 0.40
    primer_ideal_max_gc: float = 0.60
    #: Homology arm. Used as a ceiling: the shortest arm that still reaches
    #: `overlap_tm_min` is taken, because every base of overlap is a base of
    #: oligo somebody pays for and synthesises with some error rate.
    overlap_length: int = 25
    #: A floor under whatever the kit itself asks for, not a replacement
    #: for it: homology_range() already holds NEBuilder HiFi at 15-25 and
    #: Gibson at 20-40, which is what each supplier specifies.
    #:
    #: Raised to 20 for this lab. NEB's own range starts at 15, but 15 is
    #: the point at which it still works rather than the point at which
    #: it works well, and NEB themselves ask for more as soon as there
    #: are several fragments in the tube. Five more bases on one oligo is
    #: a cheaper way to buy an assembly than a repeat of the assembly.
    overlap_min_length: int = 20
    overlap_tm_min: float = 48.0
    primer_conc_nm: float = 500.0
    #: Free-ion concentrations of the reaction the Tm is quoted for.
    #: 2 mM Mg2+ is what Q5 buffer contains; the monovalent figure was
    #: fitted to 84 values from NEB's own Tm API across 18-40 nt and
    #: 20-80% GC, because NEB's calculator is what a primer gets checked
    #: against before it is ordered and disagreeing with it by 5 C is a
    #: bug however defensible the thermodynamics.
    na_mm: float = CALIBRATED_NA_MM
    mg_mm: float = CALIBRATED_MG_MM

    # --- Stage 10: route ranking weights (higher = more discouraged) -----
    cost_weights: dict = field(default_factory=lambda: {
        "assembly_reaction": 10.0,
        "pcr": 6.0,
        "primer": 1.5,
        "restriction_digest": 3.0,
        "gel_extraction": 4.0,
        "synthesis": 40.0,
        # Per kb of ordered DNA, on top of the per-fragment fee. Roughly
        # what a gBlock costs relative to a reaction somebody runs.
        "synthesis_per_kb": 90.0,
        "long_pcr_per_kb": 1.0,
        "difficult_junction": 5.0,
        "unverified_template": 8.0,
        "non_stocked_enzyme": 12.0,
        "historical_failure": 15.0,
    })

    def save(self, path: str | Path) -> None:
        p = Path(path).expanduser()
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(asdict(self), indent=2, default=list))

    @classmethod
    def load(cls, path: str | Path | None = None) -> "Config":
        if path is None:
            path = Path("~/.genesis/config.json").expanduser()
        p = Path(path).expanduser()
        if not p.exists():
            return cls()
        data = json.loads(p.read_text())
        known = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in data.items() if k in known})

    @property
    def resolved_db(self) -> Path:
        return Path(self.db_path).expanduser()
