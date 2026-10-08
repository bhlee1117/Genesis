"""Stage 0 deployment configuration and the Section 16 checklist.

Everything the specification flags as "decide before building" lives here as
an explicit, overridable setting rather than being hard-coded. The shipped
defaults are derived from what the lab share actually contains.
"""

from __future__ import annotations

import json
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
    default_aav_backbone: str = "PP072"
    default_lenti_backbone: str = ""
    include_transformation_steps: bool = True
    include_sequencing_steps: bool = True

    # --- Stage 9: primer design -----------------------------------------
    primer_tm_target: float = 60.0
    primer_tm_tolerance: float = 3.0
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
    #: Homology arm. Used as a ceiling: the shortest arm that still reaches
    #: `overlap_tm_min` is taken, because every base of overlap is a base of
    #: oligo somebody pays for and synthesises with some error rate.
    overlap_length: int = 25
    overlap_min_length: int = 15
    overlap_tm_min: float = 48.0
    primer_conc_nm: float = 500.0
    na_mm: float = 50.0
    mg_mm: float = 2.0

    # --- Stage 10: route ranking weights (higher = more discouraged) -----
    cost_weights: dict = field(default_factory=lambda: {
        "assembly_reaction": 10.0,
        "pcr": 6.0,
        "primer": 1.5,
        "restriction_digest": 3.0,
        "gel_extraction": 4.0,
        "synthesis": 40.0,
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
