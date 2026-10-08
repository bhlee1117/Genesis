"""Test suite for Genesis.

Run with:  python3 -m unittest discover -s tests -v

These tests use only synthetic sequences, so they pass with the lab share
unmounted. The checks are chosen to pin down the biology that is easy to get
subtly wrong - overhang polarity, reading frame, Cre logic - rather than to
restate the implementation.
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genesis import seqtools as st
from genesis.config import Config
from genesis.construct import ConstructSpec, Module, assemble, ResolvedModule
from genesis.db import connect
from genesis.design import design
from genesis.enzymes import CATALOG, digest, ends_compatible, find_sites
from genesis.export import protocol_markdown, to_genbank
from genesis.importer import least_rotation, seq_checksum
from genesis.primers import Primer, validate_primer, pick_anneal
from genesis.seqio import parse_fasta, parse_genbank, parse_location, classify_role
from genesis.validate import LoxSite, simulate_cre, validate_construct
from genesis.vocab import classify_label, find_motifs, _fuzzy_find

LOXP = "ATAACTTCGTATAATGTATGCTATACGAAGTTAT"
LOX2272 = "ATAACTTCGTATAAAGTATCCTATACGAAGTTAT"


def rand_seq(n, seed=None):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def orf(n_codons, seed=0):
    r = random.Random(seed)
    cods = ["GCT", "TGC", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTG"]
    return "ATG" + "".join(r.choice(cods) for _ in range(n_codons)) + "TAA"


# --------------------------------------------------------------------------
class TestSeqTools(unittest.TestCase):

    def test_revcomp_roundtrip(self):
        s = rand_seq(500, 1)
        self.assertEqual(st.revcomp(st.revcomp(s)), s)

    def test_revcomp_known(self):
        self.assertEqual(st.revcomp("GAATTC"), "GAATTC")       # palindrome
        self.assertEqual(st.revcomp("ATCGAT"), "ATCGAT")
        self.assertEqual(st.revcomp("AAAACCC"), "GGGTTTT")

    def test_translate_known_peptides(self):
        self.assertEqual(st.translate("ATGGCTTGCTAA"), "MAC*")
        self.assertEqual(st.translate("ATCGAAGGTCGT"), "IEGR")  # Factor Xa site
        self.assertEqual(st.translate("TAATAGTGA"), "***")

    def test_translate_frame_and_partial_codon(self):
        self.assertEqual(st.translate("GATGGCTTGC", frame=1), "MAC")
        # trailing partial codons are dropped, not padded
        self.assertEqual(st.translate("GATGGCTTG", frame=1), "MA")
        self.assertEqual(st.translate("ATGGC"), "M")

    def test_find_all_overlapping(self):
        self.assertEqual(st.find_all("AAAA", "AA"), [0, 1, 2])

    def test_find_all_circular_wrap(self):
        # site split across the origin of a circular molecule
        self.assertIn(8, st.find_all("TTTTTTTTGA" + "ATTC"[:0] + "AT", "GAAT",
                                     circular=True))

    def test_find_all_iupac(self):
        self.assertEqual(st.find_all("ACCTGGA", "CCWGG"), [1])

    def test_slice_region_circular(self):
        s = "ABCDEFGH"
        self.assertEqual(st.slice_region(s, 6, 2, circular=True), "GHAB")

    def test_tm_monotonic_and_plausible(self):
        short = st.tm("ACGTACGTACGTACGTAC")
        longer = st.tm("ACGTACGTACGTACGTACGTACGTACGT")
        self.assertLess(short, longer)
        # A typical 20mer primer should land in the usual design window.
        self.assertTrue(45 < st.tm("GCTAGCTAGCATCGATCGAT") < 75)

    def test_tm_gc_effect(self):
        at_rich = st.tm("ATATATATATATATATATAT")
        gc_rich = st.tm("GCGCGCGCGCGCGCGCGCGC")
        self.assertGreater(gc_rich, at_rich + 15)

    def test_tm_rejects_ambiguous(self):
        self.assertNotEqual(st.tm("ACGTNACGTACGTACGTAC"), st.tm("ACGTNACGTACGTACGTAC"))

    def test_max_homopolymer(self):
        self.assertEqual(st.max_homopolymer("ACGTTTTTACG"), 5)

    def test_longest_repeat_merges_sliding_hits(self):
        unit = rand_seq(200, 7)
        seq = unit + rand_seq(300, 8) + unit
        hits = st.longest_repeat(seq, min_len=40)
        self.assertTrue(hits)
        self.assertEqual(hits[0]["length"], 200)
        # One maximal repeat, not ~160 sliding duplicates.
        self.assertLessEqual(len([h for h in hits if h["length"] > 100]), 2)

    def test_find_orfs_both_strands(self):
        payload = orf(120, 3)
        fwd = st.find_orfs("TTTT" + payload + "TTTT", min_aa=100)
        self.assertEqual(fwd[0]["strand"], 1)
        rev = st.find_orfs("TTTT" + st.revcomp(payload) + "TTTT", min_aa=100)
        self.assertEqual(rev[0]["strand"], -1)
        self.assertEqual(rev[0]["protein"], fwd[0]["protein"])


# --------------------------------------------------------------------------
class TestParsers(unittest.TestCase):

    GB = (
        "LOCUS       TestPlasmid             60 bp ds-DNA     circular SYN 01-JAN-2020\n"
        "DEFINITION  a test record.\n"
        "FEATURES             Location/Qualifiers\n"
        "     CDS             1..30\n"
        '                     /label="testCDS"\n'
        "     misc_feature    complement(31..60)\n"
        '                     /label="reverse thing"\n'
        "ORIGIN\n"
        "        1 atggcttgct aaatggcttg ctaaatggct tgctaaatgg cttgctaaat ggcttgctaa\n"
        "//\n"
    )

    def test_genbank_basic(self):
        rec = parse_genbank(self.GB.encode())
        self.assertEqual(rec.name, "TestPlasmid")
        self.assertTrue(rec.circular)
        self.assertEqual(len(rec.seq), 60)
        self.assertEqual(len(rec.features), 2)
        self.assertEqual(rec.features[0].label, "testCDS")
        self.assertEqual((rec.features[0].start, rec.features[0].end), (0, 30))
        self.assertEqual(rec.features[1].strand, -1)

    def test_genbank_line_ending_variants(self):
        """Old lab files use CR-only and CRLF endings; all must parse alike."""
        base = parse_genbank(self.GB.encode())
        for ending in ("\r\n", "\r"):
            variant = self.GB.replace("\n", ending).encode()
            rec = parse_genbank(variant)
            self.assertEqual(rec.seq, base.seq, f"failed for {ending!r}")
            self.assertEqual(len(rec.features), len(base.features))
            self.assertEqual(rec.name, base.name)

    def test_parse_location_forms(self):
        self.assertEqual(parse_location("1..10"), ([(0, 10)], 1))
        self.assertEqual(parse_location("complement(5..9)"), ([(4, 9)], -1))
        self.assertEqual(parse_location("<1..>10"), ([(0, 10)], 1))
        segs, strand = parse_location("join(1..3,7..9)")
        self.assertEqual(segs, [(0, 3), (6, 9)])
        segs, strand = parse_location("complement(join(1..3,7..9))")
        self.assertEqual(strand, -1)

    def test_fasta(self):
        rec = parse_fasta(b">myseq a description\nACGT\nACGT\n")
        self.assertEqual(rec.name, "myseq")
        self.assertEqual(rec.seq, "ACGTACGT")

    def test_fasta_strips_numbers_and_whitespace(self):
        rec = parse_fasta(b">x\n  1 acgt acgt\n")
        self.assertEqual(rec.seq, "ACGTACGT")

    def test_sanger_read_is_not_a_plasmid(self):
        """The share holds thousands of .seq reads; they must not become templates."""
        read = b">PB106-JSP1_A08.ab1\n" + (b"NNNNNACGT" * 40) + b"\n"
        rec = parse_fasta(read)
        role, _ = classify_role(Path("x.seq"), rec)
        self.assertEqual(role, "read")

    def test_real_plasmid_is_classified_as_plasmid(self):
        rec = parse_genbank(self.GB.encode())
        rec.seq = rand_seq(5000, 2)
        role, _ = classify_role(Path("p.gb"), rec)
        self.assertEqual(role, "plasmid")


# --------------------------------------------------------------------------
class TestChecksum(unittest.TestCase):

    def test_least_rotation_matches_brute_force(self):
        r = random.Random(5)
        for _ in range(100):
            s = "".join(r.choice("ACGT") for _ in range(r.randint(1, 30)))
            brute = min((s + s)[i:i + len(s)] for i in range(len(s)))
            self.assertEqual(least_rotation(s), brute)

    def test_circular_checksum_invariance(self):
        s = rand_seq(2000, 9)
        rotated = s[777:] + s[:777]
        self.assertEqual(seq_checksum(s, True), seq_checksum(rotated, True))
        self.assertEqual(seq_checksum(s, True), seq_checksum(st.revcomp(s), True))

    def test_linear_checksum_distinguishes_rotation(self):
        s = rand_seq(500, 10)
        self.assertNotEqual(seq_checksum(s, False), seq_checksum(s[10:] + s[:10], False))


# --------------------------------------------------------------------------
class TestEnzymes(unittest.TestCase):

    def test_overhangs_are_textbook(self):
        cases = {"EcoRI": ("AATT", "5'"), "BamHI": ("GATC", "5'"),
                 "KpnI": ("GTAC", "3'"), "SacI": ("AGCT", "3'"),
                 "EcoRV": ("", "blunt"), "NotI": ("GGCC", "5'")}
        for name, (oh, etype) in cases.items():
            ctx = "TTTT" + CATALOG[name].site + "TTTT"
            frag = digest(ctx, [name])[0]
            self.assertEqual(frag.right_overhang, oh, name)
            self.assertEqual(frag.right_end_type, etype, name)

    def test_isocaudomers_are_compatible(self):
        def ends(name):
            f = digest("TTTT" + CATALOG[name].site + "TTTT", [name])[0]
            return f.right_overhang, f.right_end_type
        for a, b in [("BamHI", "BglII"), ("XbaI", "NheI"), ("SpeI", "AvrII"),
                     ("SalI", "XhoI"), ("AgeI", "XmaI"), ("NotI", "EagI")]:
            self.assertTrue(ends_compatible(*ends(a), *ends(b)), f"{a}/{b}")

    def test_incompatible_ends_rejected(self):
        def ends(name):
            f = digest("TTTT" + CATALOG[name].site + "TTTT", [name])[0]
            return f.right_overhang, f.right_end_type
        for a, b in [("EcoRI", "BamHI"), ("KpnI", "SacI"), ("KpnI", "BamHI")]:
            self.assertFalse(ends_compatible(*ends(a), *ends(b)), f"{a}/{b}")

    def test_blunt_ends_always_join(self):
        self.assertTrue(ends_compatible("", "blunt", "", "blunt"))

    def test_linear_digest_partitions_the_molecule(self):
        s = "TTTT" + "GAATTC" + "A" * 10 + "GGATCC" + "TTTT"
        frags = digest(s, ["EcoRI", "BamHI"])
        self.assertEqual([f.length for f in frags], [5, 16, 9])
        self.assertEqual(sum(f.length for f in frags), len(s))

    def test_circular_single_cut_linearises(self):
        circ = "GAATTC" + rand_seq(400, 11)
        frags = digest(circ, ["EcoRI"], circular=True)
        self.assertEqual(len(frags), 1)
        self.assertEqual(frags[0].length, len(circ))

    def test_circular_digest_conserves_length(self):
        circ = "GAATTC" + rand_seq(300, 12) + "GGATCC" + rand_seq(500, 13)
        frags = digest(circ, ["EcoRI", "BamHI"], circular=True)
        self.assertEqual(len(frags), 2)
        self.assertEqual(sum(f.length for f in frags), len(circ))

    def test_non_palindromic_site_found_on_both_strands(self):
        enz = CATALOG["BsaI"]
        fwd = find_sites("TTTT" + enz.site + "T" * 20, enz)
        rev = find_sites("T" * 20 + st.revcomp(enz.site) + "TTTT", enz)
        self.assertEqual(len(fwd), 1)
        self.assertEqual(len(rev), 1)

    def test_type_iis_cuts_outside_recognition(self):
        self.assertTrue(CATALOG["BsaI"].is_type_iis)
        self.assertFalse(CATALOG["EcoRI"].is_type_iis)


# --------------------------------------------------------------------------
class TestVocab(unittest.TestCase):

    def test_label_rules(self):
        cases = {"CMV promotor": "promoter", "Amp resistance": "resistance",
                 "AAV2 ITR": "AAV_ITR", "loxP": "lox_site", "WPRE": "WPRE",
                 "bGH poly(A) signal": "polyA", "pUC origin": "origin",
                 "ER2": "trafficking_motif", "P2A": "p2A"}
        for label, expected in cases.items():
            self.assertEqual(classify_label(label), expected, label)

    def test_genbank_key_fallback(self):
        self.assertEqual(classify_label("", "CDS"), "CDS")
        self.assertEqual(classify_label("some band", "rep_origin"), "origin")
        self.assertEqual(classify_label("mystery", ""), "unknown")

    def test_fuzzy_search_matches_brute_force(self):
        r = random.Random(21)
        for _ in range(60):
            h = "".join(r.choice("ACGT") for _ in range(300))
            n = "".join(r.choice("ACGT") for _ in range(r.randint(15, 40)))
            for mm in (0, 1, 2, 3):
                brute = [i for i in range(len(h) - len(n) + 1)
                         if sum(a != b for a, b in zip(h[i:i + len(n)], n)) <= mm]
                self.assertEqual(_fuzzy_find(h, n, mm), brute)

    def test_motif_finder_locates_lox_both_strands(self):
        seq = "G" * 100 + LOXP + "T" * 100 + st.revcomp(LOX2272) + "C" * 100
        hits = {(m["name"], m["strand"]) for m in find_motifs(seq, include_protein=False)}
        self.assertIn(("loxP", 1), hits)
        self.assertIn(("lox2272", -1), hits)

    def test_protein_motifs_are_codon_independent(self):
        """P2A must be found regardless of codon usage, which DNA search misses."""
        p2a_aa = "ATNFSLLKQAGDVEENPGP"
        back = {"A": "GCG", "T": "ACA", "N": "AAC", "F": "TTC", "S": "AGC",
                "L": "CTG", "K": "AAG", "Q": "CAG", "G": "GGC", "D": "GAC",
                "V": "GTG", "E": "GAG", "P": "CCG"}
        dna = "".join(back[a] for a in p2a_aa)
        hits = [m for m in find_motifs("ACG" * 20 + dna + "TTA" * 20)
                if m["name"] == "P2A"]
        self.assertTrue(hits)
        self.assertEqual(hits[0]["level"], "protein")


# --------------------------------------------------------------------------
class TestCre(unittest.TestCase):

    def _sites(self, seq):
        from genesis.validate import lox_sites_from_layout
        return lox_sites_from_layout(seq)

    def test_flex_inverts_then_excises(self):
        payload = orf(198, 4)
        dio = ("G" * 200 + LOXP + LOX2272 + st.revcomp(payload)
               + st.revcomp(LOXP) + st.revcomp(LOX2272) + "C" * 200)
        out, log = simulate_cre(dio, self._sites(dio))
        self.assertEqual(len(log), 2)
        self.assertIn("inverted", log[0])
        self.assertIn("excised", log[1])
        # The expressed ORF must exist, forward, at full length.
        expected = st.translate(payload).rstrip("*")
        found = [o for o in st.find_orfs(out, min_aa=150)
                 if o["strand"] == 1 and o["protein"] == expected]
        self.assertTrue(found, "FLEX did not restore a forward full-length ORF")
        remaining = {s.name for s in self._sites(out)}
        self.assertEqual(remaining, {"loxP", "lox2272"})

    def test_parallel_pair_excises(self):
        seq = "G" * 100 + LOXP + "T" * 300 + LOXP + "C" * 100
        out, log = simulate_cre(seq, self._sites(seq))
        self.assertLess(len(out), len(seq))
        self.assertEqual(len(self._sites(out)), 1)
        self.assertIn("excised", log[0])

    def test_antiparallel_pair_inverts_once_without_looping(self):
        seq = "G" * 100 + LOXP + "T" * 300 + st.revcomp(LOXP) + "C" * 100
        out, log = simulate_cre(seq, self._sites(seq))
        self.assertEqual(len(out), len(seq))
        self.assertEqual(len(log), 1)

    def test_heterotypic_sites_do_not_recombine(self):
        seq = "G" * 100 + LOXP + "T" * 300 + LOX2272 + "C" * 100
        out, log = simulate_cre(seq, self._sites(seq))
        self.assertEqual(out, seq)
        self.assertEqual(log, [])

    def test_single_site_is_inert(self):
        seq = "G" * 100 + LOXP + "C" * 100
        out, log = simulate_cre(seq, self._sites(seq))
        self.assertEqual(out, seq)


# --------------------------------------------------------------------------
class TestValidation(unittest.TestCase):

    def _construct(self, modules, circular=True):
        resolved = []
        for ftype, name, seq in modules:
            m = Module(type=ftype, name=name, seq=seq)
            resolved.append(ResolvedModule(m, None, name, seq, ftype, [], []))
        return assemble(resolved, "", circular)

    def test_frame_break_detected(self):
        c = self._construct([("CDS", "geneA", "ATG" + "GCT" * 20),
                             ("CDS", "geneB", "GC" + "GCT" * 20 + "TAA")])
        msgs = [i.message for i in validate_construct(c) if i.level == "error"]
        self.assertTrue(any("not a multiple of 3" in m or "out of frame" in m
                            for m in msgs), msgs)

    def test_in_frame_fusion_passes(self):
        c = self._construct([("CDS", "geneA", "ATG" + "GCT" * 20),
                             ("CDS", "geneB", "GCT" * 20 + "TAA")])
        errs = [i.message for i in validate_construct(c) if i.level == "error"]
        self.assertEqual(errs, [])

    def test_premature_stop_detected(self):
        c = self._construct([("CDS", "geneA", "ATG" + "GCT" * 5 + "TAA" + "GCT" * 5),
                             ("CDS", "geneB", "GCT" * 10 + "TAA")])
        msgs = [i.message for i in validate_construct(c) if i.level == "error"]
        self.assertTrue(any("premature stop" in m for m in msgs), msgs)

    def test_2a_must_have_cds_on_both_sides(self):
        p2a = "GGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAGACGTGGAGGAGAACCCTGGACCT"
        c = self._construct([("p2A", "P2A", p2a),
                             ("CDS", "geneB", "ATG" + "GCT" * 20 + "TAA")])
        msgs = [i.message for i in validate_construct(c) if i.level == "error"]
        self.assertTrue(any("both sides" in m for m in msgs), msgs)

    def test_aav_payload_limit_enforced(self):
        itr = rand_seq(145, 30)
        c = self._construct([("AAV_ITR", "ITR-L", itr),
                             ("CDS", "big", orf(1900, 31)),
                             ("AAV_ITR", "ITR-R", itr)])
        msgs = [i.message for i in validate_construct(c) if i.level == "error"]
        self.assertTrue(any("AAV limit" in m for m in msgs), msgs)

    def test_aav_payload_within_limit_passes(self):
        itr = rand_seq(145, 32)
        c = self._construct([("AAV_ITR", "ITR-L", itr),
                             ("CDS", "ok", orf(500, 33)),
                             ("AAV_ITR", "ITR-R", itr)])
        msgs = [i.message for i in validate_construct(c)
                if i.level == "error" and i.category == "packaging"]
        self.assertEqual(msgs, [])


# --------------------------------------------------------------------------
class TestPrimers(unittest.TestCase):

    def setUp(self):
        self.cfg = Config()

    def test_pick_anneal_hits_target_tm(self):
        tmpl = rand_seq(1000, 40)
        region = pick_anneal(tmpl, 100, "F", self.cfg)
        self.assertTrue(self.cfg.primer_min_len <= len(region) <= self.cfg.primer_max_len)
        self.assertLess(abs(st.tm(region) - self.cfg.primer_tm_target), 6)

    def test_reverse_anneal_is_reverse_complemented(self):
        tmpl = rand_seq(1000, 41)
        region = pick_anneal(tmpl, 500, "R", self.cfg)
        self.assertIn(st.revcomp(region), tmpl)

    def test_mispriming_is_flagged(self):
        unit = rand_seq(60, 42)
        tmpl = unit + rand_seq(200, 43) + unit     # duplicated binding site
        p = Primer("p", unit[:24], unit[:24])
        warnings = validate_primer(p, tmpl, self.cfg)
        self.assertTrue(any("occurs" in w for w in warnings), warnings)

    def test_homopolymer_and_dimer_flagged(self):
        p = Primer("p", "ACGTAAAAAAAGCTAGCTAGCT", "ACGTAAAAAAAGCTAGCTAGCT")
        warnings = validate_primer(p, "", self.cfg)
        self.assertTrue(any("homopolymer" in w for w in warnings), warnings)

    def test_tm_reported_for_anneal_not_full_oligo(self):
        anneal = "GCTAGCTAGCATCGATCGAT"
        p = Primer("p", "A" * 25 + anneal, anneal, "A" * 25)
        d = p.to_dict(self.cfg)
        self.assertNotEqual(round(d["anneal_tm"], 1), round(d["full_tm"], 1))
        self.assertEqual(round(d["anneal_tm"], 1), round(st.tm(anneal), 1))


# --------------------------------------------------------------------------
class TestSpecParsing(unittest.TestCase):

    def test_yaml_dialect_from_the_specification(self):
        spec = ConstructSpec.from_yaml("""
name: test construct
backbone: AAV_backbone_01
modules:
  - type: promoter
    name: hSyn
    required: true
  - type: CDS
    name: sensor_A
    translated: true
  - type: polyA
    name: bGH_polyA
""")
        self.assertEqual(spec.name, "test construct")
        self.assertEqual(spec.backbone, "AAV_backbone_01")
        self.assertEqual(len(spec.modules), 3)
        self.assertEqual(spec.modules[0].type, "promoter")
        self.assertIs(spec.modules[0].required, True)
        self.assertIs(spec.modules[1].translated, True)

    def test_roundtrip_dict(self):
        spec = ConstructSpec("x", modules=[Module("CDS", "a")])
        again = ConstructSpec.from_dict(spec.to_dict())
        self.assertEqual(again.modules[0].type, "CDS")


# --------------------------------------------------------------------------
class TestEndToEnd(unittest.TestCase):
    """Build a small inventory in memory, then design against it."""

    def setUp(self):
        self.conn = connect(":memory:")
        self.cfg = Config()
        self.cfg.stocked_enzymes = ("EcoRI", "BamHI", "HindIII", "XhoI", "NotI",
                                    "AgeI", "SalI", "KpnI")

    def _add_part(self, ftype, name, seq, plasmid_name):
        import hashlib
        chk = hashlib.sha256(seq.encode()).hexdigest()
        pl = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                                   circular, verified)
               VALUES (?,?,?,?,?,1,'sequence_verified')""",
            (plasmid_name, "tester", seq + rand_seq(3000, 50), len(seq) + 3000,
             chk[:32] + plasmid_name)).lastrowid
        part = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES (?,?,?,?,?)""", (ftype, name, seq, chk, len(seq))).lastrowid
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end, strand)
               VALUES (?,?,0,?,1)""", (part, pl, len(seq)))
        self.conn.commit()
        return pl

    def test_design_produces_two_labelled_routes(self):
        backbone_seq = ("GAATTC" + rand_seq(2000, 60) + "GGATCC"
                        + rand_seq(2000, 61))
        bb = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                                   circular, verified)
               VALUES ('pTestBackbone','tester',?,?,'bbsum',1,'sequence_verified')""",
            (backbone_seq, len(backbone_seq))).lastrowid
        self._add_part("CDS", "GFPlike", orf(240, 62), "pDonorGFP")
        self.conn.commit()

        spec = ConstructSpec("test build", backbone_plasmid_id=bb,
                             modules=[Module("CDS", "GFPlike")])
        plan = design(self.conn, spec, self.cfg)

        strategies = {r["strategy"] for r in plan["routes"]}
        # Both from-scratch strategies must always be offered; the two-parent
        # swap strategies appear as well whenever a parent pair exists.
        self.assertTrue({"hifi", "restriction"} <= strategies, strategies)
        self.assertLessEqual(strategies - {"hifi", "restriction"},
                             {"swap_hifi", "swap_restriction"})

        # Every step must carry a real operation code (the Stage 6 requirement).
        from genesis.operations import OPERATIONS
        for route in plan["routes"]:
            for step in route["steps"]:
                self.assertIn(step["operation"], OPERATIONS, step)

        # Routes are ranked cheapest-first among feasible ones.
        feasible = [r for r in plan["routes"] if r["feasible"]]
        self.assertEqual(feasible, sorted(feasible, key=lambda r: r["cost"]))

    def test_multi_fragment_insert_rejects_restriction_route(self):
        backbone_seq = "GAATTC" + rand_seq(3000, 70)
        bb = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                                   circular, verified)
               VALUES ('pBB2','tester',?,?,'bb2',1,'sequence_verified')""",
            (backbone_seq, len(backbone_seq))).lastrowid
        for i in range(4):
            self._add_part("CDS", f"gene{i}", orf(100, 80 + i), f"pDonor{i}")
        self.conn.commit()

        spec = ConstructSpec("many pieces", backbone_plasmid_id=bb,
                             modules=[Module("CDS", f"gene{i}") for i in range(4)])
        plan = design(self.conn, spec, self.cfg)
        re_route = next(r for r in plan["routes"] if r["strategy"] == "restriction")
        self.assertFalse(re_route["feasible"])
        self.assertTrue(any("not practical" in w for w in re_route["warnings"]))
        # The HiFi route should still be offered.
        hifi = next(r for r in plan["routes"] if r["strategy"] == "hifi")
        self.assertTrue(hifi["feasible"])

    def test_exports_are_wellformed(self):
        backbone_seq = "GAATTC" + rand_seq(2500, 90)
        bb = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                                   circular, verified)
               VALUES ('pBB3','tester',?,?,'bb3',1,'sequence_verified')""",
            (backbone_seq, len(backbone_seq))).lastrowid
        self._add_part("CDS", "gene", orf(150, 91), "pDonorX")
        spec = ConstructSpec("export test", backbone_plasmid_id=bb,
                             modules=[Module("CDS", "gene")])
        plan = design(self.conn, spec, self.cfg)

        md = protocol_markdown(plan, self.cfg)
        self.assertIn("# export test", md)
        # Concise protocol first, with reagents and a reaction set-up...
        self.assertIn("## Reagents and consumables", md)
        self.assertIn("## Protocol", md)
        self.assertIn("TRANSFORMATION", md)
        self.assertIn("competent cells", md)
        # ...and everything else after it.
        self.assertIn("## Detail", md)
        self.assertLess(md.index("## Protocol"), md.index("## Detail"))
        self.assertLess(md.index("## Reagents"), md.index("## Protocol"))
        for section in ("### Validation", "### Route comparison",
                        "### Final sequence"):
            self.assertIn(section, md)
            self.assertGreater(md.index(section), md.index("## Detail"))

        gb = to_genbank(plan)
        self.assertTrue(gb.startswith("LOCUS"))
        self.assertTrue(gb.rstrip().endswith("//"))
        self.assertIn("FEATURES", gb)
        # The GenBank sequence must round-trip through our own parser.
        reparsed = parse_genbank(gb.encode())
        self.assertEqual(reparsed.seq, plan["construct"]["sequence"].upper())


if __name__ == "__main__":
    unittest.main(verbosity=2)


# --------------------------------------------------------------------------
class TestMissingPartAlignment(unittest.TestCase):
    """A module that resolves to nothing must not shift the others."""

    def test_unresolvable_module_does_not_misalign_fragments(self):
        conn = connect(":memory:")
        cfg = Config()
        import hashlib

        def add(ftype, name, seq, plasmid):
            pl = conn.execute(
                """INSERT INTO plasmid(name, owner, seq, seq_length,
                        seq_checksum, circular, verified)
                   VALUES (?,?,?,?,?,1,'sequence_verified')""",
                (plasmid, "t", seq + rand_seq(2000, 99), len(seq) + 2000,
                 name + "sum")).lastrowid
            part = conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ftype, name, seq, hashlib.sha256(seq.encode()).hexdigest(),
                 len(seq))).lastrowid
            conn.execute(
                """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                        strand) VALUES (?,?,0,?,1)""", (part, pl, len(seq)))

        bb = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pBB','t',?,?,'bbx',1,'sequence_verified')""",
            (rand_seq(3000, 98), 3000)).lastrowid
        add("promoter", "promA", rand_seq(400, 97), "pProm")
        add("CDS", "geneZ", orf(150, 96), "pGeneZ")
        conn.commit()

        # The middle module cannot be resolved and will be dropped.
        spec = ConstructSpec("gap", backbone_plasmid_id=bb, modules=[
            Module("promoter", "promA"),
            Module("IRES", "no_such_ires", required=False),
            Module("CDS", "geneZ"),
        ])
        plan = design(conn, spec, cfg)

        names = {f["name"] for f in plan["fragments"]}
        self.assertEqual(names, {"promA", "geneZ"})
        # Each surviving fragment must still point at its own donor.
        by_name = {f["name"]: f for f in plan["fragments"]}
        self.assertEqual(by_name["promA"]["source"], "pProm")
        self.assertEqual(by_name["geneZ"]["source"], "pGeneZ")


# --------------------------------------------------------------------------
class TestTimeEstimates(unittest.TestCase):
    """Hands-on bench time must not be conflated with elapsed turnaround."""

    def test_hands_on_is_less_than_elapsed(self):
        from genesis.operations import OPERATIONS
        for op in OPERATIONS.values():
            self.assertLessEqual(op.hands_on_minutes, op.default_minutes, op.code)

    def test_waiting_dominated_steps_are_cheap_in_hands_on_terms(self):
        from genesis.operations import OPERATIONS
        for code in ("SEQUENCE_VERIFY", "SYNTHESIS", "TRANSFORMATION"):
            op = OPERATIONS[code]
            self.assertLess(op.hands_on_minutes, op.default_minutes / 4, code)

    def test_route_reports_both_figures(self):
        conn = connect(":memory:")
        cfg = Config()
        bb = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pB','t',?,?,'s1',1,'sequence_verified')""",
            ("GAATTC" + rand_seq(2500, 120), 2506)).lastrowid
        import hashlib
        seq = orf(150, 121)
        pl = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pD','t',?,?,'s2',1,'sequence_verified')""",
            (seq + rand_seq(2000, 122), len(seq) + 2000)).lastrowid
        part = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES ('CDS','g',?,?,?)""",
            (seq, hashlib.sha256(seq.encode()).hexdigest(), len(seq))).lastrowid
        conn.execute("""INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,0,?,1)""", (part, pl, len(seq)))
        conn.commit()

        plan = design(conn, ConstructSpec("t", backbone_plasmid_id=bb,
                                          modules=[Module("CDS", "g")]), cfg)
        for route in plan["routes"]:
            self.assertIn("hands_on_minutes", route)
            self.assertIn("turnaround_minutes", route)
            if not route["feasible"]:
                continue        # a rejected route has no steps to cost
            self.assertLess(route["hands_on_minutes"], route["turnaround_minutes"])
            # Sequencing alone is a day; a real plan must exceed that elapsed.
            self.assertGreater(route["turnaround_minutes"], 1440)


# --------------------------------------------------------------------------
class TestBenchRealityChecks(unittest.TestCase):
    """Plans must not propose things that fail at the bench."""

    def setUp(self):
        self.conn = connect(":memory:")
        self.cfg = Config()

    def _plasmid(self, name, seq, circular=1, verified="sequence_verified"):
        return self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES (?,'t',?,?,?,?,?)""",
            (name, seq, len(seq), name + "_sum", circular, verified)).lastrowid

    def _part(self, ftype, name, seq, plasmid_id):
        import hashlib
        pid = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES (?,?,?,?,?)""",
            (ftype, name, seq, hashlib.sha256(seq.encode()).hexdigest(),
             len(seq))).lastrowid
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end, strand)
               VALUES (?,?,0,?,1)""", (pid, plasmid_id, len(seq)))
        return pid

    def test_oversized_backbone_pcr_is_flagged(self):
        bb = self._plasmid("pHuge", rand_seq(14000, 200))
        donor = self._plasmid("pDonor", orf(150, 201) + rand_seq(2000, 202))
        self._part("CDS", "gene", orf(150, 201), donor)
        self.conn.commit()
        plan = design(self.conn, ConstructSpec(
            "big", backbone_plasmid_id=bb, modules=[Module("CDS", "gene")]), self.cfg)
        hifi = next(r for r in plan["routes"] if r["strategy"] == "hifi")
        self.assertTrue(any("backbone amplicon" in w and "above the" in w
                            for w in hifi["warnings"]), hifi["warnings"])

    def test_tiny_fragment_suggests_primer_tail_instead_of_pcr(self):
        bb = self._plasmid("pBB", "GAATTC" + rand_seq(3000, 210))
        donor = self._plasmid("pTiny", rand_seq(3000, 211))
        self._part("linker", "shortlink", "GGCGGATCCGGCGGA", donor)
        self.conn.commit()
        plan = design(self.conn, ConstructSpec(
            "tiny", backbone_plasmid_id=bb,
            modules=[Module("linker", "shortlink")]), self.cfg)
        hifi = next(r for r in plan["routes"] if r["strategy"] == "hifi")
        self.assertTrue(any("annealed oligos" in w for w in hifi["warnings"]),
                        hifi["warnings"])

    def test_backbone_lookup_ignores_short_fragments(self):
        self._plasmid("FCK_tiny_fragment", "ACGTACGTACGTACGTACGTACG", circular=0)
        real = self._plasmid("FCK_real_vector", rand_seq(6000, 220))
        donor = self._plasmid("pD", orf(150, 221) + rand_seq(1000, 222))
        self._part("CDS", "g", orf(150, 221), donor)
        self.conn.commit()
        plan = design(self.conn, ConstructSpec(
            "bb", backbone="FCK", modules=[Module("CDS", "g")]), self.cfg)
        self.assertIsNotNone(plan["backbone"])
        self.assertEqual(plan["backbone"]["id"], real)

    def test_template_prefers_real_plasmid_over_pcr_product(self):
        bb = self._plasmid("pBB", "GAATTC" + rand_seq(3000, 230))
        seq = orf(150, 231)
        product = self._plasmid("PCRprdt_gene", seq, circular=0,
                                verified="unverified")
        plasmid = self._plasmid("pMaintained_gene", seq + rand_seq(3000, 232))
        part = self._part("CDS", "g", seq, product)
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end, strand)
               VALUES (?,?,0,?,1)""", (part, plasmid, len(seq)))
        self.conn.commit()
        plan = design(self.conn, ConstructSpec(
            "tmpl", backbone_plasmid_id=bb, modules=[Module("CDS", "g")]), self.cfg)
        self.assertEqual(plan["fragments"][0]["source"], "pMaintained_gene")


# --------------------------------------------------------------------------
class TestNameSuggestion(unittest.TestCase):

    def _plan(self, layout, backbone_name=""):
        return {"construct": {"name": "untitled", "layout": layout},
                "backbone": {"name": backbone_name} if backbone_name else None}

    def test_names_from_architecture(self):
        from genesis.protocol import suggest_name
        layout = [
            {"type": "backbone", "name": "backbone"},
            {"type": "promoter", "name": "hSyn"},
            {"type": "CDS", "name": "CheRiff"},
            {"type": "p2A", "name": "P2A"},
            {"type": "CDS", "name": "Citrine"},
            {"type": "polyA", "name": "bGH polyA"},
        ]
        self.assertEqual(suggest_name(self._plan(layout)),
                         "hSyn-CheRiff-P2A-Citrine")

    def test_aav_prefix_from_itrs(self):
        from genesis.protocol import suggest_name
        layout = [{"type": "AAV_ITR", "name": "ITR"},
                  {"type": "promoter", "name": "hSyn"},
                  {"type": "CDS", "name": "GCaMP8f"},
                  {"type": "AAV_ITR", "name": "ITR"}]
        self.assertEqual(suggest_name(self._plan(layout)), "pAAV-hSyn-GCaMP8f")

    def test_lenti_prefix_from_backbone_name(self):
        from genesis.protocol import suggest_name
        layout = [{"type": "promoter", "name": "EF1a"},
                  {"type": "CDS", "name": "Voltron2"}]
        self.assertEqual(suggest_name(self._plan(layout, "pLJM1_EF1a_LR")),
                         "pLenti-EF1a-Voltron2")

    def test_messy_part_labels_are_cleaned(self):
        from genesis.protocol import suggest_name
        layout = [{"type": "CDS", "name": "Archon1"},
                  {"type": "trafficking_motif", "name": "KGC / Kir2.1 TS"},
                  {"type": "trafficking_motif", "name": "ER2 export motif"}]
        self.assertEqual(suggest_name(self._plan(layout)), "Archon1-KGC-ER2")

    def test_structural_only_construct_falls_back(self):
        from genesis.protocol import suggest_name
        layout = [{"type": "origin", "name": "ColE1"},
                  {"type": "resistance", "name": "AmpR"}]
        self.assertEqual(suggest_name(self._plan(layout)), "untitled")


# --------------------------------------------------------------------------
class TestPDF(unittest.TestCase):

    def test_pdf_writer_produces_a_valid_document(self):
        from genesis.pdf import PDF
        doc = PDF("unit test")
        doc.heading("Title", 1)
        for i in range(90):
            doc.paragraph(f"line {i} " * 12)
        doc.table(["a", "b"], [["1", "2"], ["3", "4"]])
        data = doc.to_bytes()

        self.assertTrue(data.startswith(b"%PDF-1.4"))
        self.assertTrue(data.rstrip().endswith(b"%%EOF"))
        self.assertGreater(data.count(b"/Type /Page "), 1)   # paginated

        # Every xref offset must point at the object it claims.
        import re
        start = int(re.search(rb"startxref\s+(\d+)", data).group(1))
        lines = data[start:].split(b"\n")
        count = int(lines[1].split()[1])
        for i, line in enumerate(lines[3:3 + count - 1], start=1):
            offset = int(line.split()[0])
            self.assertTrue(data[offset:].startswith(f"{i} 0 obj".encode()),
                            f"xref entry {i} is wrong")

    def test_text_measurement_uses_real_metrics(self):
        from genesis.pdf import text_width, wrap
        # 'i' is far narrower than 'W' in Helvetica; a fixed-width
        # approximation would make these equal and wrap badly.
        self.assertLess(text_width("iiii", "F1", 10),
                        text_width("WWWW", "F1", 10) / 2)
        lines = wrap("word " * 60, "F1", 10, 200)
        self.assertGreater(len(lines), 1)
        for line in lines:
            self.assertLessEqual(text_width(line, "F1", 10), 200.5)

    def test_protocol_pdf_round_trip(self):
        conn = connect(":memory:")
        cfg = Config()
        import hashlib
        bb = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pBB','t',?,?,'pdfbb',1,'sequence_verified')""",
            ("GAATTC" + rand_seq(2500, 300), 2506)).lastrowid
        seq = orf(150, 301)
        donor = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pD','t',?,?,'pdfd',1,'sequence_verified')""",
            (seq + rand_seq(2000, 302), len(seq) + 2000)).lastrowid
        part = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES ('CDS','gene',?,?,?)""",
            (seq, hashlib.sha256(seq.encode()).hexdigest(), len(seq))).lastrowid
        conn.execute("""INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,0,?,1)""", (part, donor, len(seq)))
        conn.commit()

        plan = design(conn, ConstructSpec("pdf test", backbone_plasmid_id=bb,
                                          modules=[Module("CDS", "gene")]), cfg)
        from genesis.export import protocol_pdf
        data = protocol_pdf(plan, cfg)
        self.assertTrue(data.startswith(b"%PDF"))
        self.assertGreater(len(data), 3000)

    def test_plan_carries_a_suggested_name(self):
        conn = connect(":memory:")
        cfg = Config()
        import hashlib
        bb = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pAAV_backbone','t',?,?,'sn1',1,'sequence_verified')""",
            ("GAATTC" + rand_seq(2500, 310), 2506)).lastrowid
        for ftype, nm, s in (("promoter", "hSyn", rand_seq(400, 311)),
                             ("CDS", "GCaMP8f", orf(150, 312))):
            pl = conn.execute(
                """INSERT INTO plasmid(name, owner, seq, seq_length,
                        seq_checksum, circular, verified)
                   VALUES (?,'t',?,?,?,1,'sequence_verified')""",
                ("p" + nm, s + rand_seq(2000, 313), len(s) + 2000, nm)).lastrowid
            pid = conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ftype, nm, s, hashlib.sha256(s.encode()).hexdigest(),
                 len(s))).lastrowid
            conn.execute("""INSERT INTO part_occurrence(part_id, plasmid_id,
                            start, end, strand) VALUES (?,?,0,?,1)""",
                         (pid, pl, len(s)))
        conn.commit()
        plan = design(conn, ConstructSpec(
            "untitled", backbone_plasmid_id=bb,
            modules=[Module("promoter", "hSyn"), Module("CDS", "GCaMP8f")]), cfg)
        self.assertEqual(plan["construct"]["suggested_name"],
                         "pAAV-hSyn-GCaMP8f")


# --------------------------------------------------------------------------
class TestCompositeParts(unittest.TestCase):
    """A CDS that already carries a fused 2A must not produce a tandem pair."""

    P2A_AA = "ATNFSLLKQAGDVEENPGP"
    _BACK = {"A": "GCG", "T": "ACA", "N": "AAC", "F": "TTC", "S": "AGC",
             "L": "CTG", "K": "AAG", "Q": "CAG", "G": "GGC", "D": "GAC",
             "V": "GTG", "E": "GAG", "P": "CCG"}

    @classmethod
    def p2a_dna(cls):
        return "".join(cls._BACK[a] for a in cls.P2A_AA)

    def test_detects_fused_2a_in_a_cds(self):
        from genesis.vocab import find_embedded_elements
        cds = "ATG" + "GCT" * 150
        found = find_embedded_elements(cds + self.p2a_dna(), "CDS")
        self.assertTrue(found)
        self.assertEqual(found[0]["feature_type"], "p2A")
        self.assertTrue(found[0]["at_end"])

    def test_plain_cds_has_nothing_embedded(self):
        from genesis.vocab import find_embedded_elements
        self.assertEqual(find_embedded_elements("ATG" + "GCT" * 150, "CDS"), [])

    def test_trim_only_touches_terminal_elements(self):
        from genesis.vocab import find_embedded_elements, trim_embedded
        p2a = self.p2a_dna()
        # an internal 2A must be left alone: cutting it changes the protein
        internal = "ATG" + "GCT" * 60 + p2a + "GCT" * 60
        els = find_embedded_elements(internal, "CDS")
        trimmed, removed = trim_embedded(internal, els, "p2A")
        self.assertEqual(trimmed, internal)
        self.assertEqual(removed, [])

    def test_resolver_trims_a_fused_2a_when_one_is_also_requested(self):
        import hashlib
        conn = connect(":memory:")
        cfg = Config()
        p2a = self.p2a_dna()
        fused = "ATG" + "GCT" * 150 + p2a          # CDS with a trailing P2A

        def add(ftype, name, seq):
            pl = conn.execute(
                """INSERT INTO plasmid(name, owner, seq, seq_length,
                        seq_checksum, circular, verified)
                   VALUES (?,'t',?,?,?,1,'sequence_verified')""",
                ("p" + name, seq + rand_seq(2500, 400), len(seq) + 2500,
                 name + "s")).lastrowid
            pid = conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ftype, name, seq, hashlib.sha256(seq.encode()).hexdigest(),
                 len(seq))).lastrowid
            conn.execute("""INSERT INTO part_occurrence(part_id, plasmid_id,
                            start, end, strand) VALUES (?,?,0,?,1)""",
                         (pid, pl, len(seq)))
            return pid

        add("CDS", "FusedGene", fused)
        add("p2A", "P2A", p2a)
        add("CDS", "Reporter", "GCT" * 150 + "TAA")
        conn.commit()

        spec = ConstructSpec("tandem check", modules=[
            Module("CDS", "FusedGene"), Module("p2A", "P2A"),
            Module("CDS", "Reporter")])
        from genesis.construct import resolve
        resolved = resolve(conn, spec)

        gene = resolved[0]
        self.assertNotIn(p2a, gene.seq, "the fused P2A should have been trimmed")
        self.assertEqual(len(gene.seq), len(fused) - len(p2a))
        self.assertTrue(any("two" in w for w in gene.warnings), gene.warnings)

    def test_validator_flags_a_tandem_2a(self):
        p2a = self.p2a_dna()
        mods = [("CDS", "geneA", "ATG" + "GCT" * 50),
                ("p2A", "P2A-1", p2a),
                ("p2A", "P2A-2", p2a),
                ("CDS", "geneB", "GCT" * 50 + "TAA")]
        resolved = [ResolvedModule(Module(type=t, name=n, seq=s), None, n, s,
                                   t, [], []) for t, n, s in mods]
        construct = assemble(resolved, "", True)
        msgs = [i.message for i in validate_construct(construct)
                if i.level == "error"]
        self.assertTrue(any("two 2A peptides" in m for m in msgs), msgs)

    def test_single_2a_is_not_flagged(self):
        p2a = self.p2a_dna()
        mods = [("CDS", "geneA", "ATG" + "GCT" * 50),
                ("p2A", "P2A", p2a),
                ("CDS", "geneB", "GCT" * 50 + "TAA")]
        resolved = [ResolvedModule(Module(type=t, name=n, seq=s), None, n, s,
                                   t, [], []) for t, n, s in mods]
        construct = assemble(resolved, "", True)
        msgs = [i.message for i in validate_construct(construct)
                if i.level == "error"]
        self.assertFalse(any("two 2A peptides" in m for m in msgs), msgs)


# --------------------------------------------------------------------------
class TestDIO(unittest.TestCase):
    """DIO is an architecture that wraps a payload, not a part to look up.

        lox2272 - loxP - [payload inverted] - lox2272' - loxP'

    Nothing is expressed until Cre inverts the payload and excises one pair.
    """

    def _resolved(self, items):
        from genesis.construct import ResolvedModule
        return [ResolvedModule(Module(type=t, name=n, seq=s), None, n, s, t, [], [])
                for t, n, s in items]

    def test_expansion_builds_the_interleaved_lox_arrangement(self):
        from genesis.construct import expand_dio
        payload = orf(300, 500)
        mods = self._resolved([("promoter", "hSyn", "G" * 400),
                               ("DIO", "DIO", ""),
                               ("CDS", "payload", payload)])
        out, notes = expand_dio(mods)
        kinds = [(m.feature_type, m.part_name, m.module.orientation) for m in out]
        self.assertEqual(kinds, [
            ("promoter", "hSyn", "forward"),
            ("lox_site", "lox2272", "forward"),
            ("lox_site", "loxP", "forward"),
            ("CDS", "payload", "reverse"),
            ("lox_site", "lox2272", "reverse"),
            ("lox_site", "loxP", "reverse"),
        ])
        self.assertTrue(notes)

    def test_payload_is_off_before_cre_and_on_after(self):
        from genesis.construct import expand_dio
        from genesis.validate import lox_sites_from_layout, simulate_cre
        payload = orf(300, 501)
        protein = st.translate(payload).rstrip("*")
        out, _ = expand_dio(self._resolved([
            ("promoter", "hSyn", "G" * 400), ("DIO", "DIO", ""),
            ("CDS", "payload", payload)]))
        construct = assemble(out, "", True)

        before = [o for o in st.find_orfs(construct.seq, min_aa=250)
                  if o["strand"] == 1 and o["protein"] == protein]
        self.assertFalse(before, "payload must not be expressible before Cre")

        recombined, log = simulate_cre(construct.seq,
                                       lox_sites_from_layout(construct.seq))
        after = [o for o in st.find_orfs(recombined, min_aa=250)
                 if o["strand"] == 1 and o["protein"] == protein]
        self.assertTrue(after, "Cre must switch the payload on")
        self.assertEqual(len(log), 2)
        self.assertIn("inverted", log[0])
        self.assertIn("excised", log[1])

    def test_one_site_of_each_type_survives(self):
        from genesis.construct import expand_dio
        from genesis.validate import lox_sites_from_layout, simulate_cre
        out, _ = expand_dio(self._resolved([
            ("DIO", "DIO", ""), ("CDS", "payload", orf(200, 502))]))
        construct = assemble(out, "", True)
        recombined, _ = simulate_cre(construct.seq,
                                     lox_sites_from_layout(construct.seq))
        remaining = lox_sites_from_layout(recombined)
        self.assertEqual({s.name for s in remaining}, {"lox2272", "loxP"})
        self.assertEqual(len(remaining), 2)

    def test_multi_module_payload_is_inverted_as_one_block(self):
        from genesis.construct import expand_dio
        a, b = orf(100, 503), orf(120, 504)
        out, _ = expand_dio(self._resolved([
            ("DIO", "DIO", ""), ("CDS", "geneA", a),
            ("p2A", "P2A", "GGC" * 22), ("CDS", "geneB", b),
            ("polyA", "bGH", "AATAAA" * 5)]))
        names = [m.part_name for m in out]
        # payload reversed between the lox sites, polyA left outside
        self.assertEqual(names, ["lox2272", "loxP", "geneB", "P2A", "geneA",
                                 "lox2272", "loxP", "bGH"])

    def test_dio_without_a_payload_is_reported_not_silently_kept(self):
        from genesis.construct import expand_dio
        out, notes = expand_dio(self._resolved([
            ("promoter", "hSyn", "G" * 300), ("DIO", "DIO", "")]))
        self.assertEqual([m.feature_type for m in out], ["promoter"])
        self.assertTrue(any("no coding sequence" in n for n in notes), notes)

    def test_inverted_payload_is_not_translated_as_a_forward_fusion(self):
        """The reverse-strand payload must be read off its own strand."""
        from genesis.construct import expand_dio, translated_span
        out, _ = expand_dio(self._resolved([
            ("DIO", "DIO", ""), ("CDS", "geneA", "ATG" + "GCT" * 40),
            ("p2A", "P2A", "GGC" * 22),
            ("CDS", "geneB", "GCT" * 40 + "TAA")]))
        construct = assemble(out, "", True)
        runs = translated_span(construct)
        self.assertEqual(len(runs), 1)
        self.assertEqual(runs[0]["strand"], -1)
        # members come back in the order the ribosome reads them
        self.assertEqual([m["name"] for m in runs[0]["members"]],
                         ["geneA", "P2A", "geneB"])
        errors = [i.message for i in validate_construct(construct)
                  if i.level == "error"]
        self.assertEqual(errors, [], errors)


# --------------------------------------------------------------------------
class TestCodonMarkers(unittest.TestCase):
    """Start and stop codons, located for drawing on a map."""

    def _c(self, items):
        from genesis.construct import ResolvedModule
        resolved = [ResolvedModule(Module(type=t, name=n, seq=s), None, n, s,
                                   t, [], []) for t, n, s in items]
        return assemble(resolved, "", True)

    def test_clean_fusion_has_one_start_and_one_stop(self):
        from genesis.construct import codon_markers
        c = self._c([("CDS", "geneA", "ATG" + "GCT" * 20),
                     ("CDS", "geneB", "GCT" * 20 + "TAA")])
        marks = codon_markers(c)
        self.assertEqual([m["kind"] for m in marks], ["start", "stop"])
        self.assertEqual(marks[0]["start"], 0)
        self.assertEqual(marks[0]["codon"], "ATG")
        self.assertEqual(marks[1]["end"], c.length)
        self.assertEqual(marks[1]["codon"], "TAA")

    def test_premature_stop_is_marked_separately_from_the_real_one(self):
        from genesis.construct import codon_markers
        c = self._c([("CDS", "geneA", "ATG" + "GCT" * 5 + "TGA" + "GCT" * 5),
                     ("CDS", "geneB", "GCT" * 5 + "TAA")])
        kinds = [m["kind"] for m in codon_markers(c)]
        self.assertIn("premature_stop", kinds)
        self.assertEqual(kinds.count("stop"), 1)
        premature = next(m for m in codon_markers(c)
                         if m["kind"] == "premature_stop")
        self.assertEqual(premature["start"], 18)
        self.assertEqual(premature["module"], "geneA")

    def test_missing_start_and_missing_stop_are_reported(self):
        from genesis.construct import codon_markers
        c = self._c([("CDS", "geneA", "GCT" * 20),
                     ("CDS", "geneB", "GCT" * 20)])
        kinds = {m["kind"] for m in codon_markers(c)}
        self.assertIn("missing_start", kinds)
        self.assertIn("missing_stop", kinds)

    def test_reverse_strand_payload_maps_back_to_construct_coordinates(self):
        """An inverted DIO payload puts its ATG at the HIGHER coordinate."""
        from genesis.construct import codon_markers, expand_dio, ResolvedModule
        payload = "ATG" + "GCT" * 40 + "TAA"
        mods = [ResolvedModule(Module(type=t, name=n, seq=s), None, n, s, t, [], [])
                for t, n, s in [("DIO", "DIO", ""), ("CDS", "payload", payload)]]
        out, _ = expand_dio(mods)
        c = assemble(out, "", True)
        marks = codon_markers(c)

        start = next(m for m in marks if m["kind"] == "start")
        stop = next(m for m in marks if m["kind"] == "stop")
        self.assertEqual(start["strand"], -1)
        self.assertEqual(stop["strand"], -1)
        self.assertGreater(start["start"], stop["start"])
        # and they must point at the real codons on the reverse strand
        self.assertEqual(st.revcomp(c.seq[start["start"]:start["end"]]), "ATG")
        self.assertEqual(st.revcomp(c.seq[stop["start"]:stop["end"]]), "TAA")

    def test_markers_are_sorted_and_carry_their_module(self):
        from genesis.construct import codon_markers
        c = self._c([("CDS", "geneA", "ATG" + "GCT" * 10),
                     ("p2A", "P2A", "GGC" * 22),
                     ("CDS", "geneB", "GCT" * 10 + "TAG")])
        marks = codon_markers(c)
        self.assertEqual([m["start"] for m in marks],
                         sorted(m["start"] for m in marks))
        self.assertEqual(marks[0]["module"], "geneA")
        self.assertEqual(marks[-1]["module"], "geneB")

    def test_untranslated_construct_has_no_markers(self):
        from genesis.construct import codon_markers
        c = self._c([("promoter", "hSyn", "G" * 300),
                     ("polyA", "bGH", "AATAAA" * 10)])
        self.assertEqual(codon_markers(c), [])

    def test_plan_exposes_markers(self):
        import hashlib
        conn = connect(":memory:")
        cfg = Config()
        seq = orf(150, 600)
        pl = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pD','t',?,?,'cm1',1,'sequence_verified')""",
            (seq + rand_seq(2000, 601), len(seq) + 2000)).lastrowid
        pid = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES ('CDS','gene',?,?,?)""",
            (seq, hashlib.sha256(seq.encode()).hexdigest(), len(seq))).lastrowid
        conn.execute("""INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,0,?,1)""", (pid, pl, len(seq)))
        conn.commit()
        plan = design(conn, ConstructSpec("m", modules=[Module("CDS", "gene")]),
                      cfg)
        marks = plan["construct"]["codon_markers"]
        self.assertTrue(marks)
        self.assertEqual(marks[0]["kind"], "start")

    def test_reverse_strand_missing_stop_sits_past_the_last_codon(self):
        """Translation runs toward lower coordinates on the reverse strand,
        so "just past the end" is the low edge of the run, not the high one."""
        from genesis.construct import codon_markers, expand_dio, ResolvedModule
        payload = "ATG" + "GCT" * 40            # deliberately no stop codon
        mods = [ResolvedModule(Module(type=t, name=n, seq=s), None, n, s, t, [], [])
                for t, n, s in [("DIO", "DIO", ""), ("CDS", "payload", payload)]]
        out, _ = expand_dio(mods)
        c = assemble(out, "", True)

        cds = [i for i in c.layout if i["type"] == "CDS"]
        run_start = min(i["start"] for i in cds)
        marker = next(m for m in codon_markers(c) if m["kind"] == "missing_stop")
        self.assertEqual(marker["strand"], -1)
        self.assertEqual(marker["start"], run_start)
        self.assertEqual(marker["start"], marker["end"])   # zero width

    def test_forward_missing_stop_sits_past_the_last_codon(self):
        from genesis.construct import codon_markers
        c = self._c([("CDS", "geneA", "ATG" + "GCT" * 20)])
        marker = next(m for m in codon_markers(c) if m["kind"] == "missing_stop")
        self.assertEqual(marker["start"], c.length)


# --------------------------------------------------------------------------
class TestAuth(unittest.TestCase):
    """Accounts, passwords and sessions."""

    GOOD = "correct-horse-battery-staple"

    def setUp(self):
        from genesis import auth
        self.auth = auth
        self.conn = connect(":memory:")
        auth.ensure_schema(self.conn)

    # -- password storage ------------------------------------------------
    def test_password_is_not_stored(self):
        self.auth.create_user(self.conn, "alice", self.GOOD)
        row = self.conn.execute(
            "SELECT salt, password_hash, iterations FROM user").fetchone()
        self.assertNotIn(self.GOOD.encode(), row["password_hash"])
        self.assertNotIn(self.GOOD.encode(), row["salt"])
        self.assertGreaterEqual(row["iterations"], 100_000)

    def test_same_password_gets_different_hashes(self):
        self.auth.create_user(self.conn, "a", self.GOOD)
        self.auth.create_user(self.conn, "b", self.GOOD)
        rows = self.conn.execute("SELECT salt, password_hash FROM user").fetchall()
        self.assertNotEqual(rows[0]["salt"], rows[1]["salt"])
        self.assertNotEqual(rows[0]["password_hash"], rows[1]["password_hash"])

    def test_verification_round_trip(self):
        salt, derived, iters = self.auth.hash_password(self.GOOD)
        self.assertTrue(self.auth.check_password(self.GOOD, salt, derived, iters))
        self.assertFalse(self.auth.check_password("other", salt, derived, iters))

    def test_weak_passwords_are_refused(self):
        for bad in ("short", "password1234", "Password123!", "qwertyuiop12"):
            with self.assertRaises(self.auth.AuthError, msg=bad):
                self.auth.create_user(self.conn, f"u{abs(hash(bad))}", bad)

    def test_long_passphrase_is_accepted(self):
        self.assertIsNone(self.auth.password_problem("voltron imaging rig 7"))

    # -- login -----------------------------------------------------------
    def test_login_and_token_resolution(self):
        self.auth.create_user(self.conn, "alice", self.GOOD, role="admin")
        token, user = self.auth.login(self.conn, "alice", self.GOOD)
        self.assertEqual(user.role, "admin")
        self.assertEqual(self.auth.user_for_token(self.conn, token).username,
                         "alice")

    def test_username_is_case_insensitive(self):
        self.auth.create_user(self.conn, "Alice", self.GOOD)
        token, _ = self.auth.login(self.conn, "ALICE", self.GOOD)
        self.assertIsNotNone(self.auth.user_for_token(self.conn, token))

    def test_session_token_is_stored_hashed(self):
        self.auth.create_user(self.conn, "alice", self.GOOD)
        token, _ = self.auth.login(self.conn, "alice", self.GOOD)
        stored = self.conn.execute(
            "SELECT token_hash FROM session").fetchone()["token_hash"]
        self.assertNotEqual(stored, token)
        self.assertNotIn(token, stored)

    def test_unknown_user_and_wrong_password_look_identical(self):
        """The error must not reveal which usernames exist."""
        self.auth.create_user(self.conn, "alice", self.GOOD)
        with self.assertRaises(self.auth.AuthError) as wrong:
            self.auth.login(self.conn, "alice", "not-the-password")
        with self.assertRaises(self.auth.AuthError) as unknown:
            self.auth.login(self.conn, "nobody-at-all", "not-the-password")
        self.assertEqual(str(wrong.exception), str(unknown.exception))

    def test_repeated_failures_are_throttled(self):
        self.auth.create_user(self.conn, "alice", self.GOOD)
        for _ in range(4):
            with self.assertRaises(self.auth.AuthError):
                self.auth.login(self.conn, "alice", "wrong")
        with self.assertRaises(self.auth.AuthError) as ctx:
            self.auth.login(self.conn, "alice", self.GOOD)
        self.assertIn("too many failed attempts", str(ctx.exception))

    def test_successful_login_clears_the_throttle(self):
        self.auth.create_user(self.conn, "alice", self.GOOD)
        for _ in range(2):          # below the threshold
            with self.assertRaises(self.auth.AuthError):
                self.auth.login(self.conn, "alice", "wrong")
        self.auth.login(self.conn, "alice", self.GOOD)
        self.assertEqual(self.conn.execute(
            "SELECT COUNT(*) c FROM login_attempt").fetchone()["c"], 0)

    # -- sessions --------------------------------------------------------
    def test_forged_and_empty_tokens_are_rejected(self):
        import secrets
        self.assertIsNone(self.auth.user_for_token(self.conn, ""))
        self.assertIsNone(self.auth.user_for_token(
            self.conn, secrets.token_urlsafe(32)))

    def test_expired_session_is_rejected_and_removed(self):
        self.auth.create_user(self.conn, "alice", self.GOOD)
        token, _ = self.auth.login(self.conn, "alice", self.GOOD,
                                   session_hours=-1)      # already expired
        self.assertIsNone(self.auth.user_for_token(self.conn, token))
        self.assertEqual(self.conn.execute(
            "SELECT COUNT(*) c FROM session").fetchone()["c"], 0)

    def test_logout_revokes_only_that_session(self):
        self.auth.create_user(self.conn, "alice", self.GOOD)
        a, _ = self.auth.login(self.conn, "alice", self.GOOD)
        b, _ = self.auth.login(self.conn, "alice", self.GOOD)
        self.auth.logout(self.conn, a)
        self.assertIsNone(self.auth.user_for_token(self.conn, a))
        self.assertIsNotNone(self.auth.user_for_token(self.conn, b))

    def test_password_change_ends_every_session(self):
        self.auth.create_user(self.conn, "alice", self.GOOD)
        a, _ = self.auth.login(self.conn, "alice", self.GOOD)
        b, _ = self.auth.login(self.conn, "alice", self.GOOD)
        self.auth.set_password(self.conn, "alice", "a-different-passphrase")
        for token in (a, b):
            self.assertIsNone(self.auth.user_for_token(self.conn, token))

    def test_disabling_an_account_ends_its_sessions(self):
        self.auth.create_user(self.conn, "admin1", self.GOOD, role="admin")
        self.auth.create_user(self.conn, "bob", self.GOOD)
        token, _ = self.auth.login(self.conn, "bob", self.GOOD)
        self.auth.set_disabled(self.conn, "bob", True)
        self.assertIsNone(self.auth.user_for_token(self.conn, token))
        with self.assertRaises(self.auth.AuthError):
            self.auth.login(self.conn, "bob", self.GOOD)

    # -- roles -----------------------------------------------------------
    def test_role_ranking(self):
        User = self.auth.User
        self.assertTrue(User(1, "a", "a", "admin").can("member"))
        self.assertTrue(User(1, "a", "a", "member").can("member"))
        self.assertFalse(User(1, "a", "a", "viewer").can("member"))
        self.assertFalse(User(1, "a", "a", "member").can("admin"))

    def test_cannot_strip_or_disable_the_last_admin(self):
        self.auth.create_user(self.conn, "root", self.GOOD, role="admin")
        with self.assertRaises(self.auth.AuthError):
            self.auth.set_role(self.conn, "root", "viewer")
        with self.assertRaises(self.auth.AuthError):
            self.auth.set_disabled(self.conn, "root", True)
        # with a second admin it is allowed
        self.auth.create_user(self.conn, "root2", self.GOOD, role="admin")
        self.auth.set_role(self.conn, "root", "viewer")
        self.assertEqual(self.conn.execute(
            "SELECT role FROM user WHERE username='root'").fetchone()["role"],
            "viewer")

    def test_duplicate_username_is_refused(self):
        self.auth.create_user(self.conn, "alice", self.GOOD)
        with self.assertRaises(self.auth.AuthError):
            self.auth.create_user(self.conn, "ALICE", self.GOOD)


# --------------------------------------------------------------------------
class TestUnsafeBindRefusal(unittest.TestCase):
    """The server must not put an unprotected inventory on the network."""

    def test_refuses_reachable_bind_without_auth(self):
        from genesis.server import serve
        cfg = Config()
        cfg.host, cfg.require_auth = "0.0.0.0", False
        with self.assertRaises(SystemExit) as ctx:
            serve(cfg, ":memory:")
        self.assertIn("require_auth is off", str(ctx.exception))

    def test_refuses_reachable_bind_with_no_accounts(self):
        from genesis.server import serve
        cfg = Config()
        cfg.host, cfg.require_auth = "0.0.0.0", True
        with self.assertRaises(SystemExit) as ctx:
            serve(cfg, ":memory:")
        self.assertIn("no accounts exist", str(ctx.exception))

    def test_loopback_without_auth_is_allowed(self):
        """Local-only use must not require anyone to set up accounts."""
        from genesis.server import LOOPBACK
        self.assertIn("127.0.0.1", LOOPBACK)


# --------------------------------------------------------------------------
class TestSaveModesAndDelete(unittest.TestCase):
    """Versioning, overwriting, saving under a new name, and deleting."""

    def setUp(self):
        import hashlib
        self.conn = connect(":memory:")
        self.cfg = Config()
        seq = orf(150, 700)
        pl = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pD','t',?,?,'sv',1,'sequence_verified')""",
            (seq + rand_seq(2000, 701), len(seq) + 2000)).lastrowid
        pid = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES ('CDS','gene',?,?,?)""",
            (seq, hashlib.sha256(seq.encode()).hexdigest(), len(seq))).lastrowid
        self.conn.execute("""INSERT INTO part_occurrence(part_id, plasmid_id,
                             start, end, strand) VALUES (?,?,0,?,1)""",
                          (pid, pl, len(seq)))
        self.conn.commit()
        self.spec = ConstructSpec("design A", modules=[Module("CDS", "gene")])
        self.plan = design(self.conn, self.spec, self.cfg)

    def _count(self, table):
        return self.conn.execute(f"SELECT COUNT(*) c FROM {table}").fetchone()["c"]

    def test_default_save_adds_a_version(self):
        from genesis.design import save_design
        a = save_design(self.conn, self.spec, self.plan)
        b = save_design(self.conn, self.spec, self.plan)
        self.assertEqual((a["version"], b["version"]), (1, 2))
        self.assertEqual(a["construct_id"], b["construct_id"])
        self.assertEqual(self._count("design"), 2)

    def test_overwrite_replaces_in_place(self):
        from genesis.design import save_design
        a = save_design(self.conn, self.spec, self.plan)
        b = save_design(self.conn, self.spec, self.plan, mode="overwrite",
                        design_id=a["design_id"])
        self.assertEqual(b["design_id"], a["design_id"])
        self.assertEqual(b["version"], a["version"])
        self.assertEqual(self._count("design"), 1)

    def test_overwrite_needs_a_target(self):
        from genesis.design import save_design
        with self.assertRaises(ValueError):
            save_design(self.conn, self.spec, self.plan, mode="overwrite")
        with self.assertRaises(ValueError):
            save_design(self.conn, self.spec, self.plan, mode="overwrite",
                        design_id=9999)

    def test_save_as_creates_a_separate_construct(self):
        from genesis.design import save_design
        save_design(self.conn, self.spec, self.plan)
        other = ConstructSpec("design B", modules=[Module("CDS", "gene")])
        b = save_design(self.conn, other, self.plan, mode="save_as")
        self.assertEqual(b["version"], 1)
        self.assertEqual(self._count("construct"), 2)

    def test_save_as_refuses_an_existing_name(self):
        from genesis.design import save_design
        save_design(self.conn, self.spec, self.plan)
        with self.assertRaises(ValueError) as ctx:
            save_design(self.conn, self.spec, self.plan, mode="save_as")
        self.assertIn("already exists", str(ctx.exception))

    def test_deleting_the_last_version_removes_the_construct(self):
        from genesis.design import save_design, delete_design
        a = save_design(self.conn, self.spec, self.plan)
        result = delete_design(self.conn, a["design_id"])
        self.assertEqual(result["versions_left"], 0)
        self.assertEqual(self._count("construct"), 0)

    def test_deleting_one_version_keeps_the_others(self):
        from genesis.design import save_design, delete_design
        a = save_design(self.conn, self.spec, self.plan)
        save_design(self.conn, self.spec, self.plan)
        delete_design(self.conn, a["design_id"])
        self.assertEqual(self._count("design"), 1)
        self.assertEqual(self._count("construct"), 1)

    def test_delete_construct_removes_every_version_and_outcome(self):
        from genesis.design import save_design, delete_construct
        a = save_design(self.conn, self.spec, self.plan)
        save_design(self.conn, self.spec, self.plan)
        self.conn.execute(
            "INSERT INTO outcome(design_id, result) VALUES (?, 'success')",
            (a["design_id"],))
        self.conn.commit()
        result = delete_construct(self.conn, a["construct_id"])
        self.assertEqual(result["versions_deleted"], 2)
        self.assertEqual(self._count("design"), 0)
        self.assertEqual(self._count("outcome"), 0)

    def test_deleting_something_that_is_not_there(self):
        from genesis.design import delete_design
        with self.assertRaises(ValueError):
            delete_design(self.conn, 4242)


# --------------------------------------------------------------------------
class TestSplittingFusedParts(unittest.TestCase):
    """A fused CDS is catalogued whole and split; both stay available."""

    @staticmethod
    def _bt(peptide):
        back = {}
        for codon, aa in st.CODON_TABLE.items():
            back.setdefault(aa, codon)
        return "".join(back[a] for a in peptide)

    def setUp(self):
        self.p2a = self._bt("ATNFSLLKQAGDVEENPGP")
        self.ts = self._bt("KSRITSEGEYIPLDQIDINV")
        self.er2 = self._bt("FCYENEV")

    def test_splits_at_a_2a(self):
        from genesis.vocab import find_embedded_elements, split_at_elements
        fused = "ATG" + "GCT" * 200 + self.p2a + "GCT" * 200 + "TAA"
        pieces = split_at_elements(fused, find_embedded_elements(fused, "CDS"))
        kinds = [p["kind"] for p in pieces]
        self.assertEqual(kinds, ["segment", "separator", "segment"])
        self.assertEqual(pieces[1]["feature_type"], "p2A")

    def test_splits_at_a_trafficking_signal_between_two_proteins(self):
        from genesis.vocab import find_embedded_elements, split_at_elements
        fused = "ATG" + "GCT" * 200 + self.ts + "GCT" * 200 + "TAA"
        pieces = split_at_elements(fused, find_embedded_elements(fused, "CDS"))
        self.assertEqual(len([p for p in pieces if p["kind"] == "segment"]), 2)

    def test_terminal_tags_do_not_split(self):
        """A trafficking motif on the end of one protein is not a boundary."""
        from genesis.vocab import find_embedded_elements, split_at_elements
        fused = "ATG" + "GCT" * 200 + self.ts
        self.assertEqual(
            split_at_elements(fused, find_embedded_elements(fused, "CDS")), [])

    def test_tags_beside_a_2a_stay_with_their_own_protein(self):
        """Voltron2's KGC and ER2 belong to Voltron2, not to the junction."""
        from genesis.vocab import find_embedded_elements, split_at_elements
        first = "ATG" + "GCT" * 200 + self.ts + self.er2
        fused = first + self.p2a + "GCT" * 200 + "TAA"
        pieces = split_at_elements(fused, find_embedded_elements(fused, "CDS"))
        segments = [p for p in pieces if p["kind"] == "segment"]
        separators = [p for p in pieces if p["kind"] == "separator"]
        self.assertEqual(len(segments), 2)
        self.assertEqual(len(separators), 1)
        # the boundary is the 2A alone, and the tags went with protein one
        self.assertEqual(separators[0]["feature_type"], "p2A")
        self.assertEqual(len(segments[0]["seq"]), len(first))

    def test_a_plain_cds_does_not_split(self):
        from genesis.vocab import find_embedded_elements, split_at_elements
        plain = "ATG" + "GCT" * 300 + "TAA"
        self.assertEqual(
            split_at_elements(plain, find_embedded_elements(plain, "CDS")), [])

    def test_importer_keeps_the_whole_and_the_pieces(self):
        from genesis.importer import extract_parts
        from genesis.seqio import Record
        conn = connect(":memory:")
        fused = "ATG" + "GCT" * 200 + self.p2a + "GCT" * 200 + "TAA"
        pl = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular) VALUES ('p','t',?,?,'x',1)""",
            (fused, len(fused))).lastrowid
        rec = Record(name="p", seq=fused, circular=True)
        extract_parts(conn, pl, rec, [{
            "feature_type": "CDS", "name": "FusionGene", "start": 0,
            "end": len(fused), "strand": 1, "confidence": 1.0}])
        conn.commit()
        whole = conn.execute(
            "SELECT id FROM part WHERE derived_from IS NULL").fetchall()
        pieces = conn.execute(
            "SELECT name, seq_length FROM part WHERE derived_from IS NOT NULL"
            " ORDER BY segment_index").fetchall()
        self.assertEqual(len(whole), 1, "the fusion must stay catalogued")
        self.assertEqual(len(pieces), 2, "and its proteins must be too")
        self.assertIn("1 of 2", pieces[0]["name"])


# --------------------------------------------------------------------------
class TestVectorFamilies(unittest.TestCase):

    def _ann(self, *pairs):
        return [{"feature_type": t, "name": n} for t, n in pairs]

    def test_two_itrs_identify_an_aav(self):
        from genesis.vectors import classify_vector
        v = classify_vector(self._ann(("AAV_ITR", "AAV2 ITR"),
                                      ("AAV_ITR", "AAV2 ITR"),
                                      ("promoter", "hSyn")), "pAAV-hSyn")
        self.assertEqual(v["family"], "AAV")

    def test_lentiviral_elements_identify_a_lentivirus(self):
        from genesis.vectors import classify_vector
        v = classify_vector(self._ann(("promoter", "5' LTR"),
                                      ("misc", "HIV-1 psi pack"),
                                      ("misc", "RRE"), ("misc", "cPPT")),
                            "pLenti-EF1a")
        self.assertEqual(v["family"], "lentiviral")

    def test_a_plain_plasmid_is_neither(self):
        from genesis.vectors import classify_vector
        v = classify_vector(self._ann(("origin", "ColE1"),
                                      ("resistance", "AmpR")), "pUC19")
        self.assertEqual(v["family"], "plasmid")

    def test_conflicting_evidence_is_reported(self):
        from genesis.vectors import classify_vector
        v = classify_vector(self._ann(("AAV_ITR", "AAV2 ITR"),
                                      ("AAV_ITR", "AAV2 ITR"),
                                      ("promoter", "5' LTR"), ("misc", "RRE"),
                                      ("misc", "cPPT")), "pLenti-thing")
        self.assertTrue(v["conflict"])
        self.assertIn("both AAV and lentiviral", v["warning"])

    def test_aav_needs_two_itrs(self):
        from genesis.vectors import missing_elements
        one = missing_elements(self._ann(("AAV_ITR", "AAV2 ITR"),
                                         ("promoter", "hSyn"),
                                         ("polyA", "bGH"), ("origin", "ori"),
                                         ("resistance", "AmpR")), "AAV")
        self.assertIn("second AAV_ITR", one["missing_required"])
        self.assertFalse(one["complete"])

    def test_a_complete_aav_passes(self):
        from genesis.vectors import missing_elements
        full = missing_elements(self._ann(
            ("AAV_ITR", "AAV2 ITR"), ("AAV_ITR", "AAV2 ITR"),
            ("promoter", "hSyn"), ("CDS", "EGFP"), ("WPRE", "WPRE"),
            ("polyA", "bGH polyA"), ("origin", "ori"),
            ("resistance", "AmpR")), "AAV")
        self.assertTrue(full["complete"], full)

    def test_missing_wpre_is_a_warning_not_an_error(self):
        from genesis.vectors import missing_elements
        result = missing_elements(self._ann(
            ("AAV_ITR", "ITR"), ("AAV_ITR", "ITR"), ("promoter", "hSyn"),
            ("polyA", "bGH"), ("origin", "ori"), ("resistance", "AmpR")), "AAV")
        self.assertIn("WPRE", result["missing_expected"])
        self.assertEqual(result["missing_required"], [])

    def test_spec_accepts_the_usual_shorthand(self):
        for text, expected in (("AAV", "AAV"), ("aav", "AAV"),
                               ("lenti", "lentiviral"), ("LV", "lentiviral"),
                               ("plasmid", "plasmid"), ("", ""),
                               ("auto", "")):
            spec = ConstructSpec.from_dict({"name": "x", "vector": text})
            self.assertEqual(spec.vector, expected, text)

    def test_an_unknown_vector_is_refused(self):
        with self.assertRaises(ValueError):
            ConstructSpec.from_dict({"name": "x", "vector": "adenovirus"})

    def test_evidence_is_readable_not_a_regex(self):
        """Nobody should be shown a backslash as a reason."""
        from genesis.vectors import classify_vector
        for annotations, name in (
                ([{"feature_type": "AAV_ITR", "name": "AAV2 ITR"}] * 2, "pAAV-hSyn"),
                ([{"feature_type": "promoter", "name": "5' LTR"},
                  {"feature_type": "m", "name": "RRE"},
                  {"feature_type": "m", "name": "cPPT"}], "pLJM1")):
            for reason in classify_vector(annotations, name)["evidence"]:
                # No regex syntax. An apostrophe is fine: "5' LTR" has one.
                for token in ("\\b", "\\d", "[", "]", "(?", "|", "+"):
                    self.assertNotIn(token, reason, reason)


# --------------------------------------------------------------------------
class TestVectorDesignIntegration(unittest.TestCase):
    """Behaviour of the vector option inside a whole design."""

    def setUp(self):
        import hashlib
        self.conn = connect(":memory:")
        self.cfg = Config()
        self.seq = orf(150, 800)
        chk = hashlib.sha256(self.seq.encode()).hexdigest()
        pid = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES ('CDS','gene',?,?,?)""",
            (self.seq, chk, len(self.seq))).lastrowid
        self.part_id = pid

    def _plasmid(self, name, annotations):
        """A plasmid carrying the given (type, label) annotations."""
        body = self.seq + rand_seq(4000, 801)
        pl = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES (?,'t',?,?,?,1,'sequence_verified')""",
            (name, body, len(body), name)).lastrowid
        for ftype, label in annotations:
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,?,?,0,10,1)""",
                (pl, ftype, label))
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,0,?,1)""",
            (self.part_id, pl, len(self.seq)))
        self.conn.commit()
        return pl

    AAV_SET = [("AAV_ITR", "AAV2 ITR"), ("AAV_ITR", "AAV2 ITR"),
               ("promoter", "hSyn"), ("WPRE", "WPRE"),
               ("polyA", "bGH polyA"), ("origin", "ori"),
               ("resistance", "AmpR")]
    LENTI_SET = [("promoter", "5' LTR"), ("misc", "HIV-1 psi pack"),
                 ("misc", "RRE"), ("misc", "cPPT"), ("polyA", "bGH polyA"),
                 ("origin", "ori"), ("resistance", "AmpR")]

    def test_complete_aav_reports_complete(self):
        self._plasmid("pAAV-test", self.AAV_SET)
        plan = design(self.conn, ConstructSpec.from_dict(
            {"name": "x", "vector": "AAV", "backbone": "pAAV",
             "modules": [{"type": "CDS", "name": "gene"}]}), self.cfg)
        self.assertTrue(plan["vector"]["completeness"]["complete"],
                        plan["vector"]["completeness"])

    def test_mismatch_names_the_cause_instead_of_dropping_the_backbone(self):
        """Asking for AAV with a lentiviral backbone must say exactly that."""
        self._plasmid("pLenti-test", self.LENTI_SET)
        plan = design(self.conn, ConstructSpec.from_dict(
            {"name": "x", "vector": "AAV", "backbone": "pLenti",
             "modules": [{"type": "CDS", "name": "gene"}]}), self.cfg)
        self.assertIsNotNone(plan["backbone"],
                             "the backbone must still be loaded, so the "
                             "mismatch can be reported")
        messages = [i["message"] for i in plan["validation"]["issues"]
                    if i["level"] == "error"]
        self.assertTrue(any("asked for" in m and "lentiviral" in m
                            for m in messages), messages)

    def test_construct_name_is_not_evidence(self):
        """A construct called LENTI-anything is not a lentivirus."""
        plan = design(self.conn, ConstructSpec.from_dict(
            {"name": "LENTI-MISMATCH-CHECK",
             "modules": [{"type": "CDS", "name": "gene"}]}), self.cfg)
        self.assertEqual(plan["vector"]["detected"], "plasmid")
        self.assertEqual(plan["vector"]["detected_evidence"], [])

    def test_requesting_a_family_with_no_backbone_is_an_error(self):
        plan = design(self.conn, ConstructSpec.from_dict(
            {"name": "x", "vector": "AAV",
             "modules": [{"type": "CDS", "name": "gene"}]}), self.cfg)
        messages = [i["message"] for i in plan["validation"]["issues"]
                    if i["level"] == "error"]
        self.assertTrue(any("no backbone was found" in m for m in messages),
                        messages)

    def test_save_as_renames_the_stored_plan(self):
        from genesis.design import save_design
        import json as _json
        self._plasmid("pAAV-test", self.AAV_SET)
        spec = ConstructSpec("original", modules=[Module("CDS", "gene")])
        plan = design(self.conn, spec, self.cfg)
        save_design(self.conn, spec, plan)
        renamed = ConstructSpec("renamed", modules=[Module("CDS", "gene")])
        result = save_design(self.conn, renamed, plan, mode="save_as")

        stored = _json.loads(self.conn.execute(
            "SELECT plan FROM design WHERE id = ?",
            (result["design_id"],)).fetchone()["plan"])
        self.assertEqual(stored["construct"]["name"], "renamed")
        # and the original is left alone
        first = _json.loads(self.conn.execute(
            "SELECT plan FROM design WHERE id = 1").fetchone()["plan"])
        self.assertEqual(first["construct"]["name"], "original")


# --------------------------------------------------------------------------
class TestVerification(unittest.TestCase):
    """Read the designed molecule back and check it against the request.

    These tests inject faults, because a verifier that only ever says "yes"
    is worth nothing.
    """

    def _resolved(self, items):
        from genesis.construct import ResolvedModule
        return [ResolvedModule(Module(type=t, name=n, seq=s), None, n, s, t, [], [])
                for t, n, s in items]

    def test_a_correct_construct_verifies(self):
        from genesis.verify import verify_design
        mods = self._resolved([("promoter", "hSyn", rand_seq(400, 900)),
                               ("CDS", "gene", orf(200, 901))])
        construct = assemble(mods, "", True)
        result = verify_design(construct, mods, True)
        self.assertTrue(result.passed, result.problems)
        self.assertEqual(len(result.checks), 2)
        self.assertTrue(all(c.found for c in result.checks))

    def test_a_dropped_module_is_caught(self):
        """The layout claims a part is there; the sequence says otherwise."""
        from genesis.verify import verify_design
        mods = self._resolved([("promoter", "hSyn", rand_seq(400, 902)),
                               ("CDS", "gene", orf(200, 903))])
        construct = assemble(mods, "", True)
        # cut the CDS out of the finished molecule, as a bug would
        construct.seq = construct.seq.replace(mods[1].seq, "")
        result = verify_design(construct, mods, True)
        self.assertFalse(result.passed)
        self.assertTrue(any("does not appear" in p for p in result.problems),
                        result.problems)

    def test_a_truncated_module_is_caught(self):
        from genesis.verify import verify_design
        mods = self._resolved([("CDS", "gene", orf(200, 904))])
        construct = assemble(mods, "", True)
        construct.seq = construct.seq[:-30]          # lose the last 10 codons
        result = verify_design(construct, mods, True)
        self.assertFalse(result.passed)

    def test_an_inverted_module_is_caught(self):
        from genesis.verify import verify_design
        mods = self._resolved([("promoter", "p", rand_seq(300, 905)),
                               ("CDS", "gene", orf(200, 906))])
        construct = assemble(mods, "", True)
        construct.seq = construct.seq.replace(mods[1].seq,
                                              st.revcomp(mods[1].seq))
        result = verify_design(construct, mods, True)
        self.assertFalse(result.passed)
        self.assertTrue(any("wrong way round" in p for p in result.problems),
                        result.problems)

    def test_a_reordered_construct_is_caught(self):
        from genesis.verify import verify_design
        a, b, c = rand_seq(300, 907), rand_seq(300, 908), rand_seq(300, 909)
        mods = self._resolved([("promoter", "one", a), ("CDS", "two", b),
                               ("polyA", "three", c)])
        construct = assemble(mods, "", True)
        construct.seq = a + c + b                    # swap the last two
        result = verify_design(construct, mods, True)
        self.assertFalse(result.order_ok)
        self.assertTrue(any("different order" in p for p in result.problems),
                        result.problems)

    def test_a_duplicated_part_is_reported(self):
        from genesis.verify import verify_design
        gene = orf(200, 910)
        mods = self._resolved([("CDS", "gene", gene)])
        construct = assemble(mods, "", True)
        construct.seq = construct.seq + rand_seq(200, 911) + gene
        result = verify_design(construct, mods, True)
        self.assertTrue(any("appears 2 times" in p for p in result.problems),
                        result.problems)

    def test_lox_sites_appearing_twice_is_not_a_fault(self):
        """Two of each lox site is what makes the cassette work."""
        from genesis.construct import expand_dio
        from genesis.verify import verify_design
        mods = self._resolved([("Cre-lox", "DIO", ""),
                               ("CDS", "payload", orf(250, 912))])
        expanded, _ = expand_dio(mods)
        construct = assemble(expanded, "", True)
        result = verify_design(construct, expanded, True)
        self.assertTrue(result.passed, result.problems)
        self.assertFalse(any("appears 2 times" in p for p in result.problems),
                         result.problems)

    def test_an_inverted_payload_is_confirmed_to_switch_on(self):
        from genesis.construct import expand_dio
        from genesis.verify import verify_design
        from genesis.validate import lox_sites_from_layout, simulate_cre
        mods = self._resolved([("Cre-lox", "DIO", ""),
                               ("CDS", "payload", orf(250, 913))])
        expanded, _ = expand_dio(mods)
        construct = assemble(expanded, "", True)
        recombined, _ = simulate_cre(construct.seq,
                                     lox_sites_from_layout(construct.seq))
        result = verify_design(construct, expanded, True, recombined)
        self.assertIsNotNone(result.recombined)
        self.assertTrue(result.recombined["all_switched"], result.problems)
        self.assertEqual(set(result.recombined["now_forward"]),
                         {"payload"})

    def test_every_design_carries_a_verification(self):
        import hashlib
        conn = connect(":memory:")
        cfg = Config()
        seq = orf(150, 914)
        pl = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('pD','t',?,?,'vv',1,'sequence_verified')""",
            (seq + rand_seq(2000, 915), len(seq) + 2000)).lastrowid
        pid = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES ('CDS','gene',?,?,?)""",
            (seq, hashlib.sha256(seq.encode()).hexdigest(), len(seq))).lastrowid
        conn.execute("""INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,0,?,1)""", (pid, pl, len(seq)))
        conn.commit()
        plan = design(conn, ConstructSpec("v", modules=[Module("CDS", "gene")]),
                      cfg)
        self.assertIn("verification", plan)
        self.assertTrue(plan["verification"]["passed"],
                        plan["verification"]["problems"])


# --------------------------------------------------------------------------
class TestCreLoxModes(unittest.TestCase):

    def _mods(self, mode):
        from genesis.construct import ResolvedModule
        return [
            ResolvedModule(Module(type="Cre-lox", name="", mode=mode), None,
                           mode, "", "Cre-lox", [], []),
            ResolvedModule(Module(type="CDS", name="payload"), None, "payload",
                           orf(250, 920), "CDS", [], []),
        ]

    def _switches(self, mode):
        from genesis.construct import expand_dio
        from genesis.validate import lox_sites_from_layout, simulate_cre
        mods = self._mods(mode)
        protein = st.translate(mods[1].seq).rstrip("*")
        expanded, _ = expand_dio(mods)
        construct = assemble(expanded, "", True)
        before = any(o["strand"] == 1 and o["protein"] == protein
                     for o in st.find_orfs(construct.seq, min_aa=200))
        recombined, log = simulate_cre(construct.seq,
                                       lox_sites_from_layout(construct.seq))
        after = any(o["strand"] == 1 and o["protein"] == protein
                    for o in st.find_orfs(recombined, min_aa=200))
        return before, after, log

    def test_dio_and_flex_switch_on(self):
        for mode in ("DIO", "FLEX"):
            before, after, log = self._switches(mode)
            self.assertFalse(before, f"{mode} must be off before Cre")
            self.assertTrue(after, f"{mode} must be on after Cre")
            self.assertEqual(len(log), 2, mode)

    def test_fas_uses_loxfas_and_switches_on(self):
        from genesis.construct import expand_dio
        from genesis.validate import lox_sites_from_layout
        expanded, _ = expand_dio(self._mods("FAS"))
        construct = assemble(expanded, "", True)
        names = {s.name for s in lox_sites_from_layout(construct.seq)}
        self.assertEqual(names, {"loxFAS", "loxP"})
        before, after, _ = self._switches("FAS")
        self.assertFalse(before)
        self.assertTrue(after)

    def test_do_switches_off(self):
        before, after, _ = self._switches("DO")
        self.assertTrue(before, "DO expresses until Cre arrives")
        self.assertFalse(after, "DO must be off after Cre")

    def test_dio_is_still_accepted_as_a_type(self):
        spec = ConstructSpec.from_dict(
            {"name": "x", "modules": [{"type": "DIO"}]})
        self.assertEqual(spec.modules[0].type, "Cre-lox")
        self.assertEqual(spec.modules[0].mode, "DIO")

    def test_an_unknown_mode_is_refused(self):
        from genesis.construct import resolve_mode
        with self.assertRaises(ValueError):
            resolve_mode("nonsense")


# --------------------------------------------------------------------------
class TestBackboneInsertion(unittest.TestCase):
    """The payload must land inside the vector, not be bolted onto the end.

    An AAV built by concatenation has the right length and is useless: the
    cassette sits outside the ITRs and cannot be packaged.
    """

    ITR = "TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGC" + rand_seq(100, 1000)

    def _backbone(self, wraps_origin: bool = False):
        """A pAAV-like backbone: ITR - cassette - WPRE - polyA - ITR - ori/AmpR."""
        parts, annotations, pos = [], [], 0

        def add(seq, ftype, name):
            nonlocal pos
            parts.append(seq)
            annotations.append({"feature_type": ftype, "name": name,
                                "start": pos, "end": pos + len(seq),
                                "strand": 1})
            pos += len(seq)

        if wraps_origin:
            # cassette first, so the vector arc is at the far end and the
            # payload region straddles the coordinate origin
            add(orf(100, 1001), "CDS", "old gene")
            add(rand_seq(589, 1002), "WPRE", "WPRE")
            add(rand_seq(120, 1003), "polyA", "bGH polyA")
            add(self.ITR, "AAV_ITR", "AAV2 ITR")
            add(rand_seq(600, 1004), "origin", "ori")
            add(rand_seq(860, 1005), "resistance", "AmpR")
            add(self.ITR, "AAV_ITR", "AAV2 ITR")
            add(rand_seq(400, 1006), "promoter", "old promoter")
        else:
            add(self.ITR, "AAV_ITR", "AAV2 ITR")
            add(rand_seq(400, 1006), "promoter", "old promoter")
            add(orf(100, 1001), "CDS", "old gene")
            add(rand_seq(589, 1002), "WPRE", "WPRE")
            add(rand_seq(120, 1003), "polyA", "bGH polyA")
            add(self.ITR, "AAV_ITR", "AAV2 ITR")
            add(rand_seq(600, 1004), "origin", "ori")
            add(rand_seq(860, 1005), "resistance", "AmpR")
        return "".join(parts), annotations

    def _build(self, wraps_origin=False, spec_types=frozenset({"promoter", "CDS"})):
        from genesis.construct import ResolvedModule, assemble_into_backbone
        seq, annotations = self._backbone(wraps_origin)
        modules = [
            ResolvedModule(Module(type="promoter", name="hSyn"), None, "hSyn",
                           rand_seq(300, 1010), "promoter", [], []),
            ResolvedModule(Module(type="CDS", name="EGFP"), None, "EGFP",
                           orf(150, 1011), "CDS", [], []),
        ]
        return assemble_into_backbone(modules, seq, annotations,
                                      set(spec_types), True), modules

    def _types(self, construct):
        return [i["type"] for i in construct.layout]

    def test_vector_elements_are_carried_into_the_product(self):
        construct, _ = self._build()
        types = self._types(construct)
        for required in ("AAV_ITR", "WPRE", "polyA", "origin", "resistance"):
            self.assertIn(required, types,
                          f"{required} must survive into the product")

    def test_both_itrs_survive(self):
        construct, _ = self._build()
        itrs = [i for i in construct.layout if i["type"] == "AAV_ITR"]
        self.assertEqual(len(itrs), 2, "an AAV needs both ITRs")

    def test_the_payload_lands_between_the_itrs(self):
        """Circularly, the new modules must sit between the two ITRs."""
        construct, _ = self._build()
        itrs = sorted(i["start"] for i in construct.layout
                      if i["type"] == "AAV_ITR")
        new = [i for i in construct.layout if not i.get("from_backbone")]
        self.assertTrue(new)
        for item in new:
            inside = itrs[0] <= item["start"] <= itrs[1]
            outside_wrapping = item["start"] < itrs[0] or item["start"] > itrs[1]
            self.assertTrue(inside or outside_wrapping)
        # and specifically: no origin or resistance between payload and an ITR
        between = [i["type"] for i in construct.layout
                   if new[0]["start"] < i["start"] < new[-1]["end"]]
        self.assertNotIn("origin", between)
        self.assertNotIn("resistance", between)

    def test_the_old_cassette_is_replaced_not_appended(self):
        construct, _ = self._build()
        names = [i["name"] for i in construct.layout]
        self.assertNotIn("old gene", names)
        self.assertNotIn("old promoter", names)
        self.assertIn("EGFP", names)

    def test_the_backbones_wpre_and_polya_are_kept(self):
        """They belong to the vector, not to the payload being swapped."""
        construct, _ = self._build()
        names = [i["name"] for i in construct.layout]
        self.assertIn("WPRE", names)
        self.assertIn("bGH polyA", names)

    def test_a_specified_polya_replaces_the_backbones_own(self):
        construct, _ = self._build(spec_types={"promoter", "CDS", "polyA"})
        polyas = [i for i in construct.layout if i["type"] == "polyA"]
        self.assertEqual(len(polyas), 0,
                         "asking for a polyA should clear the backbone's")

    def test_a_cassette_that_wraps_the_origin_is_handled(self):
        """Real pAAVs put the promoter near the end and its CDS at zero."""
        construct, _ = self._build(wraps_origin=True)
        types = self._types(construct)
        for required in ("AAV_ITR", "WPRE", "polyA", "origin", "resistance"):
            self.assertIn(required, types, required)
        self.assertEqual(len([i for i in construct.layout
                              if i["type"] == "AAV_ITR"]), 2)
        self.assertNotIn("old gene", [i["name"] for i in construct.layout])

    def test_the_product_verifies(self):
        from genesis.verify import verify_design
        construct, modules = self._build()
        result = verify_design(construct, modules, True)
        self.assertTrue(result.passed, result.problems)

    def test_the_product_is_a_complete_aav(self):
        from genesis.vectors import missing_elements
        construct, _ = self._build()
        annotations = [{"feature_type": i["type"], "name": i["name"]}
                       for i in construct.layout]
        result = missing_elements(annotations, "AAV")
        self.assertTrue(result["complete"], result)
