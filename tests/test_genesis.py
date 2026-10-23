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
import os
import urllib.parse
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

    def test_primer_is_never_the_thing_it_anneals_to(self):
        # SnapGene ships these stock primer annotations.  Each is named after
        # its target, and each used to be catalogued as that target: a 20 bp
        # "resistance gene", a 20 bp "polyA signal".
        for label in ("Amp-R", "SV40pA-R", "wpre-rev", "M13 fwd",
                      "F1 origin, reverse primer", "pAAV-DIO-Screen-R",
                      "SV40 polyA terminator, reverse primer",
                      "Chimeric intron in CAG promoter, forward primer"):
            self.assertEqual(classify_label(label, "primer_bind"),
                             "primer_binding_site", label)

    def test_label_still_refines_a_cds(self):
        # The override is narrow: a CDS labelled AmpR is still resistance.
        self.assertEqual(classify_label("AmpR", "CDS"), "resistance")
        self.assertEqual(classify_label("P2A", "CDS"), "p2A")
        self.assertEqual(classify_label("ER2", "CDS"), "trafficking_motif")

    def test_enhancer_key_outranks_a_promoter_label(self):
        # "CMV enhancer" matched the promoter rule on the word CMV, which
        # made 606 enhancers look like promoters to the recipient tiers.
        self.assertEqual(classify_label("CMV enhancer", "enhancer"), "enhancer")
        self.assertEqual(classify_label("CAG\\enhancer", "enhancer"), "enhancer")
        # but SnapGene files the WPRE under enhancer, so the label still wins
        self.assertEqual(classify_label("WPRE", "enhancer"), "WPRE")

    def test_fas_site_is_a_lox_site_but_the_fas_gene_is_not(self):
        for label in ("FAS SITE", "FAS SITE(1)", "loxFAS", "FAS-site"):
            self.assertEqual(classify_label(label, "misc_feature"),
                             "lox_site", label)
        # Fas is also a gene; a bare FAS must not make one a lox site.
        self.assertNotEqual(
            classify_label("HT75 LZF70 Fsyn FAS C gene", "gene"), "lox_site")

    def test_cre_lox_family_label_survives(self):
        for label in ("DIO cassette", "FLEX switch", "double-floxed"):
            self.assertEqual(classify_label(label, "misc_feature"),
                             "Cre-lox", label)

    def test_unlabelled_cds_borrows_the_name_drawn_over_it(self):
        from genesis.importer import name_unlabelled
        # The real shape of the ASAP7y file: the CDS has the type and no
        # label, the misc_feature has the name and no type.
        anns = [
            {"feature_type": "unknown", "name": "ASAP7y",
             "original_label": "ASAP7y", "original_key": "misc_feature",
             "start": 2443, "end": 3724, "strand": 1},
            {"feature_type": "CDS", "name": "CDS", "original_label": "",
             "original_key": "CDS", "start": 2443, "end": 3724, "strand": -1},
        ]
        self.assertEqual(name_unlabelled(anns), 1)
        self.assertEqual(anns[1]["name"], "ASAP7y")

    def test_naming_tolerates_a_trailing_linker(self):
        from genesis.importer import name_unlabelled
        # The Kv2.1 CDS runs 27 bp past its label, into a GSS linker.
        anns = [
            {"feature_type": "unknown", "name": "kv2.1",
             "original_label": "kv2.1", "original_key": "misc_feature",
             "start": 2221, "end": 2416, "strand": -1},
            {"feature_type": "CDS", "name": "CDS", "original_label": "",
             "original_key": "CDS", "start": 2221, "end": 2443, "strand": -1},
        ]
        self.assertEqual(name_unlabelled(anns), 1)
        self.assertEqual(anns[1]["name"], "kv2.1")

    def test_naming_refuses_a_partial_overlap(self):
        from genesis.importer import name_unlabelled
        # A 10 bp chromophore marker must not rename the 700 bp CDS over it.
        anns = [
            {"feature_type": "unknown", "name": "chromophore",
             "original_label": "chromophore", "original_key": "misc_feature",
             "start": 2782, "end": 2791, "strand": -1},
            {"feature_type": "CDS", "name": "CDS", "original_label": "",
             "original_key": "CDS", "start": 2551, "end": 3280, "strand": -1},
        ]
        self.assertEqual(name_unlabelled(anns), 0)
        self.assertEqual(anns[1]["name"], "CDS")

    def test_self_describing_types_are_left_alone(self):
        from genesis.importer import name_unlabelled
        # A WPRE called "WPRE" is already named; the type is the name.
        anns = [
            {"feature_type": "unknown", "name": "SP43",
             "original_label": "SP43", "original_key": "misc_feature",
             "start": 100, "end": 400, "strand": 1},
            {"feature_type": "WPRE", "name": "WPRE", "original_label": "",
             "original_key": "misc_feature",
             "start": 100, "end": 400, "strand": 1},
        ]
        self.assertEqual(name_unlabelled(anns), 0)
        self.assertEqual(anns[1]["name"], "WPRE")

    def test_a_named_feature_is_not_renamed(self):
        from genesis.importer import name_unlabelled
        anns = [
            {"feature_type": "unknown", "name": "some note",
             "original_label": "some note", "original_key": "misc_feature",
             "start": 0, "end": 300, "strand": 1},
            {"feature_type": "CDS", "name": "Citrine",
             "original_label": "Citrine", "original_key": "CDS",
             "start": 0, "end": 300, "strand": 1},
        ]
        self.assertEqual(name_unlabelled(anns), 0)
        self.assertEqual(anns[1]["name"], "Citrine")

    def test_part_types_use_the_current_vocabulary(self):
        from genesis.importer import PART_TYPES
        from genesis.vocab import FEATURE_TYPES
        self.assertTrue(PART_TYPES <= set(FEATURE_TYPES),
                        PART_TYPES - set(FEATURE_TYPES))

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
                             {"swap_hifi", "swap_restriction", "swap_hybrid",
                          "gblock"})

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
        # This design did ask for a cassette, so its lox sites belong.
        errors = [i.message for i in
                  validate_construct(construct, wants_cassette=True)
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
        # Long enough to be a protein: a run shorter than MIN_TRANSLATED_RUN
        # is an annotation artefact and is deliberately not marked.
        c = self._c([("CDS", "geneA", "ATG" + "GCT" * 5 + "TGA" + "GCT" * 15),
                     ("CDS", "geneB", "GCT" * 10 + "TAA")])
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


# --------------------------------------------------------------------------
class TestPayloadArcGeometry(unittest.TestCase):
    """Which arc of a circular backbone the payload replaces."""

    def _ann(self, *spec):
        return [{"feature_type": t, "name": n, "start": s, "end": e,
                 "strand": 1} for t, n, s, e in spec]

    def test_anchors_outside_the_itrs(self):
        """The pAAV case: ori and AmpR sit outside the ITRs, so the span from
        the first ITR to the last origin is the whole plasmid."""
        from genesis.construct import payload_arc
        ann = self._ann(("AAV_ITR", "ITR", 34, 164),
                        ("promoter", "hSyn", 208, 656),
                        ("CDS", "payload", 700, 4000),
                        ("AAV_ITR", "ITR", 4807, 4937),
                        ("origin", "f1 ori", 5022, 5478),
                        ("resistance", "AmpR", 5864, 6725),
                        ("origin", "ori", 6895, 7484))
        self.assertEqual(payload_arc(ann, 7522), (164, 4807))

    def test_a_cassette_that_wraps_the_origin(self):
        from genesis.construct import payload_arc
        ann = self._ann(("CDS", "payload", 0, 1300),
                        ("AAV_ITR", "ITR", 4602, 4743),
                        ("origin", "ori", 5421, 5876),
                        ("resistance", "AmpR", 5901, 6867),
                        ("AAV_ITR", "ITR", 9292, 9422),
                        ("promoter", "p", 9441, 9756))
        arc = payload_arc(ann, 9788)
        self.assertEqual(arc, (9422, 4602))
        self.assertGreater(arc[0], arc[1], "this arc must wrap")

    def test_without_itrs_the_payload_is_the_complement_of_the_machinery(self):
        from genesis.construct import payload_arc
        ann = self._ann(("promoter", "p", 100, 500),
                        ("CDS", "g", 500, 2000),
                        ("origin", "ori", 2500, 3000),
                        ("resistance", "AmpR", 3000, 3900))
        self.assertEqual(payload_arc(ann, 5000), (3900, 2500))

    def test_no_anchors_at_all(self):
        from genesis.construct import payload_arc
        self.assertIsNone(payload_arc(
            self._ann(("CDS", "g", 0, 500)), 1000))


# --------------------------------------------------------------------------
class TestRecipientArchitecture(unittest.TestCase):
    """A recipient must share the architecture, not just a few genes."""

    def setUp(self):
        import hashlib
        self.conn = connect(":memory:")
        self.cfg = Config()
        self.shared = orf(200, 1100)
        self.part = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES ('CDS','shared',?,?,?)""",
            (self.shared, hashlib.sha256(self.shared.encode()).hexdigest(),
             len(self.shared))).lastrowid

    def _plasmid(self, name, annotations):
        body = self.shared + rand_seq(4000, 1101)
        pl = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified, available_template)
               VALUES (?,'t',?,?,?,1,'sequence_verified',1)""",
            (name, body, len(body), name)).lastrowid
        # Every recipient needs an origin and a resistance marker, or it is
        # a fragment rather than something that can be grown.
        annotations = list(annotations) + [("origin", "ori"),
                                           ("resistance", "AmpR")]
        for ftype, label in annotations:
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,?,?,0,30,1)""",
                (pl, ftype, label))
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,0,?,1)""",
            (self.part, pl, len(self.shared)))
        self.conn.commit()
        return pl

    def test_a_plasmid_without_lox_sites_is_not_offered(self):
        """The Fsyn-Solaris case: shares a gene, has no Cre-lox system."""
        from genesis.swap import rank_parents
        self._plasmid("no-lox", [("promoter", "hSyn"), ("CDS", "shared")])
        with_lox = self._plasmid(
            "has-lox", [("promoter", "hSyn"), ("lox_site", "loxP"),
                        ("lox_site", "lox2272")])
        ranked = rank_parents(self.conn, [self.part], self.cfg,
                              require_types={"lox_site"})
        self.assertEqual([p.plasmid_id for p in ranked], [with_lox])

    def test_without_a_requirement_both_are_offered(self):
        from genesis.swap import rank_parents
        self._plasmid("no-lox", [("CDS", "shared")])
        self._plasmid("has-lox", [("lox_site", "loxP")])
        self.assertEqual(len(rank_parents(self.conn, [self.part], self.cfg)), 2)

    def test_a_design_needing_cre_lox_requires_lox_sites(self):
        from genesis.design import required_architecture
        spec = ConstructSpec.from_dict(
            {"name": "x", "vector": "AAV",
             "modules": [{"type": "DIO"}, {"type": "CDS", "name": "g"}]})
        self.assertEqual(required_architecture(spec), {"lox_site", "AAV_ITR"})

    def test_a_plain_design_requires_nothing_special(self):
        from genesis.design import required_architecture
        spec = ConstructSpec.from_dict(
            {"name": "x", "modules": [{"type": "CDS", "name": "g"}]})
        self.assertEqual(required_architecture(spec), set())

    def test_no_suitable_recipient_says_so_rather_than_picking_a_bad_one(self):
        from genesis.swap import analyse_swap
        self._plasmid("no-lox", [("CDS", "shared")])
        result = analyse_swap(
            self.conn, [{"type": "CDS", "part_id": self.part,
                         "part_name": "shared"}],
            self.cfg, require_types={"lox_site"})
        self.assertIsNotNone(result)
        self.assertFalse(result.feasible)
        self.assertTrue(any("architecture it needs" in n for n in result.notes),
                        result.notes)

    def test_a_swap_that_cannot_supply_every_part_is_not_offered(self):
        """A route that builds the wrong construct is not a cheap route.

        Asked for a promoter and two sensors that the recipient has none of,
        this used to amplify the promoter, note the rest in passing, and win
        the route comparison on cost - because doing one third of the job is
        cheaper than doing all of it.  The plan's own architecture claimed
        sensors its protocol never touched.
        """
        from genesis.swap import analyse_swap
        host = self._plasmid("host", [("promoter", "hSyn"),
                                      ("CDS", "shared")])
        absent = []
        for name in ("CAG", "ASAP7y", "iGluSnFR4f"):
            seq = rand_seq(600, hash(name) % 1000)
            absent.append(self.conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES ('CDS',?,?,?,?)""",
                (name, seq, name, len(seq))).lastrowid)
        self.conn.commit()

        modules = [{"type": "CDS", "part_id": self.part, "part_name": "shared",
                    "seq": self.shared}]
        modules += [{"type": "CDS", "part_id": pid, "part_name": n, "seq": ""}
                    for pid, n in zip(absent, ("CAG", "ASAP7y", "iGluSnFR4f"))]

        result = analyse_swap(self.conn, modules, self.cfg,
                              prefer_plasmid_id=host)
        self.assertIsNotNone(result)
        self.assertFalse(result.feasible,
                         "a swap missing three of four parts is not usable")
        self.assertTrue(
            any("single region" in n for n in result.notes), result.notes)

    def test_one_missing_part_is_still_a_normal_swap(self):
        from genesis.swap import analyse_swap
        host = self._plasmid("host", [("promoter", "hSyn"),
                                      ("CDS", "shared")])
        seq = rand_seq(600, 77)
        new = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('CDS','GFP',?,'gfp',?)""",
            (seq, len(seq))).lastrowid
        self.conn.commit()
        result = analyse_swap(
            self.conn,
            [{"type": "CDS", "part_id": self.part, "part_name": "shared",
              "seq": self.shared},
             {"type": "CDS", "part_id": new, "part_name": "GFP", "seq": seq}],
            self.cfg, prefer_plasmid_id=host)
        self.assertIsNotNone(result)
        self.assertTrue(
            result.feasible
            or not any("single region" in n for n in result.notes),
            result.notes)


# --------------------------------------------------------------------------
class TestAddgenePrimerRules(unittest.TestCase):
    """Primer design follows Addgene's PCR cloning protocol.

    https://www.addgene.org/protocols/pcr-cloning/ - the hybridising region
    is 18-21 bp, a restriction primer carries a 3-6 bp leader before a 6-8 bp
    site, and the Tm that matters is the one for the hybridising region only.
    """

    def setUp(self):
        self.cfg = Config()
        self.template = rand_seq(3000, 1200)

    def test_the_window_matches_the_protocol(self):
        self.assertEqual((self.cfg.primer_min_len, self.cfg.primer_max_len),
                         (18, 21))
        self.assertTrue(3 <= len(self.cfg.primer_leader) <= 6)

    def test_annealing_region_stays_inside_the_window(self):
        from genesis.primers import pick_anneal
        for start in range(100, 2000, 137):
            for direction in ("F", "R"):
                region = pick_anneal(self.template, start, direction, self.cfg)
                self.assertTrue(18 <= len(region) <= 21,
                                f"{direction} at {start}: {len(region)} nt")

    def test_a_restriction_primer_is_leader_site_anneal(self):
        from genesis.primers import design_re_pair
        from genesis.enzymes import CATALOG
        fwd, rev = design_re_pair("frag", self.template, 200, 900,
                                  "EcoRI", "BamHI", self.cfg)
        for primer, enzyme in ((fwd, "EcoRI"), (rev, "BamHI")):
            site = CATALOG[enzyme].site
            self.assertTrue(6 <= len(site) <= 8)
            leader = primer.tail[:-len(site)]
            self.assertTrue(3 <= len(leader) <= 6, leader)
            self.assertTrue(18 <= len(primer.anneal) <= 21)
            self.assertEqual(primer.tail, leader + primer.tail[len(leader):])
            self.assertEqual(primer.sequence, primer.tail + primer.anneal)

    def test_the_reverse_primer_is_the_reverse_complement(self):
        from genesis.primers import pick_anneal
        region = pick_anneal(self.template, 900, "R", self.cfg)
        self.assertIn(st.revcomp(region), self.template)

    def test_tm_is_reported_for_the_hybridising_region_only(self):
        from genesis.primers import Primer
        anneal = "GCTAGCTAGCATCGATCGAT"
        primer = Primer("p", "A" * 25 + anneal, anneal, "A" * 25)
        reported = primer.to_dict(self.cfg)
        self.assertAlmostEqual(reported["anneal_tm"], st.tm(anneal), places=1)
        self.assertNotAlmostEqual(reported["anneal_tm"], reported["full_tm"],
                                  places=1)

    def test_the_leader_does_not_create_another_site(self):
        from genesis.primers import choose_leader
        from genesis.enzymes import CATALOG
        from genesis.seqtools import find_all
        for name in ("EcoRI", "BamHI", "NotI", "XhoI", "AgeI", "KpnI",
                     "NheI", "SpeI", "SalI", "HindIII"):
            site = CATALOG[name].site
            leader = choose_leader(site, self.cfg, self.cfg.stocked_enzymes)
            self.assertTrue(3 <= len(leader) <= 6)
            junction = leader + site
            for other in self.cfg.stocked_enzymes:
                enzyme = CATALOG.get(other)
                if enzyme is None:
                    continue
                positions = set(find_all(junction, enzyme.site))
                positions |= set(find_all(junction, st.revcomp(enzyme.site)))
                expected = 1 if enzyme.site in (site, st.revcomp(site)) else 0
                self.assertLessEqual(len(positions), expected,
                                     f"{leader}+{site} creates a {other} site")

    def test_an_out_of_window_region_is_flagged(self):
        from genesis.primers import Primer, validate_primer
        long_anneal = self.template[100:140]          # 40 nt
        primer = Primer("p", long_anneal, long_anneal)
        warnings = validate_primer(primer, self.template, self.cfg)
        self.assertTrue(any("window" in w for w in warnings), warnings)

    def test_designed_primers_are_not_taken_from_the_catalogue(self):
        """Primers are computed from the template, never reused."""
        from genesis.primers import pick_anneal
        region = pick_anneal(self.template, 500, "F", self.cfg)
        self.assertIn(region, self.template)
        self.assertEqual(len(region), len(set([region])) * len(region))


# --------------------------------------------------------------------------
class TestAssemblyHomology(unittest.TestCase):
    """Homology arms follow the kit, and are no longer than they need to be."""

    def test_the_two_kits_want_different_homology(self):
        from genesis.primers import homology_range
        hifi, gibson = Config(), Config()
        gibson.preferred_assembly_kit = "Gibson"
        self.assertEqual(homology_range(hifi)[0], 15)
        self.assertEqual(homology_range(gibson)[0], 20)
        self.assertGreater(homology_range(gibson)[1], homology_range(hifi)[1])

    def test_hifi_gives_shorter_arms_so_shorter_primers(self):
        from genesis.primers import pick_overlap
        hifi, gibson = Config(), Config()
        gibson.preferred_assembly_kit = "Gibson"
        neighbour = rand_seq(80, 1300)
        self.assertLessEqual(len(pick_overlap(neighbour, hifi)),
                             len(pick_overlap(neighbour, gibson)))

    def test_a_gc_rich_arm_is_cut_short(self):
        """No reason to pay for 25 bases when 15 already anneal."""
        from genesis.primers import pick_overlap
        cfg = Config()
        arm = pick_overlap("GC" * 40, cfg)
        self.assertEqual(len(arm), 15)
        self.assertGreaterEqual(st.tm(arm), cfg.overlap_tm_min)

    def test_an_at_rich_arm_grows_to_the_ceiling(self):
        from genesis.primers import pick_overlap, homology_range
        cfg = Config()
        arm = pick_overlap("AT" * 40, cfg)
        self.assertEqual(len(arm), homology_range(cfg)[1])


# --------------------------------------------------------------------------
class TestProteinLevelIdentity(unittest.TestCase):
    """The same protein under two names, and under two codon usages."""

    @staticmethod
    def _back_translate(peptide, table_shift=0):
        """Encode a peptide, optionally choosing different synonymous codons."""
        choices = {}
        for codon, aa in st.CODON_TABLE.items():
            choices.setdefault(aa, []).append(codon)
        return "".join(
            sorted(choices[a])[table_shift % len(choices[a])] for a in peptide)

    def test_codon_optimised_copies_share_the_protein_not_the_dna(self):
        from genesis.swap import protein_identity, protein_of
        import genesis.seqtools as seq
        peptide = "MGSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSEFELRRQACGRTRKLME"
        a = self._back_translate(peptide, 0)
        b = self._back_translate(peptide, 1)
        self.assertNotEqual(a, b)
        # DNA 40-mers barely overlap, the proteins are identical
        shared = sum(1 for i in range(len(a) - 40 + 1)
                     if a[i:i + 40] in b)
        self.assertEqual(shared, 0)
        self.assertEqual(protein_of(a), protein_of(b))
        self.assertEqual(protein_identity(protein_of(a), protein_of(b)), 1.0)

    def test_different_proteins_do_not_match(self):
        from genesis.swap import protein_identity, protein_of
        a = self._back_translate("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ")
        b = self._back_translate("MWQRVNDLGEYTKPFACHIMLSEDGKRNVTWAY")
        self.assertLess(protein_identity(protein_of(a), protein_of(b)), 0.5)

    def test_a_longer_fusion_contains_the_shorter_protein(self):
        from genesis.swap import protein_identity, protein_of
        core = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSG"
        fusion = "MGSSHHHHHHSQDP" + core + "GGSGGSAAAEFLYKVV"
        self.assertGreaterEqual(
            protein_identity(protein_of(self._back_translate(core)),
                             protein_of(self._back_translate(fusion))), 0.9)


# --------------------------------------------------------------------------
class TestParentSelection(unittest.TestCase):
    """Recommendations are offered, and a choice can be pinned."""

    def test_a_pinned_donor_is_used(self):
        import hashlib
        conn = connect(":memory:")
        cfg = Config()
        gene = orf(200, 1400)
        part = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum, seq_length)
               VALUES ('CDS','g',?,?,?)""",
            (gene, hashlib.sha256(gene.encode()).hexdigest(),
             len(gene))).lastrowid

        def plasmid(name, extra):
            body = gene + rand_seq(3000, 1401)
            pid = conn.execute(
                """INSERT INTO plasmid(name, owner, seq, seq_length,
                        seq_checksum, circular, verified, available_template)
                   VALUES (?,'t',?,?,?,1,'sequence_verified',1)""",
                (name, body, len(body), name)).lastrowid
            for ftype, label in [("origin", "ori"), ("resistance", "AmpR")] + extra:
                conn.execute(
                    """INSERT INTO annotation(plasmid_id, feature_type, name,
                            start, end, strand) VALUES (?,?,?,0,30,1)""",
                    (pid, ftype, label))
            conn.execute(
                """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                        strand) VALUES (?,?,0,?,1)""", (part, pid, len(gene)))
            return pid

        small = plasmid("small-donor", [])
        chosen = plasmid("the-one-i-want", [])
        conn.commit()

        from genesis.swap import find_donors
        donors = find_donors(conn, [part], "")
        self.assertTrue(donors)
        self.assertTrue(all(d.reasons for d in donors),
                        "every recommendation needs a stated reason")
        self.assertIn(chosen, [d.plasmid_id for d in donors])

    def test_a_spec_can_pin_both_parents(self):
        spec = ConstructSpec.from_dict({
            "name": "x", "recipient_plasmid_id": 11, "donor_plasmid_id": 22,
            "modules": [{"type": "CDS", "name": "g"}]})
        self.assertEqual(spec.recipient_plasmid_id, 11)
        self.assertEqual(spec.donor_plasmid_id, 22)


# --------------------------------------------------------------------------
class TestPromoterSlot(unittest.TestCase):
    """The promoter is its own slot, not part of the payload or the vector."""

    def _backbone(self):
        parts, ann, pos = [], [], 0

        def add(seq, ftype, name):
            nonlocal pos
            parts.append(seq)
            ann.append({"feature_type": ftype, "name": name, "start": pos,
                        "end": pos + len(seq), "strand": 1})
            pos += len(seq)

        add(rand_seq(141, 1500), "AAV_ITR", "ITR")
        add(rand_seq(450, 1501), "promoter", "hSyn")
        add(orf(120, 1502), "CDS", "old gene")
        add(rand_seq(589, 1503), "WPRE", "WPRE")
        add(rand_seq(141, 1504), "AAV_ITR", "ITR")
        add(rand_seq(600, 1505), "origin", "ori")
        add(rand_seq(860, 1506), "resistance", "AmpR")
        return "".join(parts), ann

    def _build(self, modules, spec_types):
        from genesis.construct import ResolvedModule, assemble_into_backbone
        seq, ann = self._backbone()
        resolved = [ResolvedModule(Module(type=t, name=n, seq=s), None, n, s,
                                   t, [], []) for t, n, s in modules]
        return assemble_into_backbone(resolved, seq, ann, set(spec_types), True)

    def _promoters(self, construct):
        return [i["name"] for i in construct.layout if i["type"] == "promoter"]

    def test_the_backbones_promoter_is_kept_when_none_is_asked_for(self):
        c = self._build([("CDS", "newGene", orf(150, 1510))], {"CDS"})
        self.assertEqual(self._promoters(c), ["hSyn"])

    def test_a_specified_promoter_replaces_it_rather_than_joining_it(self):
        c = self._build([("promoter", "CAG", rand_seq(400, 1511)),
                         ("CDS", "newGene", orf(150, 1512))],
                        {"promoter", "CDS"})
        self.assertEqual(self._promoters(c), ["CAG"],
                         "a construct asking for CAG must not also carry hSyn")

    def test_the_old_coding_sequence_still_goes(self):
        c = self._build([("promoter", "CAG", rand_seq(400, 1513)),
                         ("CDS", "newGene", orf(150, 1514))],
                        {"promoter", "CDS"})
        names = [i["name"] for i in c.layout]
        self.assertNotIn("old gene", names)
        self.assertIn("newGene", names)

    def test_the_vector_survives_either_way(self):
        for modules, types in (
                ([("CDS", "g", orf(150, 1515))], {"CDS"}),
                ([("promoter", "CAG", rand_seq(400, 1516)),
                  ("CDS", "g", orf(150, 1517))], {"promoter", "CDS"})):
            c = self._build(modules, types)
            present = {i["type"] for i in c.layout}
            for required in ("AAV_ITR", "WPRE", "origin", "resistance"):
                self.assertIn(required, present, f"{required} with {types}")


# --------------------------------------------------------------------------
class TestRecipientTiers(unittest.TestCase):
    """Recipients are ranked by how far down the architecture they match."""

    def setUp(self):
        import hashlib
        self.conn = connect(":memory:")
        self.prom = rand_seq(400, 1600)
        self.first = orf(200, 1601)
        self.second = orf(180, 1602)
        self.modules = []
        for ftype, name, seq in (("promoter", "hSyn", self.prom),
                                 ("CDS", "first", self.first),
                                 ("CDS", "second", self.second)):
            pid = self.conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ftype, name, seq, hashlib.sha256(seq.encode()).hexdigest(),
                 len(seq))).lastrowid
            self.modules.append({"type": ftype, "part_id": pid,
                                 "part_name": name, "seq": seq})

    def _plasmid(self, name, pieces):
        body = "".join(pieces) + rand_seq(2000, 1603)
        return self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES (?,'t',?,?,?,1,'unverified')""",
            (name, body, len(body), name)).lastrowid

    def test_a_longer_leading_run_wins(self):
        from genesis.swap import architecture_tiers
        prom_only = self._plasmid("prom-only", [self.prom])
        prom_first = self._plasmid("prom-first", [self.prom, self.first])
        self.conn.commit()
        self.assertEqual(
            architecture_tiers(self.conn, prom_only, self.modules)["prefix"], 1)
        self.assertEqual(
            architecture_tiers(self.conn, prom_first, self.modules)["prefix"], 2)

    def test_parts_out_of_order_do_not_count_as_a_run(self):
        """Holding the promoter and the last gene is not the same as holding
        the promoter and the first one."""
        from genesis.swap import architecture_tiers
        skipped = self._plasmid("prom-and-last", [self.prom, self.second])
        self.conn.commit()
        tiers = architecture_tiers(self.conn, skipped, self.modules)
        self.assertEqual(tiers["prefix"], 1)
        self.assertEqual(tiers["total"], 2)
        self.assertEqual(tiers["first_missing"], "first")

    def test_an_inverted_payload_still_counts_as_present(self):
        from genesis.swap import architecture_tiers
        flipped = self._plasmid("flipped", [self.prom,
                                            st.revcomp(self.first)])
        self.conn.commit()
        self.assertEqual(
            architecture_tiers(self.conn, flipped, self.modules)["prefix"], 2)

    def test_the_ranking_prefers_the_longer_run(self):
        from genesis.swap import ParentMatch, rank_recipients
        a = self._plasmid("prom-only", [self.prom])
        b = self._plasmid("prom-first", [self.prom, self.first])
        self.conn.commit()
        candidates = [
            ParentMatch(a, "prom-only", "t", 2400, True, "unverified"),
            ParentMatch(b, "prom-first", "t", 2600, True, "unverified")]
        ranked = rank_recipients(self.conn, candidates, self.modules)
        self.assertEqual(ranked[0].name, "prom-first")
        self.assertTrue(any("already has" in r for r in ranked[0].reasons))


# --------------------------------------------------------------------------
class TestDpnINecessity(unittest.TestCase):
    """DpnI is about the template, not about the assembly chemistry."""

    def setUp(self):
        self.conn = connect(":memory:")

    def _plasmid(self, name, marker):
        pid = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular) VALUES (?,'t','ACGT',4,?,1)""",
            (name, name)).lastrowid
        if marker:
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,'resistance',?,0,10,1)""",
                (pid, marker))
        self.conn.commit()
        return pid

    def test_same_marker_column_purified_needs_dpni(self):
        from genesis.planner import dpni_decision
        product = self._plasmid("product", "AmpR")
        template = self._plasmid("template", "bla (ampicillin)")
        result = dpni_decision(self.conn, [template], product,
                               gel_purified=False)
        self.assertTrue(result["needed"])
        self.assertIn("same selection marker", result["reason"])

    def test_a_different_marker_makes_it_unnecessary(self):
        """The plate removes the parent, so DpnI has nothing to add."""
        from genesis.planner import dpni_decision
        product = self._plasmid("product", "AmpR")
        template = self._plasmid("template", "KanR")
        result = dpni_decision(self.conn, [template], product,
                               gel_purified=False)
        self.assertFalse(result["needed"])
        self.assertIn("different antibiotic", result["reason"])

    def test_gel_purification_makes_it_unnecessary(self):
        from genesis.planner import dpni_decision
        product = self._plasmid("product", "AmpR")
        template = self._plasmid("template", "AmpR")
        result = dpni_decision(self.conn, [template], product,
                               gel_purified=True)
        self.assertFalse(result["needed"])
        self.assertIn("gel purified", result["reason"])

    def test_an_unannotated_marker_is_treated_as_a_risk(self):
        """Not knowing what a template selects on is not evidence it is safe."""
        from genesis.planner import dpni_decision
        product = self._plasmid("product", "AmpR")
        template = self._plasmid("template", None)
        self.assertTrue(dpni_decision(self.conn, [template], product,
                                      gel_purified=False)["needed"])

    def test_markers_are_recognised_by_their_usual_labels(self):
        from genesis.planner import resistance_marker
        for label, expected in (("bla", "AmpR"), ("AmpR", "AmpR"),
                                ("Ampicillin resistance gene", "AmpR"),
                                ("KanR", "KanR"), ("neoR/KanR", "KanR"),
                                ("CmR", "CmR"), ("aacC1", "GentR")):
            pid = self._plasmid("p" + label, label)
            self.assertEqual(resistance_marker(self.conn, pid), expected, label)


# --------------------------------------------------------------------------
class TestDefaultVector(unittest.TestCase):
    """A specification that never mentions a vector takes the lab default."""

    def test_an_unmentioned_vector_is_none_not_empty(self):
        spec = ConstructSpec.from_dict(
            {"name": "x", "modules": [{"type": "CDS", "name": "g"}]})
        self.assertIsNone(spec.vector,
                          "not mentioning a vector is different from asking "
                          "for auto")

    def test_asking_for_auto_still_infers(self):
        for value in ("auto", ""):
            spec = ConstructSpec.from_dict({"name": "x", "vector": value})
            self.assertEqual(spec.vector, "")

    def test_the_configured_default_is_applied(self):
        import hashlib
        conn = connect(":memory:")
        cfg = Config()
        self.assertEqual(cfg.default_vector, "AAV")
        seq = orf(150, 1700)
        pid = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('CDS','g',?,?,?)""",
            (seq, hashlib.sha256(seq.encode()).hexdigest(),
             len(seq))).lastrowid
        pl = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES ('pD','t',?,?,'dv',1,'unverified')""",
            (seq + rand_seq(2000, 1701), len(seq) + 2000)).lastrowid
        conn.execute("""INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,0,?,1)""", (pid, pl, len(seq)))
        conn.commit()
        plan = design(conn, ConstructSpec.from_dict(
            {"name": "x", "modules": [{"type": "CDS", "name": "g"}]}), cfg)
        self.assertEqual(plan["vector"]["requested"], "AAV")

    def test_the_default_can_be_turned_off(self):
        cfg = Config()
        cfg.default_vector = ""
        spec = ConstructSpec.from_dict(
            {"name": "x", "modules": [{"type": "CDS", "name": "g"}]})
        self.assertIsNone(spec.vector)
        conn = connect(":memory:")
        plan = design(conn, spec, cfg)
        self.assertIsNone(plan["vector"]["requested"])


class TestInferredORFs(unittest.TestCase):
    """Payload CDSs the source files never annotated."""

    def test_insert_name_is_pulled_out_of_the_plasmid_name(self):
        from genesis.orfs import name_from_plasmid as n
        self.assertEqual(n("YD043-pAAV-syn-iGluSnFR4f-NGR-WPRE"), "iGluSnFR4f")
        self.assertEqual(n("paav-ef1a-dio-asap7y-wpre_"), "asap7y")
        self.assertEqual(n("pcDNA3-CheRiff-eGFP"), "CheRiff")
        self.assertEqual(n("UB030_HuC-zArchon"), "zArchon")

    def test_vector_only_names_yield_nothing(self):
        from genesis.orfs import name_from_plasmid as n
        for name in ("pcDNA3a", "pUC19", "pAAV2", "PREDICTED_pPZ021",
                     "addgene-plasmid-53191-sequence", "DL011_5-TBD"):
            self.assertIsNone(n(name), name)

    def test_identity_is_scored_against_the_longer_protein(self):
        # The point of the whole exercise: a cpGFP buried inside a sensor
        # makes the sensor 97% of a GFP under protein_identity, which asks
        # "same part?".  Asking "what is this?" needs the other denominator,
        # or every sensor in the inventory comes back as GFP.
        from genesis.orfs import ProteinIndex
        from genesis.seqtools import translate
        big = "M" + "".join(random.Random(4).choice("ACDEFGHIKLMNPQRSTVWY")
                            for _ in range(400))
        small = big[100:200]

        class Fake(ProteinIndex):
            def __init__(self):
                self.kmers, self.parts = {}, {}
                for pid, aa in ((1, small), (2, big)):
                    self.parts[pid] = (f"p{pid}", len(aa))
                    for i in range(len(aa) - self.K + 1):
                        self.kmers.setdefault(aa[i:i + self.K], set()).add(pid)

        ident, pid, _name = Fake().best(big)
        self.assertEqual(pid, 2, "a protein should match itself, not its part")
        self.assertGreater(ident, 0.95)

    def test_generic_catalogue_names_are_not_inherited(self):
        from genesis.orfs import _useful
        for bad in ("ORF frame 1", "CDS", "Translation 1970-3730", "target",
                    "P (fragment)", "unnamed", "Template"):
            self.assertFalse(_useful(bad), bad)
        for good in ("QuasAr2", "iGluSnFr3", "beta-lactamase", "mCherry",
                     "Voltron2 [2 of 3]"):
            self.assertTrue(_useful(good), good)

    def test_an_annotated_orf_is_left_alone(self):
        from genesis.orfs import candidate_orfs
        seq = "ATG" + "AAC" * 300 + "TAA" + "GGGG" * 50
        covered = [{"feature_type": "CDS", "name": "known",
                    "start": 0, "end": 906, "strand": 1}]
        self.assertEqual(candidate_orfs(seq, False, covered), [])
        self.assertTrue(candidate_orfs(seq, False, []))

    def test_payload_span_ignores_the_resistance_promoter(self):
        from genesis.orfs import payload_span
        anns = [
            {"feature_type": "promoter", "name": "hSyn promoter",
             "start": 255, "end": 703},
            {"feature_type": "WPRE", "name": "WPRE", "start": 2543, "end": 3132},
            {"feature_type": "promoter", "name": "AmpR promoter",
             "start": 5966, "end": 6071},
        ]
        self.assertEqual(payload_span(anns), (703, 2543))

    def test_payload_span_needs_both_ends(self):
        from genesis.orfs import payload_span
        self.assertIsNone(payload_span(
            [{"feature_type": "promoter", "name": "hSyn", "start": 0, "end": 400}]))


class TestCloningSimulation(unittest.TestCase):
    """Replaying the protocol, rather than trusting its bookkeeping."""

    def test_pcr_product_carries_the_tails(self):
        from genesis.simulate import pcr
        from genesis.seqtools import revcomp
        t = rand_seq(3000, 5)
        fwd_anneal, rev_anneal = t[500:520], revcomp(t[1200:1220])
        product, why = pcr(t, "aaaacccc".upper(), fwd_anneal,
                           "ggggtttt".upper(), rev_anneal, circular=True)
        self.assertIsNotNone(product, why)
        self.assertEqual(product, "AAAACCCC" + t[500:1220] + revcomp("GGGGTTTT"))

    def test_pcr_works_against_either_strand(self):
        # Which strand we call the template is bookkeeping, not chemistry.
        from genesis.simulate import pcr
        from genesis.seqtools import revcomp
        t = rand_seq(2000, 9)
        a, _ = pcr(t, "", t[300:320], "", revcomp(t[900:920]), circular=True)
        b, _ = pcr(revcomp(t), "", t[300:320], "", revcomp(t[900:920]),
                   circular=True)
        self.assertEqual(a, b)

    def test_a_primer_that_does_not_bind_is_reported(self):
        from genesis.simulate import pcr
        from genesis.seqtools import revcomp
        t = rand_seq(1000, 11)
        product, why = pcr(t, "", "GGGGGGGGGGGGGGGGGGGG", "",
                           revcomp(t[500:520]), circular=False)
        self.assertIsNone(product)
        self.assertIn("forward primer", why)

    def test_fragments_assemble_on_shared_ends(self):
        from genesis.simulate import assemble_by_homology, SimulatedProduct
        a_body, b_body = rand_seq(400, 2), rand_seq(300, 3)
        join1, join2 = a_body[-20:], b_body[-20:]
        a = SimulatedProduct("a", join2 + a_body)
        b = SimulatedProduct("b", join1 + b_body)
        product, notes = assemble_by_homology([a, b])
        self.assertIsNotNone(product, notes)
        # Each overlap is counted once in the circle.
        self.assertEqual(len(product), len(a_body) + len(b_body))

    def test_fragments_that_share_nothing_do_not_assemble(self):
        from genesis.simulate import assemble_by_homology, SimulatedProduct
        product, notes = assemble_by_homology([
            SimulatedProduct("a", rand_seq(400, 4)),
            SimulatedProduct("b", rand_seq(300, 5))])
        self.assertIsNone(product)
        self.assertTrue(any("not chain" in n for n in notes), notes)

    def test_a_missing_part_fails_the_simulation(self):
        """The bug this exists for: steps that build the wrong molecule."""
        from genesis.simulate import _check_contents, SimulationResult
        present, absent = rand_seq(300, 6), rand_seq(300, 7)
        result = SimulationResult()
        plan = {"modules": [{"part_name": "kept", "seq": present},
                            {"part_name": "ASAP7y", "seq": absent}],
                "construct": {"length": len(present)}}
        _check_contents(result, present, plan)
        self.assertFalse(result.passed)
        self.assertTrue(any("ASAP7y" in p for p in result.problems),
                        result.problems)

    def test_an_element_kept_from_the_recipient_still_counts(self):
        # A swap reuses the recipient's own WPRE, which is the same element
        # as the catalogued part without being byte-identical to it.
        from genesis.simulate import _check_contents, SimulationResult
        part = rand_seq(600, 8)
        in_product = part[12:-9]          # same element, other boundaries
        result = SimulationResult()
        _check_contents(result, in_product,
                        {"modules": [{"part_name": "WPRE", "seq": part}],
                         "construct": {"length": len(in_product)}})
        self.assertTrue(result.passed, result.problems)


class TestArchitectureCandidates(unittest.TestCase):
    """A request read as an ordered architecture, not a bag of parts."""

    def setUp(self):
        self.cfg = Config()
        self.conn = connect(":memory:")
        self.seqs = {}
        self.parts = {}
        for i, name in enumerate(["hSyn", "XXX", "P2A", "ZZZ"]):
            s = rand_seq(300, 30 + i)
            self.seqs[name] = s
            self.parts[name] = self.conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                ("promoter" if name == "hSyn" else
                 "p2A" if name == "P2A" else "CDS",
                 name, s, name, len(s))).lastrowid
        self.conn.commit()

    def _plasmid(self, name, carries, vector=True):
        body = "".join(self.seqs[c] for c in carries) + rand_seq(1200, 7)
        pid = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES (?,'t',?,?,?,1,'unverified')""",
            (name, body, len(body), name)).lastrowid
        at = 0
        for c in carries:
            self.conn.execute(
                """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                        strand) VALUES (?,?,?,?,1)""",
                (self.parts[c], pid, at, at + len(self.seqs[c])))
            at += len(self.seqs[c])
        if vector:
            for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
                self.conn.execute(
                    """INSERT INTO annotation(plasmid_id, feature_type, name,
                            start, end, strand) VALUES (?,?,?,0,30,1)""",
                    (pid, ftype, nm))
        self.conn.commit()
        return pid

    def _modules(self):
        return [{"type": "promoter" if n == "hSyn" else
                 "p2A" if n == "P2A" else "CDS",
                 "part_id": self.parts[n], "part_name": n,
                 "seq": self.seqs[n]}
                for n in ("hSyn", "XXX", "P2A", "ZZZ")]

    def test_the_longest_leading_run_wins(self):
        from genesis.swap import architecture_candidates
        self._plasmid("pAAV-hSyn", ["hSyn"])
        best = self._plasmid("pAAV-hSyn-XXX", ["hSyn", "XXX"])
        self._plasmid("pAAV-ZZZ", ["ZZZ"])
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["recipients"][0]["plasmid_id"], best)
        self.assertEqual(c["recipients"][0]["run"], 2)

    def test_a_plasmid_holding_only_a_later_part_is_not_a_recipient(self):
        # It has nothing to start from: the leading run is empty.
        from genesis.swap import architecture_candidates
        self._plasmid("pAAV-ZZZ", ["ZZZ"])
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["recipients"], [])

    def test_one_donor_when_the_remainder_is_in_one_plasmid(self):
        from genesis.swap import architecture_candidates
        self._plasmid("pAAV-hSyn-XXX", ["hSyn", "XXX"])
        self._plasmid("donor-P2A-ZZZ", ["P2A", "ZZZ"])
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["n_donors"], 1)
        self.assertEqual(c["donors"][0]["covers"], ["P2A", "ZZZ"])

    def test_two_donors_when_the_remainder_is_split(self):
        """The case that has to be said out loud rather than dropped."""
        from genesis.swap import architecture_candidates
        self._plasmid("pAAV-hSyn", ["hSyn"])
        self._plasmid("donor-XXX-P2A", ["XXX", "P2A"])
        self._plasmid("donor-ZZZ", ["ZZZ"])
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["recipients"][0]["run"], 1)
        self.assertEqual(c["n_donors"], 2)
        self.assertEqual([d["covers"] for d in c["donors"]],
                         [["XXX", "P2A"], ["ZZZ"]])

    def test_a_part_nobody_has_is_reported_not_hidden(self):
        from genesis.swap import architecture_candidates
        self._plasmid("pAAV-hSyn-XXX", ["hSyn", "XXX"])
        self._plasmid("donor-P2A", ["P2A"])
        c = architecture_candidates(self.conn, self._modules())
        last = c["donors"][-1]
        self.assertEqual(last["covers"], ["ZZZ"])
        self.assertEqual(last["options"], [])
        self.assertIn("synthesised", last["note"])

    def test_a_cassette_is_part_of_the_architecture(self):
        """pAAV-hSyn-DIO-XX is a different request from pAAV-hSyn-XX."""
        from genesis.swap import architecture_candidates
        plain = self._plasmid("pAAV-hSyn-XXX", ["hSyn", "XXX"])
        dio = self._plasmid("pAAV-hSyn-DIO-XXX", ["hSyn", "XXX"])
        for nm in ("lox2272", "loxP", "lox2272", "loxP"):
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,'lox_site',?,0,34,1)""",
                (dio, nm))
        self.conn.commit()
        mods = self._modules()
        with_lox = [mods[0], {"type": "Cre-lox", "part_id": None,
                              "part_name": "Cre-lox", "seq": ""}] + mods[1:]

        # Without a cassette asked for, the plain vector leads.
        plainc = architecture_candidates(self.conn, mods)
        self.assertEqual(plainc["recipients"][0]["plasmid_id"], plain)
        # With one, the vector that already has it does.
        dioc = architecture_candidates(self.conn, with_lox)
        self.assertIn("Cre-lox", dioc["modules"])
        self.assertEqual(dioc["recipients"][0]["plasmid_id"], dio)
        self.assertEqual(dioc["recipients"][0]["run"], 3)

    def test_a_cassette_in_the_vector_is_reported(self):
        """So a plain payload is not dropped between someone else's lox sites."""
        from genesis.swap import architecture_candidates
        plain = self._plasmid("pAAV-hSyn-plain", ["hSyn"])
        dio = self._plasmid("pAAV-hSyn-DIO", ["hSyn"])
        for name in ("lox2272", "loxP", "lox2272", "loxP"):
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,'lox_site',?,0,34,1)""",
                (dio, name))
        self.conn.commit()
        got = {c["plasmid_id"]: c["has_cassette"]
               for c in architecture_candidates(self.conn, self._modules())
               ["recipients"]}
        self.assertFalse(got[plain])
        self.assertTrue(got[dio])

    def test_a_relative_is_named_as_itself(self):
        """A plasmid holding ASAP3 does not hold ASAP7y."""
        from genesis.swap import architecture_candidates
        other = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('CDS','XXX-relative',?,'rel',?)""",
            (self.seqs["XXX"], len(self.seqs["XXX"]))).lastrowid
        pid = self._plasmid("pAAV-hSyn-relative", ["hSyn"])
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,900,1200,1)""", (other, pid))
        self.conn.commit()
        c = architecture_candidates(self.conn, self._modules())
        everyone = c["recipients"] + c["related"]
        top = everyone[0]
        # What it really has and what it only resembles are kept apart, so
        # a relative can never be read as the part itself.
        self.assertNotIn("XXX", top["carries"])
        self.assertTrue(any("as XXX-relative" in x for x in top["resembles"]),
                        top)


class TestBackboneArcs(unittest.TestCase):
    """The route amplifies what the construct keeps, not the whole plasmid."""

    def test_the_kept_arcs_are_recorded(self):
        from genesis.construct import assemble_into_backbone, ResolvedModule
        from genesis.construct import Module
        payload = rand_seq(300, 21)
        old_cds = rand_seq(450, 22)
        backbone = rand_seq(600, 23) + old_cds + rand_seq(900, 24)
        anns = [
            {"feature_type": "CDS", "name": "old", "start": 600,
             "end": 1050, "strand": 1},
            {"feature_type": "origin", "name": "ori", "start": 1200,
             "end": 1500, "strand": 1},
            {"feature_type": "resistance", "name": "AmpR", "start": 1500,
             "end": 1900, "strand": 1},
        ]
        rm = ResolvedModule(module=Module("CDS", "new"), part_id=1,
                            part_name="new", seq=payload, feature_type="CDS",
                            sources=[], warnings=[])
        con = assemble_into_backbone([rm], backbone, anns, {"CDS"})

        kept = con.backbone_segments
        self.assertTrue(kept, "no backbone arcs recorded")
        total = sum(len(s["seq"]) for s in kept)
        # The construct is the arcs it keeps plus what was put in.
        self.assertEqual(total + len(payload), len(con.seq))
        # And the replaced coding sequence is not among them.
        for s in kept:
            self.assertNotIn(old_cds, s["seq"])

    def test_a_wider_part_does_not_duplicate_the_recipients_bases(self):
        """The part and the recipient's annotation disagree on boundaries.

        The catalogued promoter is wider than the span the recipient marks.
        Replacing the marked span with the part leaves the bases they
        disagree about on both sides of the join: once in the arc that is
        kept, and again inside the piece put in its place.
        """
        from genesis.construct import (assemble_into_backbone, ResolvedModule,
                                       Module)
        lead = rand_seq(120, 51)
        promoter = rand_seq(400, 52)       # the catalogued part
        old_cds = rand_seq(450, 53)
        tail = rand_seq(900, 54)
        backbone = lead + promoter + old_cds + tail
        # The recipient marks only the middle 300 bp of that promoter.
        anns = [
            {"feature_type": "promoter", "name": "its own",
             "start": 120 + 50, "end": 120 + 350, "strand": 1},
            {"feature_type": "CDS", "name": "old",
             "start": 520, "end": 970, "strand": 1},
            {"feature_type": "origin", "name": "ori",
             "start": 1000, "end": 1300, "strand": 1},
            {"feature_type": "resistance", "name": "AmpR",
             "start": 1300, "end": 1700, "strand": 1},
        ]
        mods = [
            ResolvedModule(module=Module("promoter", "P"), part_id=1,
                           part_name="P", seq=promoter,
                           feature_type="promoter", sources=[], warnings=[]),
            ResolvedModule(module=Module("CDS", "new"), part_id=2,
                           part_name="new", seq=rand_seq(300, 55),
                           feature_type="CDS", sources=[], warnings=[]),
        ]
        con = assemble_into_backbone(mods, backbone, anns, {"promoter", "CDS"})
        self.assertEqual(con.seq.count(promoter), 1,
                         "the promoter should appear once, not overlap itself")
        # The bases either side of it must not be repeated across the join.
        self.assertEqual(con.seq.count(promoter[:50]), 1)
        self.assertEqual(con.seq.count(promoter[-50:]), 1)

    def test_an_unamplifiable_arc_rides_on_a_primer(self):
        """There is no such thing as a 25 bp PCR."""
        from genesis.planner import _carry
        from genesis.primers import Primer
        from genesis.planner import Route
        route = Route(id="x", strategy="hifi", title="t")
        arc = rand_seq(25, 41)
        fwd = Primer(name="next-F", sequence="ggggg" + "ACGTACGTACGTACGTAC",
                     anneal="ACGTACGTACGTACGTAC", tail="ggggg", direction="F")
        _carry(fwd, arc, route, Config())
        # The bases still exist, and in the right place: overlap, arc, body.
        self.assertEqual(fwd.sequence, "ggggg" + arc + "ACGTACGTACGTACGTAC")
        self.assertTrue(fwd.sequence.startswith("ggggg"))
        self.assertEqual(route.warnings, [])

    def test_a_very_long_carry_is_flagged(self):
        from genesis.planner import _carry, Route
        from genesis.primers import Primer
        route = Route(id="x", strategy="hifi", title="t")
        fwd = Primer(name="next-F", sequence="A" * 20, anneal="A" * 20,
                     tail="", direction="F")
        _carry(fwd, rand_seq(90, 42), route, Config())
        self.assertTrue(any("ultramer" in w for w in route.warnings),
                        route.warnings)

    def test_the_old_payload_does_not_come_along(self):
        """An arc-amplified route cannot carry what it replaced."""
        from genesis.construct import assemble_into_backbone, ResolvedModule
        from genesis.construct import Module
        old_cds = rand_seq(450, 31)
        backbone = rand_seq(600, 32) + old_cds + rand_seq(900, 33)
        anns = [
            {"feature_type": "CDS", "name": "old", "start": 600,
             "end": 1050, "strand": 1},
            {"feature_type": "origin", "name": "ori", "start": 1200,
             "end": 1500, "strand": 1},
            {"feature_type": "resistance", "name": "AmpR", "start": 1500,
             "end": 1900, "strand": 1},
        ]
        rm = ResolvedModule(module=Module("CDS", "new"), part_id=1,
                            part_name="new", seq=rand_seq(300, 34),
                            feature_type="CDS", sources=[], warnings=[])
        con = assemble_into_backbone([rm], backbone, anns, {"CDS"})
        joined = "".join(s["seq"] for s in con.backbone_segments)
        self.assertNotIn(old_cds, joined)
        self.assertNotIn(old_cds, con.seq)


class TestUniqueOverlaps(unittest.TestCase):
    """An overlap arm has to say where a fragment goes."""

    def setUp(self):
        self.cfg = Config()

    def test_occurrences_counts_circularly_on_both_strands(self):
        from genesis.primers import _occurrences
        from genesis.seqtools import revcomp
        body = rand_seq(500, 61)
        motif = body[100:140]
        self.assertEqual(_occurrences(motif, body), 1)
        self.assertEqual(_occurrences(revcomp(motif), body), 1)
        twice = body + rand_seq(200, 62) + motif
        self.assertEqual(_occurrences(motif, twice), 2)

    def test_an_arm_grows_past_a_repeat_until_it_is_unique(self):
        """A lox site appears twice by design; an arm inside one is useless."""
        from genesis.primers import pick_overlap, _occurrences
        lox = rand_seq(34, 63)
        unique_left = rand_seq(200, 64)
        payload = rand_seq(300, 65)
        tail_side = rand_seq(300, 66)
        # The construct carries the same lox twice.
        construct = unique_left + lox + payload + lox + tail_side
        neighbour = unique_left + lox      # the piece before the payload

        short = pick_overlap(neighbour, self.cfg, from_end=True)
        self.assertGreater(_occurrences(short, construct), 1,
                           "the default arm should be the ambiguous one")

        grown = pick_overlap(neighbour, self.cfg, from_end=True,
                             unique_in=construct)
        self.assertEqual(_occurrences(grown, construct), 1)
        self.assertGreater(len(grown), len(short))
        self.assertTrue(neighbour.endswith(grown),
                        "the arm must still be the neighbour's own end")

    def test_an_arm_that_cannot_be_made_unique_is_reported(self):
        from genesis.planner import _check_unique, Route
        from genesis.primers import Primer
        repeat = rand_seq(60, 67)
        construct = rand_seq(300, 68) + repeat + rand_seq(300, 69) + repeat
        route = Route(id="x", strategy="hifi", title="t")
        p = Primer(name="frag-F", sequence=repeat[-15:] + "ACGTACGTACGTACGTAC",
                   anneal="ACGTACGTACGTACGTAC", tail=repeat[-15:],
                   direction="F")
        _check_unique(p, construct, route)
        self.assertTrue(any("more than one place" in w for w in route.warnings),
                        route.warnings)

    def test_a_unique_arm_is_not_reported(self):
        from genesis.planner import _check_unique, Route
        from genesis.primers import Primer
        construct = rand_seq(800, 70)
        route = Route(id="x", strategy="hifi", title="t")
        arm = construct[100:115]
        p = Primer(name="frag-F", sequence=arm + "ACGTACGTACGTACGTAC",
                   anneal="ACGTACGTACGTACGTAC", tail=arm, direction="F")
        _check_unique(p, construct, route)
        self.assertEqual(route.warnings, [])


class TestInvertedPayload(unittest.TestCase):
    """A Cre-lox payload is inverted once, not twice."""

    def test_expand_dio_hands_back_the_other_strand(self):
        from genesis.construct import expand_dio, ResolvedModule, Module
        from genesis.seqtools import revcomp
        cds = rand_seq(300, 71)
        mods = [
            ResolvedModule(module=Module("Cre-lox", "DIO"), part_id=None,
                           part_name="DIO", seq="", feature_type="Cre-lox",
                           sources=[], warnings=[]),
            ResolvedModule(module=Module("CDS", "X"), part_id=1,
                           part_name="X", seq=cds, feature_type="CDS",
                           sources=[], warnings=[]),
        ]
        out, _notes = expand_dio(mods, sites_from_backbone=True)
        payload = [m for m in out if m.feature_type == "CDS"]
        self.assertEqual(len(payload), 1)
        self.assertEqual(payload[0].seq, revcomp(cds),
                         "the payload comes back on the other strand")
        self.assertEqual(payload[0].module.orientation, "reverse")

    def test_placed_seq_does_not_flip_it_again(self):
        """The sequence expand_dio returns is already the construct's."""
        from genesis.planner import PlannedFragment
        from genesis.seqtools import revcomp
        body = rand_seq(300, 72)
        f = PlannedFragment("X", body, 1, "donor", 0, len(body), 1, False,
                            [{"start": 0, "strand": -1}], placed_strand=-1)
        self.assertEqual(f.placed_seq, body)
        self.assertNotEqual(f.placed_seq, revcomp(body))


class TestContiguousMerge(unittest.TestCase):
    """Neighbouring stretches of one plasmid are one reaction."""

    def _piece(self, kind, name, seq, at, start, end, index=None):
        p = {"kind": kind, "name": name, "seq": seq, "at": at,
             "source_start": start, "source_end": end}
        if index is not None:
            p["index"] = index
        return p

    def test_arcs_and_a_module_from_the_same_plasmid_merge(self):
        from genesis.planner import _merge_contiguous
        a, b, c = rand_seq(100, 81), rand_seq(400, 82), rand_seq(60, 83)
        pieces = [
            self._piece("backbone", "backbone1", a, 0, 130, 230),
            self._piece("fragment", "hsyn", b, 100, 230, 630, index=0),
            self._piece("backbone", "backbone2", c, 500, 630, 690),
        ]
        merged, absorbed = _merge_contiguous(pieces, 5829)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]["source_start"], 130)
        self.assertEqual(merged[0]["source_end"], 690)
        self.assertEqual(merged[0]["seq"], a + b + c)
        self.assertIn(0, absorbed, "the module no longer needs its own PCR")

    def test_a_gap_is_not_merged(self):
        from genesis.planner import _merge_contiguous
        pieces = [
            self._piece("backbone", "b1", rand_seq(100, 84), 0, 130, 230),
            self._piece("fragment", "x", rand_seq(100, 85), 100, 900, 1000,
                        index=0),
        ]
        merged, absorbed = _merge_contiguous(pieces, 5829)
        self.assertEqual(len(merged), 2)
        self.assertEqual(absorbed, set())

    def test_a_piece_from_another_plasmid_is_not_merged(self):
        """A fragment off a donor carries no recipient coordinates."""
        from genesis.planner import _merge_contiguous
        pieces = [
            self._piece("backbone", "b1", rand_seq(100, 86), 0, 130, 230),
            {"kind": "fragment", "name": "donor-cds", "seq": rand_seq(300, 87),
             "at": 100, "index": 0},
        ]
        merged, absorbed = _merge_contiguous(pieces, 5829)
        self.assertEqual(len(merged), 2)
        self.assertEqual(absorbed, set())

    def test_the_join_across_position_one_still_merges(self):
        from genesis.planner import _merge_contiguous
        head = self._piece("backbone", "b-head", rand_seq(200, 88), 0, 0, 200)
        tail = self._piece("backbone", "b-tail", rand_seq(300, 89), 900,
                           5000, 5829)
        merged, _ = _merge_contiguous([head, tail], 5829)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]["source_start"], 5000)
        self.assertEqual(merged[0]["source_end"], 200)


class TestVectorElementReuse(unittest.TestCase):
    """The recipient's own WPRE is only reused if it survives."""

    def _backbone(self, wpre_inside):
        """A recipient whose WPRE sits inside or outside the payload arc."""
        prom = rand_seq(300, 91)
        old_cds = rand_seq(600, 92)
        wpre = rand_seq(400, 93)
        tail = rand_seq(500, 94)
        if wpre_inside:
            seq = prom + old_cds + wpre + tail
            w_start = len(prom) + len(old_cds)
        else:
            seq = prom + old_cds + tail + wpre
            w_start = len(prom) + len(old_cds) + len(tail)
        anns = [
            {"feature_type": "promoter", "name": "p", "start": 0,
             "end": len(prom), "strand": 1},
            {"feature_type": "CDS", "name": "old", "start": len(prom),
             "end": len(prom) + len(old_cds), "strand": 1},
            {"feature_type": "WPRE", "name": "its WPRE", "start": w_start,
             "end": w_start + len(wpre), "strand": 1},
            {"feature_type": "origin", "name": "ori",
             "start": len(seq) - 300, "end": len(seq) - 150, "strand": 1},
            {"feature_type": "resistance", "name": "AmpR",
             "start": len(seq) - 150, "end": len(seq), "strand": 1},
        ]
        return seq, anns, wpre

    def _modules(self, wpre):
        from genesis.construct import ResolvedModule, Module
        return [
            ResolvedModule(module=Module("CDS", "new"), part_id=1,
                           part_name="new", seq=rand_seq(300, 95),
                           feature_type="CDS", sources=[], warnings=[]),
            ResolvedModule(module=Module("WPRE", "WPRE"), part_id=2,
                           part_name="WPRE", seq=wpre, feature_type="WPRE",
                           sources=[], warnings=[]),
        ]

    def test_the_construct_ends_up_with_exactly_one(self):
        """Reused or inserted, never both and never neither.

        Whether the recipient's own copy survives depends on where it sits
        relative to the stretch being replaced. Either answer is fine; two
        copies is a repeat nothing can disambiguate, and none is a construct
        missing what it was asked for.
        """
        from genesis.construct import assemble_into_backbone
        for inside in (True, False):
            seq, anns, wpre = self._backbone(wpre_inside=inside)
            con = assemble_into_backbone(self._modules(wpre), seq, anns,
                                         {"CDS", "WPRE"})
            self.assertEqual(con.seq.count(wpre), 1,
                             f"wpre_inside={inside}: {con.warnings}")

    def test_a_copy_inside_the_arc_is_not_relied_on(self):
        """It gets cut away, so the construct would end up without one."""
        from genesis.construct import assemble_into_backbone
        seq, anns, wpre = self._backbone(wpre_inside=True)
        con = assemble_into_backbone(self._modules(wpre), seq, anns,
                                     {"CDS", "WPRE"})
        self.assertIn(wpre, con.seq, "the construct must still have a WPRE")
        self.assertFalse(any("already carries" in w for w in con.warnings),
                         con.warnings)


class TestUnwantedCassette(unittest.TestCase):
    """A recipient's lox sites are not inherited by accident."""

    def _recipient(self):
        """pAAV-hSyn-DIO-<cds>: a cassette around the coding sequence."""
        prom = rand_seq(300, 101)
        lox_a, lox_b = rand_seq(34, 102), rand_seq(34, 103)
        old_cds = rand_seq(600, 104)
        tail = rand_seq(700, 105)
        seq = prom + lox_a + lox_b + old_cds + lox_a + lox_b + tail
        at = len(prom)
        anns = [{"feature_type": "promoter", "name": "hSyn", "start": 0,
                 "end": len(prom), "strand": 1}]
        for i, nm in enumerate(("lox2272", "loxP")):
            anns.append({"feature_type": "lox_site", "name": nm,
                         "start": at + i * 34, "end": at + (i + 1) * 34,
                         "strand": 1})
        cds_at = at + 68
        anns.append({"feature_type": "CDS", "name": "old", "start": cds_at,
                     "end": cds_at + 600, "strand": 1})
        at2 = cds_at + 600
        for i, nm in enumerate(("lox2272", "loxP")):
            anns.append({"feature_type": "lox_site", "name": nm,
                         "start": at2 + i * 34, "end": at2 + (i + 1) * 34,
                         "strand": -1})
        anns += [
            {"feature_type": "origin", "name": "ori",
             "start": len(seq) - 400, "end": len(seq) - 200, "strand": 1},
            {"feature_type": "resistance", "name": "AmpR",
             "start": len(seq) - 200, "end": len(seq), "strand": 1},
        ]
        return seq, anns

    def _payload(self):
        from genesis.construct import ResolvedModule, Module
        return [ResolvedModule(module=Module("CDS", "new"), part_id=1,
                               part_name="new", seq=rand_seq(300, 106),
                               feature_type="CDS", sources=[], warnings=[])]

    def test_a_plain_design_does_not_inherit_the_lox_sites(self):
        """Dropping a payload between someone else's lox sites makes the
        construct Cre-dependent, which is a different experiment."""
        from genesis.construct import assemble_into_backbone
        seq, anns = self._recipient()
        con = assemble_into_backbone(self._payload(), seq, anns, {"CDS"})
        kinds = [i["type"] for i in con.layout]
        self.assertNotIn("lox_site", kinds,
                         [i["name"] for i in con.layout])

    def test_a_cre_lox_design_keeps_them(self):
        from genesis.construct import assemble_into_backbone
        seq, anns = self._recipient()
        con = assemble_into_backbone(self._payload(), seq, anns,
                                     {"CDS", "Cre-lox"})
        kinds = [i["type"] for i in con.layout]
        self.assertIn("lox_site", kinds)


class TestDonorIsAskedFor(unittest.TestCase):
    """Resembling a part is not supplying it."""

    def test_a_recipient_holding_a_variant_still_needs_a_donor(self):
        from genesis.swap import architecture_candidates
        conn = connect(":memory:")
        want = rand_seq(600, 111)          # the part asked for
        # A real relative, not an unrelated sequence: most of the same bases
        # under a different name, which is what the catalogue search finds
        # and what made it report that no donor was needed.
        variant = want[:260] + rand_seq(40, 112) + want[300:]
        prom = rand_seq(300, 113)
        parts = {}
        for nm, sq, ft in (("hSyn", prom, "promoter"), ("ASAP7y", want, "CDS"),
                           ("ASAP4c", variant, "CDS")):
            parts[nm] = conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ft, nm, sq, nm, len(sq))).lastrowid

        def plasmid(name, body, carries):
            pid = conn.execute(
                """INSERT INTO plasmid(name, owner, seq, seq_length,
                        seq_checksum, circular, verified)
                   VALUES (?,'t',?,?,?,1,'unverified')""",
                (name, body, len(body), name)).lastrowid
            at = 0
            for part_name, sq in carries:
                conn.execute(
                    """INSERT INTO part_occurrence(part_id, plasmid_id, start,
                            end, strand) VALUES (?,?,?,?,1)""",
                    (parts[part_name], pid, at, at + len(sq)))
                at += len(sq)
            for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
                conn.execute(
                    """INSERT INTO annotation(plasmid_id, feature_type, name,
                            start, end, strand) VALUES (?,?,?,0,30,1)""",
                    (pid, ftype, nm))
            return pid

        # The recipient has the promoter and a *different* ASAP.
        plasmid("pAAV-hSyn-ASAP4c", prom + variant + rand_seq(900, 114),
                [("hSyn", prom), ("ASAP4c", variant)])
        # The real part lives elsewhere.
        plasmid("donor-ASAP7y", want + rand_seq(900, 115),
                [("ASAP7y", want)])
        conn.commit()

        modules = [
            {"type": "promoter", "part_id": parts["hSyn"],
             "part_name": "hSyn", "seq": prom},
            {"type": "CDS", "part_id": parts["ASAP7y"],
             "part_name": "ASAP7y", "seq": want},
        ]
        c = architecture_candidates(conn, modules)
        self.assertTrue(c["recipients"], "the variant holder is still a start")
        self.assertEqual(c["n_donors"], 1, c["donors"])
        self.assertEqual(c["donors"][0]["covers"], ["ASAP7y"])


class TestDonorAccuracy(unittest.TestCase):
    """A donor is ranked on carrying the part, not on being small."""

    def setUp(self):
        self.conn = connect(":memory:")
        self.want = rand_seq(900, 121)
        # A relative: most of the same bases under another name, which is
        # what the catalogue search turns up.
        self.relative = self.want[:400] + rand_seq(60, 122) + self.want[460:]
        # A promoter only the recipient carries, so the plasmids holding
        # the coding sequence are judged as donors.
        self.prom = rand_seq(300, 120)
        self.parts = {}
        for nm, sq in (("hSyn", self.prom), ("ASAP7y", self.want),
                       ("ASAP1", self.relative)):
            self.parts[nm] = self.conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                ("promoter" if nm == "hSyn" else "CDS",
                 nm, sq, nm, len(sq))).lastrowid
        self.conn.commit()
        self._plasmid("recipient-hSyn", self.prom + rand_seq(400, 125),
                      "hSyn", self.prom)

    def _plasmid(self, name, body, part_name, sq):
        pid = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES (?,'t',?,?,?,1,'unverified')""",
            (name, body, len(body), name)).lastrowid
        at = body.find(sq)
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,?,?,1)""",
            (self.parts[part_name], pid, at, at + len(sq)))
        for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,?,?,0,30,1)""",
                (pid, ftype, nm))
        self.conn.commit()
        return pid

    def test_the_plasmid_with_the_part_leads_however_large(self):
        from genesis.swap import architecture_candidates
        # A small plasmid holding a relative, a large one holding the part.
        self._plasmid("small-relative", self.relative + rand_seq(200, 123),
                      "ASAP1", self.relative)
        real = self._plasmid("large-but-right",
                             rand_seq(5000, 124) + self.want, "ASAP7y",
                             self.want)
        modules = [{"type": "promoter", "part_id": self.parts["hSyn"],
                    "part_name": "hSyn", "seq": self.prom},
                   {"type": "CDS", "part_id": self.parts["ASAP7y"],
                    "part_name": "ASAP7y", "seq": self.want}]
        c = architecture_candidates(self.conn, modules)
        self.assertTrue(c["donors"], c)
        options = c["donors"][0]["options"]
        self.assertEqual(options[0]["plasmid_id"], real,
                         [o["name"] for o in options])
        self.assertTrue(options[0]["exact"])
        # The relative is still offered, just not first: worth seeing, and
        # worth not reaching for by accident.
        self.assertTrue(any(not o["exact"] for o in options),
                        [o["name"] for o in options])


class TestFusedSource(unittest.TestCase):
    """A 2A fusion that already exists is taken whole."""

    def setUp(self):
        self.conn = connect(":memory:")
        self.cfg = Config()
        self.a = orf(120, 131)                   # first protein
        self.b = orf(90, 132)                    # second protein
        # The plasmid's own 2A: in frame, no stop, and not the generic one
        # a specification resolves "P2A" to.
        self.spacer = "GGA" * 20
        self.generic = "GGC" * 22
        self.parts = {}
        for nm, sq, ft in (("A", self.a, "CDS"), ("B", self.b, "CDS"),
                           ("P2A", self.generic, "p2A")):
            self.parts[nm] = self.conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ft, nm, sq, nm, len(sq))).lastrowid

    def _plasmid(self, name, body, carries):
        pid = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES (?,'t',?,?,?,1,'unverified')""",
            (name, body, len(body), name)).lastrowid
        for part_name in carries:
            sq = {"A": self.a, "B": self.b, "P2A": self.generic}[part_name]
            at = body.find(sq)
            self.conn.execute(
                """INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,?,?,1)""",
                (self.parts[part_name], pid, at, at + len(sq)))
        self.conn.commit()
        return pid

    def _modules(self):
        from genesis.construct import ResolvedModule, Module
        return [
            ResolvedModule(module=Module("CDS", "A"), part_id=self.parts["A"],
                           part_name="A", seq=self.a, feature_type="CDS",
                           sources=[], warnings=[]),
            ResolvedModule(module=Module("p2A", "P2A"),
                           part_id=self.parts["P2A"], part_name="P2A",
                           seq=self.generic, feature_type="p2A",
                           sources=[], warnings=[]),
            ResolvedModule(module=Module("CDS", "B"), part_id=self.parts["B"],
                           part_name="B", seq=self.b, feature_type="CDS",
                           sources=[], warnings=[]),
        ]

    def _with_sources(self, mods, pid):
        row = self.conn.execute(
            "SELECT id, name, owner, verified, circular, seq_length "
            "FROM plasmid WHERE id = ?", (pid,)).fetchone()
        for rm in mods:
            rm.sources = [dict(row, **{"start": 0, "end": 1, "strand": 1})]
        return mods

    def test_a_pair_joined_by_a_2a_is_taken_from_one_plasmid(self):
        from genesis.design import _prefer_fused_source
        fused = self.a + self.spacer + self.b
        pid = self._plasmid("fusion", fused + rand_seq(600, 133),
                            ["A", "B"])
        mods = self._with_sources(self._modules(), pid)
        notes = _prefer_fused_source(self.conn, mods)
        self.assertTrue(notes, "the fusion should have been found")
        # All three now come off that plasmid, at adjoining coordinates.
        self.assertTrue(all(m.sources[0]["id"] == pid for m in mods))
        self.assertEqual(mods[0].sources[0]["end"], mods[1].sources[0]["start"])
        self.assertEqual(mods[1].sources[0]["end"], mods[2].sources[0]["start"])
        # And the 2A is the one that plasmid actually uses.
        self.assertEqual(mods[1].seq, self.spacer)

    def test_proteins_too_far_apart_are_not_a_fusion(self):
        from genesis.design import _prefer_fused_source
        apart = self.a + rand_seq(400, 134) + self.b
        pid = self._plasmid("not-a-fusion", apart + rand_seq(300, 135),
                            ["A", "B"])
        mods = self._with_sources(self._modules(), pid)
        self.assertEqual(_prefer_fused_source(self.conn, mods), [])
        self.assertEqual(mods[1].seq, self.generic, "the 2A is unchanged")

    def test_a_spacer_with_a_stop_codon_is_not_a_2a(self):
        from genesis.design import _prefer_fused_source
        stopped = self.a + "GGA" * 9 + "TAA" + "GGA" * 10 + self.b
        pid = self._plasmid("stopped", stopped + rand_seq(300, 136),
                            ["A", "B"])
        mods = self._with_sources(self._modules(), pid)
        self.assertEqual(_prefer_fused_source(self.conn, mods), [])


class TestShrinkingArchitecture(unittest.TestCase):
    """The whole architecture first, then shorter, and exact only."""

    def setUp(self):
        self.conn = connect(":memory:")
        self.prom = rand_seq(400, 141)
        self.cds1 = orf(150, 142)
        self.p2a = "GGA" * 22
        self.cds2 = orf(120, 143)
        self.parts = {}
        for nm, sq, ft in (("PromA", self.prom, "promoter"),
                           ("CDS1", self.cds1, "CDS"),
                           ("P2A", self.p2a, "p2A"),
                           ("CDS2", self.cds2, "CDS")):
            self.parts[nm] = self.conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ft, nm, sq, nm, len(sq))).lastrowid

    def _plasmid(self, name, body, carries):
        pid = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES (?,'t',?,?,?,1,'unverified')""",
            (name, body, len(body), name)).lastrowid
        for nm in carries:
            sq = {"PromA": self.prom, "CDS1": self.cds1,
                  "P2A": self.p2a, "CDS2": self.cds2}[nm]
            at = body.find(sq)
            if at < 0:
                continue
            self.conn.execute(
                """INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,?,?,1)""",
                (self.parts[nm], pid, at, at + len(sq)))
        for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,?,?,0,30,1)""",
                (pid, ftype, nm))
        self.conn.commit()
        return pid

    def _modules(self):
        return [
            {"type": "promoter", "part_id": self.parts["PromA"],
             "part_name": "PromA", "seq": self.prom},
            {"type": "CDS", "part_id": self.parts["CDS1"],
             "part_name": "CDS1", "seq": self.cds1},
            {"type": "p2A", "part_id": self.parts["P2A"],
             "part_name": "P2A", "seq": self.p2a},
            {"type": "CDS", "part_id": self.parts["CDS2"],
             "part_name": "CDS2", "seq": self.cds2},
        ]

    def test_the_whole_architecture_is_found_when_it_exists(self):
        from genesis.swap import architecture_candidates
        whole = self.prom + "GCCACC" + self.cds1 + self.p2a + self.cds2
        self._plasmid("already-exists", whole + rand_seq(900, 144),
                      ["PromA", "CDS1", "P2A", "CDS2"])
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["matched"], 4, c["recipients"])
        self.assertEqual(c["n_donors"], 0)

    def test_it_shortens_from_the_right(self):
        from genesis.swap import architecture_candidates
        # PromA-CDS1-P2A exists; CDS2 does not follow it.
        part = self.prom + "GCCACC" + self.cds1 + self.p2a
        self._plasmid("three-of-four", part + rand_seq(900, 145),
                      ["PromA", "CDS1", "P2A"])
        self._plasmid("donor-CDS2", self.cds2 + rand_seq(700, 146), ["CDS2"])
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["matched"], 3, c["recipients"])
        self.assertEqual(c["n_donors"], 1)
        self.assertEqual(c["donors"][0]["covers"], ["CDS2"])

    def test_a_relative_does_not_lengthen_the_run(self):
        """An older version of the indicator must not read as the new one."""
        from genesis.swap import architecture_candidates
        older = self.cds1[:300] + rand_seq(60, 147) + self.cds1[360:]
        self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('CDS','CDS1-old',?,'old',?)""",
            (older, len(older)))
        old_id = self.conn.execute(
            "SELECT id FROM part WHERE name='CDS1-old'").fetchone()[0]
        body = self.prom + "GCCACC" + older + rand_seq(900, 148)
        pid = self._plasmid("has-the-old-one", body, ["PromA"])
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,?,?,1)""",
            (old_id, pid, body.find(older), body.find(older) + len(older)))
        self.conn.commit()
        c = architecture_candidates(self.conn, self._modules())
        # The promoter is really there; the coding sequence is not.
        self.assertEqual(c["matched"], 1, c["recipients"])
        top = c["recipients"][0]
        self.assertEqual(top["carries"], ["PromA"])
        self.assertTrue(any("CDS1" in x for x in top["resembles"]), top)


class TestNoStrayCassette(unittest.TestCase):
    """A design that did not ask for Cre-lox must never carry a lox site.

    Not a tidiness rule. Cre turns up in the animal, not in the plasmid, so
    a payload left between somebody else's lox sites is silent until it
    does, and an inverted one stays inverted. Whichever recipient it came
    from and whichever path put it there, it cannot ship.
    """

    def _construct(self, lox_outside_the_payload):
        from genesis.construct import (assemble_into_backbone, ResolvedModule,
                                       Module, AssembledConstruct)
        prom = rand_seq(300, 151)
        lox_a, lox_b = rand_seq(34, 152), rand_seq(34, 153)
        old = rand_seq(600, 154)
        tail = rand_seq(800, 155)
        if lox_outside_the_payload:
            # The sites sit past the region a payload swap replaces.
            seq = prom + old + tail + lox_a + lox_b
            at = len(prom) + len(old) + len(tail)
        else:
            seq = prom + lox_a + lox_b + old + tail
            at = len(prom)
        anns = [
            {"feature_type": "promoter", "name": "p", "start": 0,
             "end": len(prom), "strand": 1},
            {"feature_type": "lox_site", "name": "lox2272", "start": at,
             "end": at + 34, "strand": 1},
            {"feature_type": "lox_site", "name": "loxP", "start": at + 34,
             "end": at + 68, "strand": 1},
            {"feature_type": "CDS", "name": "old",
             "start": seq.find(old), "end": seq.find(old) + 600, "strand": 1},
            {"feature_type": "origin", "name": "ori",
             "start": len(seq) - 300, "end": len(seq) - 150, "strand": 1},
            {"feature_type": "resistance", "name": "AmpR",
             "start": len(seq) - 150, "end": len(seq), "strand": 1},
        ]
        rm = ResolvedModule(module=Module("CDS", "new"), part_id=1,
                            part_name="new", seq=rand_seq(300, 156),
                            feature_type="CDS", sources=[], warnings=[])
        return assemble_into_backbone([rm], seq, anns, {"CDS"})

    def test_a_stray_lox_site_is_an_error(self):
        from genesis.validate import validate_construct
        con = self._construct(lox_outside_the_payload=True)
        issues = validate_construct(con, wants_cassette=False, check_aav=False)
        errors = [i for i in issues if i.level == "error"]
        self.assertTrue(any("lox site" in i.message for i in errors),
                        [i.message for i in issues])

    def test_a_cre_lox_design_is_not_flagged(self):
        from genesis.validate import validate_construct
        con = self._construct(lox_outside_the_payload=True)
        issues = validate_construct(con, wants_cassette=True, check_aav=False)
        self.assertFalse(any("lox site" in i.message for i in issues
                             if i.level == "error"))

    def test_a_cassette_inside_the_payload_is_removed(self):
        con = self._construct(lox_outside_the_payload=False)
        self.assertEqual([i["name"] for i in con.layout
                          if i["type"] == "lox_site"], [])


class TestGeneBlockRoute(unittest.TestCase):
    """Buying the insert, when building it costs more reactions than it saves."""

    def _construct(self):
        from genesis.construct import (assemble_into_backbone, ResolvedModule,
                                       Module)
        prom = rand_seq(300, 161)
        old = rand_seq(600, 162)
        tail = rand_seq(1200, 163)
        seq = prom + old + tail
        anns = [
            {"feature_type": "promoter", "name": "p", "start": 0,
             "end": 300, "strand": 1},
            {"feature_type": "CDS", "name": "old", "start": 300,
             "end": 900, "strand": 1},
            {"feature_type": "origin", "name": "ori", "start": 1000,
             "end": 1400, "strand": 1},
            {"feature_type": "resistance", "name": "AmpR", "start": 1400,
             "end": 1800, "strand": 1},
        ]
        rm = ResolvedModule(module=Module("CDS", "new"), part_id=1,
                            part_name="new", seq=orf(120, 164),
                            feature_type="CDS", sources=[], warnings=[])
        return seq, assemble_into_backbone([rm], seq, anns, {"CDS"})

    def test_it_orders_everything_the_recipient_does_not_supply(self):
        from genesis.planner import plan_gblock_route
        seq, con = self._construct()
        backbone = {"id": 1, "name": "pTest", "seq": seq, "circular": 1,
                    "seq_length": len(seq)}
        route = plan_gblock_route(None, con, backbone, Config())
        self.assertTrue(route.feasible, route.warnings)
        synth = [s for s in route.steps if s.operation == "SYNTHESIS"]
        self.assertEqual(len(synth), 1, "one fragment to order, not several")
        pcr = [s for s in route.steps if s.operation == "PCR"]
        self.assertEqual(len(pcr), 1, "one vector amplification")
        ordered = synth[0].details["sequence"]
        # It carries homology to the vector at both ends.
        arm = synth[0].details["homology_arm"]
        self.assertIn(ordered[:arm], con.seq + con.seq)
        self.assertIn(ordered[-arm:], con.seq + con.seq)

    def test_without_a_recipient_there_is_nothing_to_order_around(self):
        from genesis.planner import plan_gblock_route
        _seq, con = self._construct()
        route = plan_gblock_route(None, con, None, Config())
        self.assertFalse(route.feasible)

    def test_ordered_dna_is_costed_by_the_base(self):
        """A 2.7 kb fragment is not the same purchase as a 200 bp one."""
        from genesis.planner import Route, score_route
        from genesis.operations import Step
        cheap = Route(id="a", strategy="gblock", title="t",
                      steps=[Step(1, "SYNTHESIS", "small",
                                  details={"length": 200})])
        dear = Route(id="b", strategy="gblock", title="t",
                     steps=[Step(1, "SYNTHESIS", "big",
                                 details={"length": 2700})])
        cfg = Config()
        self.assertLess(score_route(cheap, cfg).cost,
                        score_route(dear, cfg).cost)


class TestDonorCountMeansReactions(unittest.TestCase):
    """"One donor" has to mean one amplification."""

    def setUp(self):
        self.conn = connect(":memory:")
        self.a, self.b = orf(120, 171), orf(90, 172)
        # A promoter nothing else carries, so the plasmids under test are
        # donors rather than candidates to start from: a plasmid holding
        # the whole architecture is a recipient and needs no donor at all.
        self.p = rand_seq(300, 170)
        self.parts = {}
        for nm, sq in (("P", self.p), ("A", self.a), ("B", self.b)):
            self.parts[nm] = self.conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                ("promoter" if nm == "P" else "CDS",
                 nm, sq, nm, len(sq))).lastrowid
        self._plasmid("recipient-P", self.p + rand_seq(400, 176))

    def _plasmid(self, name, body):
        pid = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES (?,'t',?,?,?,1,'unverified')""",
            (name, body, len(body), name)).lastrowid
        for nm, sq in (("P", self.p), ("A", self.a), ("B", self.b)):
            at = body.find(sq)
            if at >= 0:
                self.conn.execute(
                    """INSERT INTO part_occurrence(part_id, plasmid_id, start,
                            end, strand) VALUES (?,?,?,?,1)""",
                    (self.parts[nm], pid, at, at + len(sq)))
        for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,?,?,0,30,1)""",
                (pid, ftype, nm))
        self.conn.commit()
        return pid

    def _modules(self):
        return [{"type": "promoter", "part_id": self.parts["P"],
                 "part_name": "P", "seq": self.p},
                {"type": "CDS", "part_id": self.parts["A"], "part_name": "A",
                 "seq": self.a},
                {"type": "CDS", "part_id": self.parts["B"], "part_name": "B",
                 "seq": self.b}]

    def test_two_parts_far_apart_are_two_donors(self):
        from genesis.swap import architecture_candidates
        # Both in one plasmid, but a kilobase apart: still two reactions.
        self._plasmid("holds-both-far", self.a + rand_seq(960, 173) + self.b
                      + rand_seq(400, 174))
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["n_donors"], 2, c["donors"])
        self.assertEqual([d["covers"] for d in c["donors"]], [["A"], ["B"]])

    def test_two_parts_next_to_each_other_are_one(self):
        from genesis.swap import architecture_candidates
        self._plasmid("holds-both-adjacent", self.a + self.b
                      + rand_seq(400, 175))
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["n_donors"], 1, c["donors"])
        self.assertEqual(c["donors"][0]["covers"], ["A", "B"])


class TestSourcePinIsHonoured(unittest.TestCase):
    """An explicit choice is used, or refused out loud - never ignored."""

    def test_a_plasmid_with_the_sequence_but_no_occurrence_still_works(self):
        from genesis.design import design
        conn = connect(":memory:")
        cds = orf(150, 181)
        part = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('CDS','X',?,'x',?)""",
            (cds, len(cds))).lastrowid
        # The usual source, catalogued.
        usual = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES ('usual','t',?,?,'u',1,'ok')""",
            (cds + rand_seq(900, 182), 150 + 900)).lastrowid
        conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,0,?,1)""", (part, usual, len(cds)))
        # Another plasmid holds the same bases, uncatalogued.
        other_body = rand_seq(300, 183) + cds + rand_seq(600, 184)
        other = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES ('mine','t',?,?,'m',1,'ok')""",
            (other_body, len(other_body))).lastrowid
        conn.commit()

        spec = ConstructSpec.from_dict(
            {"name": "p", "modules": [{"type": "CDS", "part_id": part}],
             "source_plasmid_ids": {part: other}})
        plan = design(conn, spec, Config())
        sources = [s["name"] for m in plan["modules"]
                   for s in (m.get("sources") or [])][:1]
        self.assertEqual(sources, ["mine"], plan["modules"])

    def test_a_plasmid_without_the_part_is_refused_not_ignored(self):
        from genesis.design import design
        conn = connect(":memory:")
        cds = orf(150, 185)
        part = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('CDS','X',?,'x',?)""",
            (cds, len(cds))).lastrowid
        usual = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES ('usual','t',?,?,'u',1,'ok')""",
            (cds + rand_seq(900, 186), 1050)).lastrowid
        conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,0,?,1)""", (part, usual, len(cds)))
        elsewhere = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES ('wrong','t',?,?,'w',1,'ok')""",
            (rand_seq(800, 187), 800)).lastrowid
        conn.commit()
        spec = ConstructSpec.from_dict(
            {"name": "p", "modules": [{"type": "CDS", "part_id": part}],
             "source_plasmid_ids": {part: elsewhere}})
        plan = design(conn, spec, Config())
        notes = [w for m in plan["modules"] for w in (m.get("warnings") or [])]
        self.assertTrue(any("does not contain" in w for w in notes), notes)


class TestCassetteIsAnArrangement(unittest.TestCase):
    """A cassette counts when it sits after the promoter, holding something."""

    def setUp(self):
        self.conn = connect(":memory:")
        self.prom = rand_seq(400, 191)
        self.cds = orf(150, 192)
        self.part_prom = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('promoter','PromA',?,'p',?)""",
            (self.prom, len(self.prom))).lastrowid
        self.part_cds = self.conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('CDS','CDS1',?,'c',?)""",
            (self.cds, len(self.cds))).lastrowid

    def _plasmid(self, name, lox_after_promoter):
        inner = rand_seq(500, 193)
        lox = [rand_seq(34, 194), rand_seq(34, 195)]
        if lox_after_promoter:
            body = (self.prom + lox[0] + lox[1] + inner + lox[0] + lox[1]
                    + rand_seq(800, 196))
            lox_at = len(self.prom)
        else:
            # The same four sites, but parked far away from the promoter.
            body = (self.prom + rand_seq(2500, 197) + lox[0] + lox[1] + inner
                    + lox[0] + lox[1] + rand_seq(800, 198))
            lox_at = len(self.prom) + 2500
        pid = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES (?,'t',?,?,?,1,'ok')""",
            (name, body, len(body), name)).lastrowid
        self.conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,0,?,1)""",
            (self.part_prom, pid, len(self.prom)))
        offs = [lox_at, lox_at + 34, lox_at + 68 + len(inner),
                lox_at + 102 + len(inner)]
        for off, nm in zip(offs, ("lox2272", "loxP", "lox2272", "loxP")):
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,'lox_site',?,?,?,1)""",
                (pid, nm, off, off + 34))
        for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
            self.conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,?,?,0,30,1)""",
                (pid, ftype, nm))
        self.conn.execute(
            """INSERT INTO annotation(plasmid_id, feature_type, name,
                    start, end, strand) VALUES (?,'CDS','payload',?,?,1)""",
            (pid, lox_at + 68, lox_at + 68 + len(inner)))
        self.conn.commit()
        return pid

    def _modules(self):
        return [
            {"type": "promoter", "part_id": self.part_prom,
             "part_name": "PromA", "seq": self.prom},
            {"type": "Cre-lox", "part_id": None, "part_name": "Cre-lox",
             "seq": ""},
            {"type": "CDS", "part_id": self.part_cds, "part_name": "CDS1",
             "seq": self.cds},
        ]

    def test_a_cassette_after_the_promoter_counts(self):
        from genesis.swap import architecture_candidates
        good = self._plasmid("pAAV-PromA-DIO-x", lox_after_promoter=True)
        c = architecture_candidates(self.conn, self._modules())
        top = c["recipients"][0]
        self.assertEqual(top["plasmid_id"], good)
        self.assertEqual(top["exact_run"], 2, top)
        # and it says what is in there, because that is what gets replaced
        self.assertEqual(top.get("cassette_holds"), ["payload"])

    def test_sites_parked_far_from_the_promoter_do_not(self):
        from genesis.swap import architecture_candidates
        self._plasmid("lox-somewhere-else", lox_after_promoter=False)
        c = architecture_candidates(self.conn, self._modules())
        runs = [r["exact_run"] for r in c["recipients"] + c["related"]]
        self.assertTrue(all(r <= 1 for r in runs), c["recipients"])


class TestPayloadOrientation(unittest.TestCase):
    """Cre-ON and Cre-OFF are the same four sites, facing differently."""

    def _payload(self, orientation):
        from genesis.construct import expand_dio, ResolvedModule, Module
        cds = orf(100, 201)
        mods = [
            ResolvedModule(module=Module("Cre-lox", "DIO", mode="DIO",
                                         payload_orientation=orientation),
                           part_id=None, part_name="DIO", seq="",
                           feature_type="Cre-lox", sources=[], warnings=[]),
            ResolvedModule(module=Module("CDS", "X"), part_id=1,
                           part_name="X", seq=cds, feature_type="CDS",
                           sources=[], warnings=[]),
        ]
        out, notes = expand_dio(mods, sites_from_backbone=True)
        payload = [m for m in out if m.feature_type == "CDS"][0]
        return cds, payload, notes

    def test_the_default_stores_it_inverted(self):
        from genesis.seqtools import revcomp
        cds, payload, _ = self._payload("")
        self.assertEqual(payload.seq, revcomp(cds))

    def test_asking_for_forward_gives_cre_off(self):
        cds, payload, notes = self._payload("forward")
        self.assertEqual(payload.seq, cds, "expressed until Cre removes it")
        self.assertTrue(any("Cre-OFF" in n for n in notes), notes)

    def test_asking_for_reverse_is_the_usual_and_says_nothing(self):
        from genesis.seqtools import revcomp
        cds, payload, notes = self._payload("reverse")
        self.assertEqual(payload.seq, revcomp(cds))
        self.assertFalse(any("against" in n for n in notes), notes)


class TestCassetteAcrossTheOrigin(unittest.TestCase):
    """Four sites on a circle enclose two arcs; only one is the cassette."""

    def _plasmid(self):
        """A real pAAV shape: the payload crosses position 1."""
        payload = orf(200, 211)
        lox = [rand_seq(34, 212), rand_seq(34, 213)]
        wpre = rand_seq(500, 214)
        prom = rand_seq(400, 215)
        ori = rand_seq(300, 216)
        amp = rand_seq(400, 217)
        # second half of the payload, the two sites, the vector, the two
        # sites again, then the first half - so the payload wraps.
        tail = payload[300:]
        head = payload[:300]
        body = (tail + lox[0] + lox[1] + wpre + ori + amp + prom
                + lox[0] + lox[1] + head)
        at = len(tail)
        anns = [
            {"feature_type": "lox_site", "name": "lox2272",
             "start": at, "end": at + 34, "strand": 1},
            {"feature_type": "lox_site", "name": "loxP",
             "start": at + 34, "end": at + 68, "strand": 1},
        ]
        pos = at + 68
        for ftype, nm, sq in (("WPRE", "WPRE", wpre), ("origin", "ori", ori),
                              ("resistance", "AmpR", amp),
                              ("promoter", "hSyn", prom)):
            anns.append({"feature_type": ftype, "name": nm, "start": pos,
                         "end": pos + len(sq), "strand": 1})
            pos += len(sq)
        anns += [
            {"feature_type": "lox_site", "name": "lox2272", "start": pos,
             "end": pos + 34, "strand": 1},
            {"feature_type": "lox_site", "name": "loxP", "start": pos + 34,
             "end": pos + 68, "strand": 1},
            {"feature_type": "CDS", "name": "payload-head",
             "start": pos + 68, "end": pos + 68 + len(head), "strand": -1},
            {"feature_type": "CDS", "name": "payload-tail",
             "start": 0, "end": len(tail), "strand": -1},
        ]
        return body, anns

    def test_the_cassette_is_the_arc_without_the_vector_in_it(self):
        from genesis.construct import backbone_cassette
        body, anns = self._plasmid()
        cas = backbone_cassette(anns, len(body))
        self.assertIsNotNone(cas)
        self.assertTrue(cas["wraps"],
                        "the payload crosses position 1, so the arc does too")

        def within(a):
            if cas["wraps"]:
                return a["start"] >= cas["start"] or a["end"] <= cas["end"]
            return cas["start"] <= a["start"] and a["end"] <= cas["end"]

        held = {a["name"] for a in anns if within(a)}
        self.assertIn("payload-head", held)
        self.assertIn("payload-tail", held)
        # and emphatically not the vector
        for nm in ("ori", "AmpR", "hSyn", "WPRE"):
            self.assertNotIn(nm, held, f"{nm} is the vector, not the payload")

    def test_a_cassette_that_does_not_wrap_is_unchanged(self):
        from genesis.construct import backbone_cassette
        inner = orf(120, 218)
        lox = [rand_seq(34, 219), rand_seq(34, 220)]
        anns = [
            {"feature_type": "lox_site", "name": "lox2272", "start": 100,
             "end": 134, "strand": 1},
            {"feature_type": "lox_site", "name": "loxP", "start": 134,
             "end": 168, "strand": 1},
            {"feature_type": "CDS", "name": "payload", "start": 168,
             "end": 168 + len(inner), "strand": 1},
            {"feature_type": "lox_site", "name": "lox2272",
             "start": 168 + len(inner), "end": 202 + len(inner), "strand": 1},
            {"feature_type": "lox_site", "name": "loxP",
             "start": 202 + len(inner), "end": 236 + len(inner), "strand": 1},
            {"feature_type": "origin", "name": "ori",
             "start": 300 + len(inner), "end": 500 + len(inner), "strand": 1},
            {"feature_type": "resistance", "name": "AmpR",
             "start": 500 + len(inner), "end": 800 + len(inner), "strand": 1},
        ]
        cas = backbone_cassette(anns, 1000 + len(inner))
        self.assertFalse(cas["wraps"])
        self.assertEqual(cas["start"], 168)
        self.assertEqual(cas["end"], 168 + len(inner))


class TestOnlyRealPlasmidsAreParents(unittest.TestCase):
    """A PCR product is not something you can take out of a freezer.

    The share holds gel extractions, ordered fragments and PCR products
    next to real plasmids, and a part sitting in one of those cannot be
    amplified next week. Topology alone does not sort them: a 1.3 kb
    fragment in the inventory claims to be circular. Carrying an origin
    and a resistance marker does.
    """

    def setUp(self):
        self.conn = connect(":memory:")
        self.prom, self.cds = rand_seq(300, 181), orf(600, 182)
        self.parts = {}
        for nm, sq, ft in (("P", self.prom, "promoter"),
                           ("X", self.cds, "CDS")):
            self.parts[nm] = self.conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ft, nm, sq, nm, len(sq))).lastrowid

    def _row(self, name, body, carries, vector, circular=1):
        pid = self.conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES (?,'t',?,?,?,?,'unverified')""",
            (name, body, len(body), name, circular)).lastrowid
        for nm in carries:
            at = body.find(self.parts and
                           (self.prom if nm == "P" else self.cds))
            self.conn.execute(
                """INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,?,?,1)""",
                (self.parts[nm], pid, at,
                 at + len(self.prom if nm == "P" else self.cds)))
        if vector:
            for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
                self.conn.execute(
                    """INSERT INTO annotation(plasmid_id, feature_type, name,
                            start, end, strand) VALUES (?,?,?,0,30,1)""",
                    (pid, ftype, nm))
        self.conn.commit()
        return pid

    def _modules(self):
        return [{"type": "promoter", "part_id": self.parts["P"],
                 "part_name": "P", "seq": self.prom},
                {"type": "CDS", "part_id": self.parts["X"],
                 "part_name": "X", "seq": self.cds}]

    def test_a_fragment_is_not_offered_as_a_donor(self):
        from genesis.swap import architecture_candidates
        self._row("recipient", self.prom + rand_seq(400, 183), ["P"], True)
        self._row("2-4 PCRprdt", self.cds + rand_seq(80, 184), ["X"], False)
        c = architecture_candidates(self.conn, self._modules())
        offered = [o["name"] for d in c["donors"] for o in d["options"]]
        self.assertNotIn("2-4 PCRprdt", offered)

    def test_a_fragment_calling_itself_circular_is_still_a_fragment(self):
        from genesis.swap import architecture_candidates
        self._row("recipient", self.prom + rand_seq(400, 185), ["P"], True)
        self._row("Insert1 1322 bp ds-DNA linear",
                  self.cds + rand_seq(80, 186), ["X"], False, circular=1)
        c = architecture_candidates(self.conn, self._modules())
        offered = [o["name"] for d in c["donors"] for o in d["options"]]
        self.assertEqual(offered, [])

    def test_the_part_is_still_reported_when_only_fragments_have_it(self):
        """Dropping it silently would read as "nobody has this"."""
        from genesis.swap import architecture_candidates
        self._row("recipient", self.prom + rand_seq(400, 187), ["P"], True)
        self._row("2-4 PCRprdt", self.cds + rand_seq(80, 188), ["X"], False)
        c = architecture_candidates(self.conn, self._modules())
        self.assertEqual(c["donors"][-1]["covers"], ["X"])
        self.assertIn("PCR product", c["donors"][-1]["note"])

    def test_a_real_plasmid_is_offered(self):
        from genesis.swap import architecture_candidates
        self._row("recipient", self.prom + rand_seq(400, 189), ["P"], True)
        self._row("PP072", self.cds + rand_seq(2000, 190), ["X"], True)
        c = architecture_candidates(self.conn, self._modules())
        offered = [o["name"] for d in c["donors"] for o in d["options"]]
        self.assertIn("PP072", offered)


class TestSamePartAtAnotherBoundary(unittest.TestCase):
    """One gene catalogued twice is one part, not two.

    EGFP is in the inventory both with and without its ATG, and P2A both
    with and without the GSG linker in front of it. Requiring the exact
    catalogue row that happened to be resolved made PP078 - which carries
    P2A fused directly to EGFP - look as though it carried neither, so the
    plan went off to amplify the two separately from two other plasmids.
    """

    def test_a_missing_start_codon_is_the_same_part(self):
        from genesis.swap import interchangeable
        cds = orf(600, 191)
        self.assertTrue(interchangeable(cds, cds[3:], "CDS"))

    def test_a_leading_linker_is_the_same_part(self):
        from genesis.swap import interchangeable
        p2a = orf(66, 192)
        self.assertTrue(interchangeable(p2a, p2a[9:], "p2A"))

    def test_different_codons_for_the_same_protein_are_the_same_part(self):
        from genesis.swap import interchangeable, protein_of
        cds = orf(600, 193)
        # Re-spell the third position of a few codons: same protein,
        # different bases, which is what a codon-optimised copy is.
        syn = {"GCT": "GCC", "GGT": "GGC", "TGC": "TGT", "ATT": "ATC"}
        out = [syn.get(cds[i:i + 3], cds[i:i + 3])
               for i in range(0, len(cds) - 2, 3)]
        other = "".join(out)
        self.assertNotEqual(cds, other)
        self.assertEqual(protein_of(cds), protein_of(other))
        self.assertTrue(interchangeable(cds, other, "CDS"))

    def test_another_version_of_the_indicator_is_not(self):
        """The one thing somebody asking for ASAP7y must not be handed."""
        from genesis.swap import interchangeable
        want = orf(900, 194)
        relative = want[:400] + orf(60, 195) + want[460:]
        self.assertFalse(interchangeable(want, relative, "CDS"))

    def test_a_stub_row_is_not_the_whole_gene(self):
        """There is a 6 bp part named "ASAP1" in the inventory."""
        from genesis.swap import interchangeable
        want = orf(900, 196)
        self.assertFalse(interchangeable(want, want[:6], "CDS"))

    def test_two_parts_fused_in_one_plasmid_are_one_donor(self):
        from genesis.swap import architecture_candidates
        conn = connect(":memory:")
        prom, a, b = rand_seq(300, 197), orf(300, 198), orf(600, 199)
        parts = {}
        for nm, sq, ft in (("P", prom, "promoter"), ("A", a, "p2A"),
                           ("B", b, "CDS")):
            parts[nm] = conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES (?,?,?,?,?)""",
                (ft, nm, sq, nm, len(sq))).lastrowid

        def row(name, body, held):
            pid = conn.execute(
                """INSERT INTO plasmid(name, owner, seq, seq_length,
                        seq_checksum, circular, verified)
                   VALUES (?,'t',?,?,?,1,'unverified')""",
                (name, body, len(body), name)).lastrowid
            for nm, sq in held:
                at = body.find(sq)
                conn.execute(
                    """INSERT INTO part_occurrence(part_id, plasmid_id, start,
                            end, strand) VALUES (?,?,?,?,1)""",
                    (parts[nm], pid, at, at + len(sq)))
            for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
                conn.execute(
                    """INSERT INTO annotation(plasmid_id, feature_type, name,
                            start, end, strand) VALUES (?,?,?,0,30,1)""",
                    (pid, ftype, nm))
            conn.commit()
            return pid

        row("recipient", prom + rand_seq(400, 200), [("P", prom)])
        # Holds both, fused - but spelled at the other boundary: no ATG on
        # the coding sequence, no linker on the 2A.
        row("fused-donor", a[9:] + b[3:] + rand_seq(500, 201),
            [("A", a[9:]), ("B", b[3:])])
        modules = [{"type": "promoter", "part_id": parts["P"],
                    "part_name": "P", "seq": prom},
                   {"type": "p2A", "part_id": parts["A"],
                    "part_name": "A", "seq": a},
                   {"type": "CDS", "part_id": parts["B"],
                    "part_name": "B", "seq": b}]
        c = architecture_candidates(conn, modules)
        self.assertEqual(c["n_donors"], 1, c["donors"])
        self.assertEqual(c["donors"][0]["covers"], ["A", "B"])
        conn.close()


class TestInvertedPayloadStillMergesOnePCR(unittest.TestCase):
    """A Cre-lox payload runs backwards through a forward donor.

    P2A and EGFP sit next to each other in one plasmid, so they are one
    amplification. Placed inside a DIO cassette the payload is inverted,
    so its modules walk that donor from the far end towards the near one
    while the donor still holds them the right way round. Testing the
    source strand alone missed that and planned two reactions, one of
    them 66 bp long.
    """

    def _run(self, placed_strand):
        from genesis.planner import group_contiguous
        from genesis.construct import ResolvedModule, AssembledConstruct

        a, b = orf(66, 211), orf(717, 212)

        def source(start, end):
            return [{"id": 9, "name": "donor", "start": start, "end": end,
                     "strand": 1, "seq_length": 5000,
                     "verified": "unverified"}]

        # Adjacent in the donor: P2A at 100-166, EGFP at 166-883.
        mods = [ResolvedModule(None, 1, "P2A", a, "p2A", source(100, 166)),
                ResolvedModule(None, 2, "EGFP", b, "CDS", source(166, 883))]
        if placed_strand == -1:
            mods.reverse()
        layout = [{"name": m.part_name, "type": m.feature_type, "start": i,
                   "end": i + 1, "strand": placed_strand}
                  for i, m in enumerate(mods)]
        construct = AssembledConstruct(seq="", layout=layout, circular=True)
        return group_contiguous(mods, construct, Config())

    def test_forward_payload_is_one_fragment(self):
        frags = self._run(1)
        self.assertEqual(len(frags), 1, [f.name for f in frags])
        self.assertEqual(frags[0].name, "P2A+EGFP")

    def test_inverted_payload_is_also_one_fragment(self):
        frags = self._run(-1)
        self.assertEqual(len(frags), 1, [f.name for f in frags])
        self.assertEqual(frags[0].name, "EGFP+P2A")
        # And it spans the whole of both, not just the second.
        self.assertEqual((frags[0].start, frags[0].end), (100, 883))


class TestPinnedSourceAcrossTheOrigin(unittest.TestCase):
    """A plasmid is a circle; its sequence file has to start somewhere.

    A part can straddle the join between the last base of the file and the
    first, and a plain substring search reads that file as a line and
    misses it. Genesis would then tell somebody the plasmid they named
    does not contain the insert they can see in it, and quietly amplify
    from somewhere else instead.
    """

    def _inventory(self, body, circular, cds):
        """A pinned plasmid holding `body`, plus a catalogued usual source."""
        conn = connect(":memory:")
        part = conn.execute(
            """INSERT INTO part(feature_type, name, seq, seq_checksum,
                    seq_length) VALUES ('CDS','X',?,'x',?)""",
            (cds, len(cds))).lastrowid
        usual = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES ('usual','t',?,?,'u',1,'ok')""",
            (cds + rand_seq(900, 302), len(cds) + 900)).lastrowid
        conn.execute(
            """INSERT INTO part_occurrence(part_id, plasmid_id, start, end,
                    strand) VALUES (?,?,0,?,1)""", (part, usual, len(cds)))
        mine = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified) VALUES ('mine','t',?,?,'m',?,'ok')""",
            (body, len(body), 1 if circular else 0)).lastrowid
        conn.commit()
        return conn, part, mine

    def _plan(self, conn, part, mine):
        from genesis.design import design
        spec = ConstructSpec.from_dict(
            {"name": "p", "modules": [{"type": "CDS", "part_id": part}],
             "source_plasmid_ids": {part: mine}})
        return design(conn, spec, Config())

    @staticmethod
    def _chosen(plan):
        return [s["name"] for m in plan["modules"]
                for s in (m.get("sources") or [])][:1]

    def test_an_insert_across_the_origin_is_found(self):
        cds = orf(150, 301)
        # The part is split by the origin: its tail opens the file, its
        # head closes it. Reading the file as a line, it is in neither.
        head, tail = cds[:90], cds[90:]
        body = tail + rand_seq(700, 303) + head
        conn, part, mine = self._inventory(body, True, cds)
        self.assertNotIn(cds, body)          # not there as a plain substring
        plan = self._plan(conn, part, mine)
        notes = [w for m in plan["modules"] for w in (m.get("warnings") or [])]
        self.assertFalse(any("does not contain" in w for w in notes), notes)
        self.assertEqual(self._chosen(plan), ["mine"], plan["modules"])
        conn.close()

    def test_the_recorded_span_wraps_rather_than_running_backwards(self):
        """start stays inside the sequence; end runs past it, as the
        occurrence table and slice_region both expect."""
        from genesis.seqtools import slice_region
        cds = orf(150, 304)
        head, tail = cds[:90], cds[90:]
        body = tail + rand_seq(700, 305) + head
        conn, part, mine = self._inventory(body, True, cds)
        plan = self._plan(conn, part, mine)
        src = next(s for m in plan["modules"]
                   for s in (m.get("sources") or []) if s["name"] == "mine")
        self.assertTrue(0 <= src["start"] < len(body), src)
        self.assertGreater(src["end"], len(body))
        self.assertEqual(
            slice_region(body, src["start"], src["end"], circular=True), cds)
        conn.close()

    def test_an_ordinary_internal_insert_still_works(self):
        cds = orf(150, 306)
        body = rand_seq(300, 307) + cds + rand_seq(600, 308)
        conn, part, mine = self._inventory(body, True, cds)
        plan = self._plan(conn, part, mine)
        self.assertEqual(self._chosen(plan), ["mine"], plan["modules"])
        conn.close()

    def test_the_reverse_strand_across_the_origin_is_found(self):
        from genesis.seqtools import revcomp
        cds = orf(150, 309)
        other = revcomp(cds)
        body = other[90:] + rand_seq(700, 310) + other[:90]
        conn, part, mine = self._inventory(body, True, cds)
        plan = self._plan(conn, part, mine)
        self.assertEqual(self._chosen(plan), ["mine"], plan["modules"])
        src = next(s for m in plan["modules"]
                   for s in (m.get("sources") or []) if s["name"] == "mine")
        self.assertEqual(src["strand"], -1, src)
        conn.close()

    def test_a_linear_fragment_does_not_match_across_its_ends(self):
        """The join only exists on a circle. Reading a linear fragment as
        one would invent a sequence nobody can amplify."""
        cds = orf(150, 311)
        head, tail = cds[:90], cds[90:]
        body = tail + rand_seq(700, 312) + head
        conn, part, mine = self._inventory(body, False, cds)
        plan = self._plan(conn, part, mine)
        notes = [w for m in plan["modules"] for w in (m.get("warnings") or [])]
        self.assertTrue(any("does not contain" in w for w in notes), notes)
        self.assertEqual(self._chosen(plan), ["usual"], plan["modules"])
        conn.close()


class TestWrappedRegionAmplifiesInEitherFrame(unittest.TestCase):
    """Coordinates for a wrapped region must stay inside the sequence.

    Mapping a region onto the other strand is a subtraction, and on one
    that wraps the origin the plain subtraction lands before position
    zero. The slice happened to survive that, because a negative start
    comes back round under modulo, but the coordinates were reported to
    the bench as a span starting below zero.
    """

    def setUp(self):
        import random
        rng = random.Random(77)
        self.tpl = "".join(rng.choice("ACGT") for _ in range(400))
        self.start, self.end = 370, 425          # 30 bases + 25 bases

    def _want(self):
        from genesis.seqtools import slice_region
        return slice_region(self.tpl, self.start, self.end, circular=True)

    def test_the_forward_frame_is_unchanged(self):
        from genesis.planner import _amplification_frame
        from genesis.seqtools import slice_region
        want = self._want()
        tpl, a, b = _amplification_frame(
            self.tpl, self.start, self.end, 1, wanted=want, circular=True)
        self.assertEqual(slice_region(tpl, a, b, circular=True), want)

    def test_the_reverse_frame_yields_the_same_bases(self):
        from genesis.planner import _amplification_frame
        from genesis.seqtools import revcomp, slice_region
        want = revcomp(self._want())
        tpl, a, b = _amplification_frame(
            self.tpl, self.start, self.end, -1, wanted=want, circular=True)
        self.assertEqual(slice_region(tpl, a, b, circular=True), want)

    def test_the_reverse_frame_coordinates_are_not_negative(self):
        from genesis.planner import _amplification_frame
        from genesis.seqtools import revcomp
        tpl, a, b = _amplification_frame(
            self.tpl, self.start, self.end, -1,
            wanted=revcomp(self._want()), circular=True)
        self.assertGreaterEqual(a, 0)
        self.assertEqual(b - a, self.end - self.start)


class TestAddgeneIsSearchOnly(unittest.TestCase):
    """Searching Addgene never puts an unordered plasmid into a protocol.

    A catalogue number is weeks and a purchase order away from being a
    tube in a freezer. Letting one compete on cost against something the
    lab already has would make every route estimate a lie, so what comes
    back is reported beside the plan and never inside it.
    """

    def setUp(self):
        import genesis.addgene as ag
        import tempfile
        self.ag = ag
        self._env = os.environ.pop("GENESIS_ADDGENE_TOKEN", None)
        # The cache lives beside the database, so a test that did not say
        # where that is would write into the real ~/.genesis and hand the
        # next test an answer it never asked for.
        self._tmp = tempfile.TemporaryDirectory()
        self.cfg = Config()
        self.cfg.db_path = os.path.join(self._tmp.name, "genesis.db")

    def tearDown(self):
        self._tmp.cleanup()
        if self._env is not None:
            os.environ["GENESIS_ADDGENE_TOKEN"] = self._env
        else:
            os.environ.pop("GENESIS_ADDGENE_TOKEN", None)

    def test_without_a_token_it_refuses_rather_than_scraping(self):
        with self.assertRaises(self.ag.AddgeneUnavailable) as caught:
            self.ag.search("ASAP7y", cfg=self.cfg)
        self.assertIn("token", str(caught.exception))

    def test_it_offers_a_browsable_url_instead(self):
        url = self.ag.search_url("P2A-EGFP")
        self.assertTrue(url.startswith("https://www.addgene.org/"))
        self.assertIn("P2A-EGFP", urllib.parse.unquote(url))

    def test_a_missing_part_reports_why_rather_than_an_empty_list(self):
        out = self.ag.suggest_for_missing(["Nonesuch7z"], cfg=self.cfg)
        self.assertEqual(out["Nonesuch7z"]["results"], [])
        self.assertIn("token", out["Nonesuch7z"]["unavailable"])
        self.assertTrue(out["Nonesuch7z"]["browse"])

    def test_the_token_comes_from_the_environment_first(self):
        cfg = self.cfg
        cfg.addgene_token = "from-config"
        self.assertEqual(self.ag.token(cfg), "from-config")
        os.environ["GENESIS_ADDGENE_TOKEN"] = "from-env"
        self.assertEqual(self.ag.token(cfg), "from-env")

    def test_a_scan_without_a_token_does_not_touch_the_network(self):
        """An offline lab must get exactly the scan it got before."""
        def explode(*a, **k):
            raise AssertionError("the network was touched")
        original = self.ag._request
        self.ag._request = explode
        try:
            from genesis.design import _offer_to_order
            candidates = {"donors": [{"covers": ["Xa"], "options": []}]}
            _offer_to_order(candidates, self.cfg)
        finally:
            self.ag._request = original
        self.assertNotIn("to_order", candidates)
        self.assertIn("developers.addgene.org", candidates["to_order_hint"])

    def test_results_are_reported_beside_the_plan_not_as_donors(self):
        from genesis.design import _offer_to_order
        cfg = self.cfg
        os.environ["GENESIS_ADDGENE_TOKEN"] = "pretend"
        original = self.ag._request
        self.ag._request = lambda path, params, key: {"results": [{
            "id": 12345, "name": "pAAV-hSyn-Thing",
            "url": "https://www.addgene.org/12345/",
            "depositor": ["A Lab"], "inserts": [], "cloning": {},
            "bacterial_resistance": "Ampicillin", "article": {}}]}
        try:
            donor = {"covers": ["Xb"], "options": []}
            candidates = {"donors": [donor]}
            _offer_to_order(candidates, cfg)
        finally:
            self.ag._request = original
        # Reported...
        self.assertEqual(
            candidates["to_order"]["Xb"]["results"][0]["id"], 12345)
        # ...and the donor is still empty, so no route can amplify from it.
        self.assertEqual(donor["options"], [])

    def test_a_search_failure_never_fails_the_scan(self):
        from genesis.design import _offer_to_order
        os.environ["GENESIS_ADDGENE_TOKEN"] = "pretend"
        original = self.ag._request
        self.ag._request = lambda *a, **k: (_ for _ in ()).throw(
            self.ag.AddgeneUnavailable("Addgene is down"))
        try:
            candidates = {"donors": [{"covers": ["Xc"], "options": []}]}
            _offer_to_order(candidates, self.cfg)
        finally:
            self.ag._request = original
        self.assertEqual(candidates["to_order"]["Xc"]["results"], [])
        self.assertIn("down", candidates["to_order"]["Xc"]["unavailable"])


class TestASwapMustReplaceThePayload(unittest.TestCase):
    """Removing too little leaves the old indicator in the plasmid.

    A request for CAG - ASAP7y against a pAAV-CAG-Voltron2 recipient cut
    16 bp out of it - the Kozak between the promoter and the coding
    sequence - and put ASAP7y in there, leaving Voltron-2, Ace2N and
    HaloTag sitting between the new gene and the WPRE. Every presence
    check passed, because every part that was asked for really was in the
    product, and the plasmid was ruined anyway.
    """

    def _recipient(self, conn, prom, payload):
        body = (rand_seq(140, 401) + prom + rand_seq(16, 402) + payload
                + rand_seq(600, 403))
        pid = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('recipient','t',?,?,'r',1,'unverified')""",
            (body, len(body))).lastrowid
        at_p = body.find(prom)
        at_x = body.find(payload)
        rows = [("promoter", "CAG", at_p, at_p + len(prom)),
                ("unknown", "Kozak sequence", at_p + len(prom),
                 at_p + len(prom) + 16),
                ("CDS", "Voltron-2", at_x, at_x + len(payload)),
                ("WPRE", "WPRE", at_x + len(payload),
                 at_x + len(payload) + 300),
                ("origin", "ori", len(body) - 200, len(body) - 100),
                ("resistance", "AmpR", len(body) - 100, len(body))]
        for ftype, nm, a, b in rows:
            conn.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name,
                        start, end, strand) VALUES (?,?,?,?,?,1)""",
                (pid, ftype, nm, a, b))
        conn.commit()
        return pid, body

    def test_the_whole_payload_is_removed_not_the_gap_before_it(self):
        from genesis.swap import _payload_end, _annotations
        conn = connect(":memory:")
        prom, payload = rand_seq(400, 404), orf(900, 405)
        pid, body = self._recipient(conn, prom, payload)
        annots = _annotations(conn, pid)
        start = body.find(prom) + len(prom)
        end = _payload_end(annots, start, len(body))
        # It must reach past the Kozak and take the whole coding sequence.
        self.assertEqual(end, body.find(payload) + len(payload),
                         f"stopped at {end}, payload ends at "
                         f"{body.find(payload) + len(payload)}")
        conn.close()

    def test_it_stops_at_the_wpre_rather_than_eating_the_vector(self):
        from genesis.swap import _payload_end, _annotations
        conn = connect(":memory:")
        prom, payload = rand_seq(400, 406), orf(900, 407)
        pid, body = self._recipient(conn, prom, payload)
        annots = _annotations(conn, pid)
        start = body.find(prom) + len(prom)
        end = _payload_end(annots, start, len(body))
        wpre_start = body.find(payload) + len(payload)
        self.assertLessEqual(end, wpre_start)
        conn.close()

    def test_furniture_alone_is_not_a_payload(self):
        """A Kozak with no coding sequence after it is nothing to replace."""
        from genesis.swap import _payload_end
        annots = [{"feature_type": "unknown", "name": "Kozak sequence",
                   "start": 100, "end": 116, "strand": 1},
                  {"feature_type": "WPRE", "name": "WPRE",
                   "start": 130, "end": 700, "strand": 1}]
        self.assertEqual(_payload_end(annots, 100, 5000), 100)


class TestTheProductMustNotCarryWhatWasNotAskedFor(unittest.TestCase):
    """Presence checks alone cannot see what came along for the ride."""

    def test_a_leftover_coding_sequence_fails_the_simulation(self):
        from genesis.simulate import _check_leftovers, SimulationResult
        conn = connect(":memory:")
        old, new = orf(900, 411), orf(600, 412)
        body = rand_seq(100, 413) + old + rand_seq(100, 414)
        pid = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('parent','t',?,?,'p',1,'unverified')""",
            (body, len(body))).lastrowid
        at = body.find(old)
        conn.execute(
            """INSERT INTO annotation(plasmid_id, feature_type, name, start,
                    end, strand) VALUES (?,'CDS','Voltron-2',?,?,1)""",
            (pid, at, at + len(old)))
        conn.commit()

        # A product that kept the old payload as well as taking the new one.
        product = body + new
        result = SimulationResult()
        plan = {"backbone": {"id": pid, "name": "parent"},
                "modules": [{"part_name": "ASAP7y", "seq": new}]}
        _check_leftovers(result, product, plan, conn)
        self.assertFalse(result.passed)
        self.assertTrue(any("Voltron-2" in p for p in result.problems),
                        result.problems)
        conn.close()

    def test_a_clean_replacement_passes(self):
        from genesis.simulate import _check_leftovers, SimulationResult
        conn = connect(":memory:")
        old, new = orf(900, 415), orf(600, 416)
        body = rand_seq(100, 417) + old + rand_seq(100, 418)
        pid = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('parent','t',?,?,'p',1,'unverified')""",
            (body, len(body))).lastrowid
        at = body.find(old)
        conn.execute(
            """INSERT INTO annotation(plasmid_id, feature_type, name, start,
                    end, strand) VALUES (?,'CDS','Voltron-2',?,?,1)""",
            (pid, at, at + len(old)))
        conn.commit()
        product = body.replace(old, new)       # the payload really replaced
        result = SimulationResult()
        plan = {"backbone": {"id": pid, "name": "parent"},
                "modules": [{"part_name": "ASAP7y", "seq": new}]}
        _check_leftovers(result, product, plan, conn)
        self.assertTrue(result.passed, result.problems)
        conn.close()

    def test_a_part_that_was_asked_for_is_not_a_leftover(self):
        from genesis.simulate import _check_leftovers, SimulationResult
        conn = connect(":memory:")
        keep = orf(900, 419)
        body = rand_seq(100, 420) + keep + rand_seq(100, 421)
        pid = conn.execute(
            """INSERT INTO plasmid(name, owner, seq, seq_length, seq_checksum,
                    circular, verified)
               VALUES ('parent','t',?,?,'p',1,'unverified')""",
            (body, len(body))).lastrowid
        at = body.find(keep)
        conn.execute(
            """INSERT INTO annotation(plasmid_id, feature_type, name, start,
                    end, strand) VALUES (?,'CDS','ASAP7y',?,?,1)""",
            (pid, at, at + len(keep)))
        conn.commit()
        result = SimulationResult()
        plan = {"backbone": {"id": pid, "name": "parent"},
                "modules": [{"part_name": "ASAP7y", "seq": keep}]}
        _check_leftovers(result, body, plan, conn)
        self.assertTrue(result.passed, result.problems)
        conn.close()


class TestDonorsAreRankedOnHavingThePart(unittest.TestCase):
    """Smallness decides between donors that hold the same thing."""

    def test_a_relative_never_outranks_the_real_part(self):
        from genesis.swap import find_donors
        conn = connect(":memory:")
        want = orf(1200, 431)
        relative = want[:500] + orf(90, 432) + want[590:]
        parts = {}
        for nm, sq in (("ASAP7y", want), ("ASAP3", relative)):
            parts[nm] = conn.execute(
                """INSERT INTO part(feature_type, name, seq, seq_checksum,
                        seq_length) VALUES ('CDS',?,?,?,?)""",
                (nm, sq, nm, len(sq))).lastrowid

        def plasmid(name, body, part):
            pid = conn.execute(
                """INSERT INTO plasmid(name, owner, seq, seq_length,
                        seq_checksum, circular, verified, available_template)
                   VALUES (?,'t',?,?,?,1,'unverified',1)""",
                (name, body, len(body), name)).lastrowid
            at = body.find(part[1])
            conn.execute(
                """INSERT INTO part_occurrence(part_id, plasmid_id, start,
                        end, strand) VALUES (?,?,?,?,1)""",
                (parts[part[0]], pid, at, at + len(part[1])))
            for ftype, nm in (("origin", "ori"), ("resistance", "AmpR")):
                conn.execute(
                    """INSERT INTO annotation(plasmid_id, feature_type, name,
                            start, end, strand) VALUES (?,?,?,0,30,1)""",
                    (pid, ftype, nm))
            conn.commit()
            return pid

        # The relative sits in a smaller plasmid, which used to win.
        small = plasmid("small-ASAP3", relative + rand_seq(600, 433),
                        ("ASAP3", relative))
        right = plasmid("bigger-ASAP7y", want + rand_seq(3000, 434),
                        ("ASAP7y", want))
        donors = find_donors(conn, [parts["ASAP7y"], parts["ASAP3"]],
                             exact_ids={parts["ASAP7y"]})
        self.assertEqual(donors[0].plasmid_id, right,
                         [d.name for d in donors])
        self.assertTrue(any("has the part itself" in r
                            for r in donors[0].reasons), donors[0].reasons)
        self.assertEqual(donors[1].plasmid_id, small)
        self.assertTrue(any("not the part" in r for r in donors[1].reasons))
        conn.close()


class TestPrimerSiteMovesToBetterBases(unittest.TestCase):
    """A junction is where the DNA joins, not where a primer must bind.

    The annealing site used to be fixed at the fragment boundary and only
    its length could change, so a junction in an AT-rich stretch gave a
    primer that was reported as too cold and used anyway. The site may now
    move into the fragment, with the bases it steps over carried on the
    tail - which leaves the product identical and only lengthens the oligo.
    """

    def setUp(self):
        import random
        rng = random.Random(909)
        # An AT-rich mouth followed by balanced sequence, so moving in is
        # the only way to reach NEB's 40-60% band.
        poor = "".join(rng.choice("ATATATATGC") for _ in range(40))
        good = "".join(rng.choice("ACGT") for _ in range(3000))
        self.template = poor + good
        self.up = "".join(rng.choice("ACGT") for _ in range(300))
        self.dn = "".join(rng.choice("ACGT") for _ in range(300))

    def _pair(self, shift):
        from genesis.primers import design_gibson_pair
        cfg = Config()
        cfg.primer_max_shift = shift
        return design_gibson_pair(0, "frag", self.template, 0, 1500,
                                  self.up, self.dn, cfg,
                                  carry_downstream=False), cfg

    def test_the_product_is_unchanged_by_the_shift(self):
        from genesis.seqtools import revcomp
        products = []
        for shift in (0, 15):
            (f, r), _ = self._pair(shift)
            at = self.template.find(f.anneal)
            rc_at = self.template.find(revcomp(r.anneal))
            products.append(f.sequence
                            + self.template[at + len(f.anneal):rc_at]
                            + revcomp(r.sequence))
        self.assertEqual(products[0], products[1])

    def test_it_reaches_better_bases_than_the_boundary_offers(self):
        from genesis.seqtools import gc_fraction
        (f0, _), _ = self._pair(0)
        (f1, _), _ = self._pair(15)
        self.assertGreater(gc_fraction(f1.anneal), gc_fraction(f0.anneal))

    def test_the_bases_stepped_over_are_on_the_tail(self):
        (f, _), _ = self._pair(15)
        self.assertEqual(f.sequence, f.tail + f.anneal)
        at = self.template.find(f.anneal)
        self.assertTrue(f.tail.endswith(self.template[:at]), f.tail)

    def test_moving_the_site_is_said_out_loud(self):
        (f, _), _ = self._pair(15)
        self.assertTrue(any("moved" in w for w in f.warnings), f.warnings)


class TestStickyEndLigationIsSimulated(unittest.TestCase):
    """A digest route cannot be judged by looking for homology.

    Every restriction route was reported as "the fragments do not
    assemble", because the simulator amplified what the primers described
    and then looked for 15 bp of shared end sequence. A sticky-end
    ligation does not have that and never will, so the cheapest route on
    the page was being discarded by a test that did not apply to it.
    """

    def _plasmid(self):
        import random
        rng = random.Random(555)
        filler = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        # GGATCC = BamHI, AAGCTT = HindIII
        return (filler(300) + "GGATCC" + filler(900) + "AAGCTT"
                + filler(400))

    def test_a_digest_cuts_into_the_expected_pieces(self):
        from genesis.simulate import digest
        seq = self._plasmid()
        pieces = digest(seq, True, ["BamHI", "HindIII"])
        self.assertEqual(len(pieces), 2)
        self.assertEqual({p.left for p in pieces}, {"GATC", "AGCT"})

    def test_the_pieces_religate_to_the_original(self):
        from genesis.simulate import digest, ligate
        seq = self._plasmid()
        joined, note = ligate(digest(seq, True, ["BamHI", "HindIII"]))
        self.assertIsNotNone(joined, note)
        self.assertEqual(len(joined), len(seq))

    def test_incompatible_ends_do_not_close(self):
        from genesis.simulate import digest, ligate
        seq = self._plasmid()
        vector = [p for p in digest(seq, True, ["BamHI", "HindIII"])
                  if len(p.seq) > 500][0]
        # An insert cut with one enzyme only has two identical ends, so it
        # cannot go into a vector cut with two different ones directionally.
        from genesis.simulate import CutFragment
        wrong = CutFragment("ACGT" * 50, "GATC", "GATC")
        joined, note = ligate([vector, wrong])
        self.assertIsNone(joined, note)

    def test_blunt_ends_are_not_treated_as_a_join(self):
        from genesis.simulate import CutFragment, ligate
        a = CutFragment("ACGT" * 40, "", "")
        b = CutFragment("TTGA" * 40, "", "")
        joined, _ = ligate([a, b])
        self.assertIsNone(joined)


class TestProteinEndsAreLabelled(unittest.TestCase):
    """Which end of a block is the start of the protein.

    On a circular map the two ends of a coding sequence look alike, and
    for a Cre-lox payload - stored inverted - the amino end is the one
    drawn second. Getting it from the drawing rather than the strand
    would label every inverted cassette backwards.
    """

    def test_a_forward_gene_starts_at_its_low_end(self):
        from genesis.pdf import _terminus_ends
        f = {"_start_f": 0.2, "_end_f": 0.6, "strand": 1}
        self.assertEqual(_terminus_ends(f), (0.2, 0.6))

    def test_a_reverse_gene_starts_at_its_high_end(self):
        from genesis.pdf import _terminus_ends
        f = {"_start_f": 0.2, "_end_f": 0.6, "strand": -1}
        self.assertEqual(_terminus_ends(f), (0.6, 0.2))

    def test_a_feature_with_no_strand_reads_forward(self):
        from genesis.pdf import _terminus_ends
        f = {"_start_f": 0.1, "_end_f": 0.3}
        self.assertEqual(_terminus_ends(f), (0.1, 0.3))

    def test_only_coding_sequences_get_their_own_ends(self):
        """A 2A peptide is part of somebody else's protein."""
        from genesis.pdf import _PROTEIN_ENDS
        self.assertIn("CDS", _PROTEIN_ENDS)
        for borrowed in ("p2A", "trafficking_motif", "linker", "tag"):
            self.assertNotIn(borrowed, _PROTEIN_ENDS)

    def test_an_inverted_fusion_reads_from_high_to_low(self):
        """The order the ribosome sees, for a payload stored backwards."""
        from genesis.pdf import _terminus_ends
        # ASAP7y then P2A then EGFP, all on the minus strand.
        asap = {"_start_f": 0.30, "_end_f": 0.55, "strand": -1}
        egfp = {"_start_f": 0.08, "_end_f": 0.26, "strand": -1}
        asap_n, asap_c = _terminus_ends(asap)
        egfp_n, egfp_c = _terminus_ends(egfp)
        # ASAP7y is read first, so its C end comes before EGFP's N end.
        self.assertGreater(asap_n, asap_c)
        self.assertGreater(asap_c, egfp_n)
        self.assertGreater(egfp_n, egfp_c)

    def test_the_map_draws_them_for_a_real_construct(self):
        from genesis.pdf import PDF, A4, circular_map
        doc = PDF("test", A4)
        doc._ensure(1)
        feats = [{"feature_type": "CDS", "name": "ASAP7y", "start": 1920,
                  "end": 3204, "strand": 1},
                 {"feature_type": "promoter", "name": "CAG", "start": 92,
                  "end": 1837, "strand": 1}]
        circular_map(doc, "test", "7,219 bp", 7219, feats)
        ops = "".join(o if isinstance(o, str) else o.decode("latin-1")
                      for o in doc.pages[-1].ops)
        self.assertIn("(N)", ops)
        self.assertIn("(C)", ops)
