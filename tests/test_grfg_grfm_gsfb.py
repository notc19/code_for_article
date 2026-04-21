"""
pytest for:
  - get_regions_for_mutations (get_regions_for_muts.py)
  - get_regions_for_genes     (get_regions_for_genes.py)
  - get_sequences_from_bed    (get_sequences_from_bed.py)
"""

from __future__ import annotations

import logging
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

@pytest.fixture
def logger():
    log = logging.getLogger("test")
    log.setLevel(logging.DEBUG)
    return log


def _make_annot_for_muts_plus():
    """Single-exon gene on + strand, chr1, start_codon at 1000."""
    return {
        "ENST00000000001.1": [
            "chr1",
            [["900", "1200"]],
            "1000",
            "+",
        ]
    }

def _make_annot_for_muts_minus():
    """Single-exon gene on - strand, chr2, start_codon at 2000."""
    return {
        "ENST00000000002.1": [
            "chr2",
            [["1800", "2100"]],
            "2000",
            "-",
        ]
    }

def _make_annot_for_muts_multi_exon():
    """Two-exon gene on + strand, chr3."""
    return {
        "ENST00000000003.2": [
            "chr3",
            [["100", "200"], ["300", "400"]],
            "100",
            "+",
        ]
    }

def _make_annot_for_genes():
    return {
        "ENSG00000000001": [
            "chr1",                          # [0] chrom
            ["1000", "5000"],                # [1] gene_coords
            [[1100, 1300], [2000, 2200]],    # [2] exon_list (ints)
            [["1050", "1099"]],              # [3] utr_list
            ["ENSE00000000001", "ENSE00000000002"],  # [4] exon_id_list
            ["ENST00000000001.1"],           # [5] mane_select
            [],                              # [6] mane_plus_clinical
            ["ENST00000000001.1", "ENST00000000001.2"],  # [7] transcript_list
            "ENST00000000001.1",             # [8] selected_tid
        ],
        "ENSG00000000002": [
            "chr2",
            ["500", "3000"],
            [[600, 800], [1000, 1200], [1500, 1700]],
            [],
            ["ENSE00000000003", "ENSE00000000004", "ENSE00000000005"],
            [],
            [],
            ["ENST00000000002.1"],
            "ENST00000000002.1",
        ],
    }

def _make_gene_name_dict():
    return {
        "GENE1": "ENSG00000000001",
        "GENE2": "ENSG00000000002",
    }

# 1.  get_regions_for_mutations
@pytest.fixture(autouse=True, scope="session")
def patch_constants():
    """Provide a minimal constants module so imports don't fail."""
    import types
    fake_constants = types.ModuleType("constants")
    fake_constants.genome_versions_chromosome_dict = {
        "GRCh38": {
            "NC_000001.11": "chr1",
            "NC_000002.12": "chr2",
            "NC_000023.11": "chrX",
        },
        "GRCh37": {},
    }
    fake_constants.rev_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    sys.modules["constants"] = fake_constants
    yield

@pytest.fixture(scope="session")
def grfm(patch_constants):
    """Import get_regions_for_mutations module."""
    import importlib
    import types
    mafc = types.ModuleType("merge_and_filter_coordinates")
    sys.modules["merge_and_filter_coordinates"] = mafc
    src_dir = Path(__file__).parent.parent
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))

    import importlib.util, types
    spec = importlib.util.spec_from_file_location(
        "get_regions_for_mutations", src_dir / "get_regions_for_mutations.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["get_regions_for_mutations"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def grfg(grfm):
    """Import get_regions_for_genes module (depends on grfm being registered)."""
    import importlib.util
    src_dir = Path(__file__).parent.parent
    spec = importlib.util.spec_from_file_location(
        "get_regions_for_genes", src_dir / "get_regions_for_genes.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["get_regions_for_genes"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def gsfb(patch_constants):
    """Import get_sequences_from_bed module."""
    import importlib.util
    src_dir = Path(__file__).parent.parent
    spec = importlib.util.spec_from_file_location(
        "get_sequences_from_bed", src_dir / "get_sequences_from_bed.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["get_sequences_from_bed"] = mod
    spec.loader.exec_module(mod)
    return mod

def parse_bed(line: str) -> dict:
    parts = line.split("\t")
    return {
        "chrom": parts[0],
        "start": int(parts[1]),
        "end": int(parts[2]),
        "ref": parts[3],
        "name": parts[4],
        "type": parts[5],
        "orig_chrom": parts[6] if len(parts) > 6 else None,
    }

class TestGetCoordsForVcf:
    def test_basic_vcf(self, grfm, logger):
        vcf = ["chr1", "rs123", "1000", "A", "G"]
        result = grfm.get_coords_for_vcf(vcf, flank_length=10)
        b = parse_bed(result)
        assert b["chrom"] == "chr1"
        assert b["start"] == 990
        assert b["end"] == 1010
        assert b["ref"] == "A"
        assert b["name"] == "rs123"
        assert b["type"] == "snp"

    def test_flank_zero(self, grfm, logger):
        vcf = ["chrX", "rs999", "500", "C", "T"]
        result = grfm.get_coords_for_vcf(vcf, flank_length=0)
        b = parse_bed(result)
        assert b["start"] == 500
        assert b["end"] == 500

    def test_large_flank(self, grfm, logger):
        vcf = ["chr2", "rs1", "5000", "G", "A"]
        result = grfm.get_coords_for_vcf(vcf, flank_length=200)
        b = parse_bed(result)
        assert b["start"] == 4800
        assert b["end"] == 5200


class TestGetCoordsForHgvsG:
    def test_basic_hgvsg(self, grfm, logger):
        hgvsg = ["chr1", "g.1000A>T"]
        result = grfm.get_coords_for_hgvs_g(hgvsg, flank_length=5)
        b = parse_bed(result)
        assert b["chrom"] == "chr1"
        assert b["start"] == 995
        assert b["end"] == 1005
        assert b["ref"] == "A"
        assert b["type"] == "snp"

    def test_hgvsg_with_version_prefix(self, grfm, logger):
        # NC_000001.11:g.12345C>G
        hgvsg = ["NC_000001.11", "g.12345C>G"]
        result = grfm.get_coords_for_hgvs_g(hgvsg, flank_length=10)
        b = parse_bed(result)
        assert b["start"] == 12335
        assert b["end"] == 12355
        assert b["ref"] == "C"

    def test_hgvsg_name_field_contains_full_mutation(self, grfm, logger):
        hgvsg = ["chr3", "g.9999T>A"]
        result = grfm.get_coords_for_hgvs_g(hgvsg, flank_length=1)
        b = parse_bed(result)
        assert "chr3" in b["name"]
        assert "g.9999T>A" in b["name"]

class TestGetCoordsForHgvsC:
    def test_exact_transcript_match_plus_strand(self, grfm, logger):
        """Versioned transcript ID present in annot_for_muts → exact match."""
        afm = _make_annot_for_muts_plus()
        hgvsc = ["ENST00000000001.1", "c.1A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 5, afm, logger)
        assert result is not None
        b = parse_bed(result)
        assert b["chrom"] == "chr1"
        assert b["type"] == "snp"
        assert b["start"] < b["end"]

    def test_unversioned_transcript_match(self, grfm, logger):
        """Unversioned input ID resolved via unversioned_index."""
        afm = _make_annot_for_muts_plus()
        hgvsc = ["ENST00000000001", "c.1A>T"]   # no version suffix
        result = grfm.get_coords_for_hgvs_c(hgvsc, 5, afm, logger)
        assert result is not None

    def test_unknown_transcript_returns_none(self, grfm, logger):
        """Transcript not in annotation → returns None, emits warning."""
        afm = _make_annot_for_muts_plus()
        hgvsc = ["ENST99999999999.1", "c.1A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 5, afm, logger)
        assert result is None

    def test_mane_fallback(self, grfm, logger):
        """
        Input transcript not in annot_for_muts but gene is known →
        MANE_Select transcript is used as fallback.
        """
        afm = _make_annot_for_muts_plus()
        afg = _make_annot_for_genes()
        # ENST00000000001.2 is in transcript_list but NOT in afm
        hgvsc = ["ENST00000000001.2", "c.1A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 5, afm, logger, annot_for_genes=afg)
        # MANE_Select for ENSG00000000001 is ENST00000000001.1 which IS in afm
        assert result is not None

    def test_intronic_plus_offset(self, grfm, logger):
        """c.1+5A>T → intronic position downstream of exon boundary."""
        afm = _make_annot_for_muts_plus()
        hgvsc = ["ENST00000000001.1", "c.1+5A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 0, afm, logger)
        assert result is not None
        b = parse_bed(result)
        assert b["type"] == "snp"

    def test_intronic_minus_offset(self, grfm, logger):
        """c.1-5A>T → intronic position upstream of exon boundary."""
        afm = _make_annot_for_muts_plus()
        hgvsc = ["ENST00000000001.1", "c.1-5A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 0, afm, logger)
        assert result is not None

    def test_utr5_negative_position(self, grfm, logger):
        """c.-10A>T → 5' UTR position."""
        afm = _make_annot_for_muts_plus()
        hgvsc = ["ENST00000000001.1", "c.-10A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 0, afm, logger)
        assert result is not None

    def test_utr3_star_position(self, grfm, logger):
        """c.*5A>T → 3' UTR position."""
        afm = _make_annot_for_muts_plus()
        hgvsc = ["ENST00000000001.1", "c.*5A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 0, afm, logger)
        assert result is not None

    def test_minus_strand_transcript(self, grfm, logger):
        """Coordinate calculation on - strand transcript."""
        afm = _make_annot_for_muts_minus()
        hgvsc = ["ENST00000000002.1", "c.1A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 5, afm, logger)
        assert result is not None
        b = parse_bed(result)
        assert b["chrom"] == "chr2"

    def test_bed_name_contains_gene_id_when_known(self, grfm, logger):
        """BED name field should be GENE_ID|ENST_unversioned:c.notation."""
        afm = _make_annot_for_muts_plus()
        afg = _make_annot_for_genes()
        hgvsc = ["ENST00000000001.1", "c.1A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 5, afm, logger, annot_for_genes=afg)
        assert result is not None
        b = parse_bed(result)
        assert "ENSG00000000001" in b["name"]
        assert "ENST00000000001" in b["name"]
        assert "c.1A>T" in b["name"]

    def test_pos_c_exceeds_transcript_length_returns_none(self, grfm, logger):
        """pos_c larger than transcript coding length → ValueError caught → None."""
        afm = {
            "ENST00000000001.1": [
                "chr1",
                [["1000", "1010"]],   #11 bp exon
                "1000",
                "+",
            ]
        }
        hgvsc = ["ENST00000000001.1", "c.9999A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 5, afm, logger)
        assert result is None

    def test_unsupported_notation_returns_none(self, grfm, logger):
        """Notation that doesn't start with digit, -, or * → None."""
        afm = _make_annot_for_muts_plus()
        hgvsc = ["ENST00000000001.1", "c.?100A>T"]
        result = grfm.get_coords_for_hgvs_c(hgvsc, 5, afm, logger)
        assert result is None

class TestFixChromosomeName:
    def test_already_chr_prefix(self, grfm, logger):
        coords = ["chr1\t100\t200\tA\trs1\tsnp"]
        result = grfm.fix_chromosome_name_and_check_coords(coords, "GRCh38", logger)
        assert len(result) == 1
        parts = result[0].split("\t")
        assert parts[0] == "chr1"
        assert parts[6] == "chr1"   # original chrom preserved in col 7

    def test_numeric_chrom_gets_chr_prefix(self, grfm, logger):
        coords = ["1\t100\t200\tA\trs1\tsnp"]
        result = grfm.fix_chromosome_name_and_check_coords(coords, "GRCh38", logger)
        assert len(result) == 1
        assert result[0].split("\t")[0] == "chr1"
        assert result[0].split("\t")[6] == "1"   # original preserved

    def test_X_Y_get_chr_prefix(self, grfm, logger):
        coords = ["X\t100\t200\tA\trs1\tsnp", "Y\t100\t200\tA\trs2\tsnp"]
        result = grfm.fix_chromosome_name_and_check_coords(coords, "GRCh38", logger)
        chroms = [r.split("\t")[0] for r in result]
        assert "chrX" in chroms
        assert "chrY" in chroms

    def test_nc_accession_mapped(self, grfm, logger):
        coords = ["NC_000001.11\t100\t200\tA\trs1\tsnp"]
        result = grfm.fix_chromosome_name_and_check_coords(coords, "GRCh38", logger)
        assert len(result) == 1
        assert result[0].split("\t")[0] == "chr1"

    def test_invalid_chrom_skipped(self, grfm, logger):
        coords = ["INVALID\t100\t200\tA\trs1\tsnp"]
        result = grfm.fix_chromosome_name_and_check_coords(coords, "GRCh38", logger)
        assert len(result) == 0

    def test_negative_coords_skipped(self, grfm, logger):
        coords = ["chr1\t-10\t200\tA\trs1\tsnp"]
        result = grfm.fix_chromosome_name_and_check_coords(coords, "GRCh38", logger)
        assert len(result) == 0

    def test_seven_columns_output(self, grfm, logger):
        """Output must always have 7 tab-separated columns."""
        coords = ["chr1\t100\t200\tA\trs1\tsnp"]
        result = grfm.fix_chromosome_name_and_check_coords(coords, "GRCh38", logger)
        assert len(result[0].split("\t")) == 7

    def test_missing_ref_allele_skipped(self, grfm, logger):
        """Lines with fewer than 6 columns (no ref/type) are skipped with warning."""
        coords = ["chr1\t100\t200"]
        result = grfm.fix_chromosome_name_and_check_coords(coords, "GRCh38", logger)
        assert len(result) == 0

class TestGetRegionsForMuts:
    def test_vcf_mutation_processed(self, grfm, logger, tmp_path):
        afm = _make_annot_for_muts_plus()
        muts = ["chr1\trs1\t1000\tA\tG"]
        result = grfm.get_regions_for_muts(
            muts, afm, 10, False, tmp_path, logger, "GRCh38"
        )
        assert len(result) == 1
        b = parse_bed(result[0])
        assert b["chrom"] == "chr1"
        assert b["type"] == "snp"

    def test_hgvsg_mutation_processed(self, grfm, logger, tmp_path):
        afm = _make_annot_for_muts_plus()
        muts = ["chr1:g.1000A>T"]
        result = grfm.get_regions_for_muts(
            muts, afm, 5, False, tmp_path, logger, "GRCh38"
        )
        assert len(result) == 1

    def test_hgvsc_mutation_processed(self, grfm, logger, tmp_path):
        afm = _make_annot_for_muts_plus()
        muts = ["ENST00000000001.1:c.1A>T"]
        result = grfm.get_regions_for_muts(
            muts, afm, 5, False, tmp_path, logger, "GRCh38"
        )
        assert len(result) == 1

    def test_comment_lines_skipped(self, grfm, logger, tmp_path):
        afm = _make_annot_for_muts_plus()
        muts = ["# this is a comment", "chr1\trs1\t1000\tA\tG"]
        result = grfm.get_regions_for_muts(
            muts, afm, 5, False, tmp_path, logger, "GRCh38"
        )
        assert len(result) == 1

    def test_mixed_formats(self, grfm, logger, tmp_path):
        afm = _make_annot_for_muts_plus()
        muts = [
            "chr1\trs1\t1000\tA\tG",
            "chr1:g.2000C>T",
            "ENST00000000001.1:c.1A>T",
        ]
        result = grfm.get_regions_for_muts(
            muts, afm, 5, False, tmp_path, logger, "GRCh38"
        )
        assert len(result) == 3

    def test_write_intermediate_creates_file(self, grfm, logger, tmp_path):
        afm = _make_annot_for_muts_plus()
        muts = ["chr1\trs1\t1000\tA\tG"]
        grfm.get_regions_for_muts(
            muts, afm, 5, True, tmp_path, logger, "GRCh38"
        )
        assert (tmp_path / "INPUT_coords_for_mutations.bed").exists()

    def test_numeric_chrom_normalised_in_output(self, grfm, logger, tmp_path):
        afm = _make_annot_for_muts_plus()
        muts = ["1\trs1\t1000\tA\tG"]
        result = grfm.get_regions_for_muts(
            muts, afm, 5, False, tmp_path, logger, "GRCh38"
        )
        assert len(result) == 1
        assert result[0].split("\t")[0] == "chr1"

    def test_hgvsc_with_annot_for_genes_mane_fallback(self, grfm, logger, tmp_path):
        """MANE fallback: unrecognised transcript version resolved via gene."""
        afm = _make_annot_for_muts_plus()
        afg = _make_annot_for_genes()
        muts = ["ENST00000000001.2:c.1A>T"]   # .2 not in afm, but .1 is (MANE)
        result = grfm.get_regions_for_muts(
            muts, afm, 5, False, tmp_path, logger, "GRCh38",
            annot_for_genes=afg,
        )
        assert len(result) == 1

    def test_n_notation_processed(self, grfm, logger, tmp_path):
        """n. notation (non-coding) handled same as c."""
        afm = _make_annot_for_muts_plus()
        muts = ["ENST00000000001.1:n.1A>T"]
        result = grfm.get_regions_for_muts(
            muts, afm, 5, False, tmp_path, logger, "GRCh38"
        )
        assert len(result) == 1


class TestIndexHelpers:
    def test_build_unversioned_index(self, grfm):
        afm = _make_annot_for_muts_plus()
        idx = grfm.build_unversioned_muts_index(afm)
        assert "ENST00000000001" in idx
        assert idx["ENST00000000001"] == "ENST00000000001.1"

    def test_get_gene_id_for_transcript(self, grfm):
        afg = _make_annot_for_genes()
        tid_to_gene = grfm.get_gene_id_for_transcript(afg)
        assert tid_to_gene["ENST00000000001"] == "ENSG00000000001"
        assert tid_to_gene["ENST00000000002"] == "ENSG00000000002"

    def test_resolve_transcript_exact(self, grfm, logger):
        afm = _make_annot_for_muts_plus()
        idx = grfm.build_unversioned_muts_index(afm)
        afg = _make_annot_for_genes()
        tid_to_gene = grfm.get_gene_id_for_transcript(afg)
        tid, gene_id, note = grfm._resolve_transcript(
            "ENST00000000001.1", afm, idx, tid_to_gene, afg, logger
        )
        assert tid == "ENST00000000001.1"
        assert note == "exact"

    def test_resolve_transcript_unversioned(self, grfm, logger):
        afm = _make_annot_for_muts_plus()
        idx = grfm.build_unversioned_muts_index(afm)
        afg = _make_annot_for_genes()
        tid_to_gene = grfm.get_gene_id_for_transcript(afg)
        tid, gene_id, note = grfm._resolve_transcript(
            "ENST00000000001", afm, idx, tid_to_gene, afg, logger
        )
        assert tid == "ENST00000000001.1"
        assert note == "unversioned"

    def test_resolve_transcript_mane_fallback(self, grfm, logger):
        afm = _make_annot_for_muts_plus()
        idx = grfm.build_unversioned_muts_index(afm)
        afg = _make_annot_for_genes()
        tid_to_gene = grfm.get_gene_id_for_transcript(afg)
        # .2 is in transcript_list but not in afm; MANE_Select .1 is in afm
        tid, gene_id, note = grfm._resolve_transcript(
            "ENST00000000001.2", afm, idx, tid_to_gene, afg, logger
        )
        assert tid == "ENST00000000001.1"

    def test_resolve_transcript_not_found(self, grfm, logger):
        afm = _make_annot_for_muts_plus()
        idx = grfm.build_unversioned_muts_index(afm)
        tid, gene_id, note = grfm._resolve_transcript(
            "ENST99999999999.1", afm, idx, {}, {}, logger
        )
        assert tid is None
        assert note is None


# 2.  get_regions_for_genes

class TestGetRegionsForGenes:

    def _run(self, grfg, grfm, genes, gene_regions, intron_depth, length, tmp_path, logger,
             genome_version="GRCh38"):
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        return grfg.get_regions_for_genes(
            genes, afg, gnd, gene_regions, intron_depth,
            False, tmp_path, logger, genome_version, length
        )

    def test_full_region(self, grfg, grfm, logger, tmp_path):
        result = self._run(grfg, grfm, ["GENE1"], "full", 0, 100, tmp_path, logger)
        assert len(result) == 1
        b = parse_bed(result[0])
        assert b["chrom"] == "chr1"
        # gene_coords = [1000, 5000], depth=0 → start=1000, end=5000
        assert b["start"] == 1000
        assert b["end"] == 5000
        assert b["type"] == "full"

    def test_full_region_with_intron_depth(self, grfg, grfm, logger, tmp_path):
        result = self._run(grfg, grfm, ["GENE1"], "full", 50, 100, tmp_path, logger)
        b = parse_bed(result[0])
        assert b["start"] == 950
        assert b["end"] == 5050

    def test_exon_region_count(self, grfg, grfm, logger, tmp_path):
        """GENE1 has 2 exons → 2 BED lines."""
        result = self._run(grfg, grfm, ["GENE1"], "exon", 0, 0, tmp_path, logger)
        assert len(result) == 2

    def test_exon_region_coords(self, grfg, grfm, logger, tmp_path):
        result = self._run(grfg, grfm, ["GENE1"], "exon", 0, 0, tmp_path, logger)
        b0 = parse_bed(result[0])
        # exon_list[0] = [1100, 1300], depth=0
        assert b0["start"] == 1100
        assert b0["end"] == 1300

    def test_exon_region_with_intron_depth(self, grfg, grfm, logger, tmp_path):
        result = self._run(grfg, grfm, ["GENE1"], "exon", 10, 0, tmp_path, logger)
        b0 = parse_bed(result[0])
        assert b0["start"] == 1090
        assert b0["end"] == 1310

    def test_exon_minimum_length_padding(self, grfg, grfm, logger, tmp_path):
        """When exon_length < length, region is padded to at least `length` bp."""
        # exon[0] = [1100,1300] → length=200 bp; request length=500 → should pad
        result = self._run(grfg, grfm, ["GENE1"], "exon", 0, 500, tmp_path, logger)
        b0 = parse_bed(result[0])
        assert (b0["end"] - b0["start"]) >= 500

    def test_specific_exon_by_number(self, grfg, grfm, logger, tmp_path):
        """'1' selects only the first exon."""
        result = self._run(grfg, grfm, ["GENE1"], "1", 0, 0, tmp_path, logger)
        assert len(result) == 1
        b = parse_bed(result[0])
        assert b["start"] == 1100
        assert b["end"] == 1300

    def test_exon_range(self, grfg, grfm, logger, tmp_path):
        """'1-2' selects exons 1 and 2 of GENE2 (3 exons total)."""
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        result = grfg.get_regions_for_genes(
            ["GENE2"], afg, gnd, "1-2", 0, False, tmp_path, logger, "GRCh38", 100
        )
        assert len(result) == 2

    def test_exon_by_id(self, grfg, grfm, logger, tmp_path):
        """Exon ID string selects the matching exon."""
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        result = grfg.get_regions_for_genes(
            ["GENE1"], afg, gnd, "exon", 0, False, tmp_path, logger, "GRCh38", 100
        )
        assert len(result) == 2
        b = parse_bed(result[0])
        assert b["name"] == "ENSE00000000001"

    def test_unknown_gene_skipped(self, grfg, grfm, logger, tmp_path):
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        result = grfg.get_regions_for_genes(
            ["UNKNOWN_GENE"], afg, gnd, "exon", 0, False, tmp_path, logger, "GRCh38", 100
        )
        assert len(result) == 0

    def test_gene_by_ensembl_id(self, grfg, grfm, logger, tmp_path):
        """Gene can be specified by Ensembl ID directly."""
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        result = grfg.get_regions_for_genes(
            ["ENSG00000000001"], afg, gnd, "exon", 0, False, tmp_path, logger, "GRCh38", 100
        )
        assert len(result) == 2

    def test_per_gene_region_override(self, grfg, grfm, logger, tmp_path):
        """Gene line 'GENE1 full' overrides the default gene_regions='exon'."""
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        result = grfg.get_regions_for_genes(
            ["GENE1 full"], afg, gnd, "exon", 0, False, tmp_path, logger, "GRCh38", 100
        )
        assert len(result) == 1
        b = parse_bed(result[0])
        assert b["type"] == "full"

    def test_per_gene_depth_override(self, grfg, grfm, logger, tmp_path):
        """Gene line 'GENE1 exon 20' overrides intron_depth."""
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        result = grfg.get_regions_for_genes(
            ["GENE1 exon 20"], afg, gnd, "exon", 0, False, tmp_path, logger, "GRCh38", 100
        )
        b0 = parse_bed(result[0])
        assert b0["start"] == 1080   # 1100 - 20
        assert b0["end"] == 1320     # 1300 + 20

    def test_write_intermediate_creates_file(self, grfg, grfm, logger, tmp_path):
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        grfg.get_regions_for_genes(
            ["GENE1"], afg, gnd, "exon", 0, True, tmp_path, logger, "GRCh38", 100
        )
        assert (tmp_path / "INPUT_coords_for_genes.bed").exists()

    def test_output_has_seven_columns(self, grfg, grfm, logger, tmp_path):
        """fix_chromosome_name_and_check_coords appends original chrom → 7 cols."""
        result = self._run(grfg, grfm, ["GENE1"], "exon", 0, 0, tmp_path, logger)
        for line in result:
            assert len(line.split("\t")) == 7, f"Expected 7 cols: {line}"

    def test_selected_tid_in_full_region_name(self, grfg, grfm, logger, tmp_path):
        """For 'full' region, col 5 (name) should be the selected_tid."""
        result = self._run(grfg, grfm, ["GENE1"], "full", 0, 0, tmp_path, logger)
        b = parse_bed(result[0])
        assert b["name"] == "ENST00000000001.1"

    def test_exon_id_in_exon_region_name(self, grfg, grfm, logger, tmp_path):
        """For 'exon' region, col 5 (name) should be the exon stable ID."""
        result = self._run(grfg, grfm, ["GENE1"], "exon", 0, 0, tmp_path, logger)
        names = [parse_bed(r)["name"] for r in result]
        assert "ENSE00000000001" in names
        assert "ENSE00000000002" in names

    def test_comma_separated_exon_list(self, grfg, grfm, logger, tmp_path):
        """'1,2' selects exons 1 and 2."""
        afg = _make_annot_for_genes()
        gnd = _make_gene_name_dict()
        result = grfg.get_regions_for_genes(
            ["GENE2 1,2"], afg, gnd, "exon", 0, False, tmp_path, logger, "GRCh38", 100
        )
        assert len(result) == 2


# 3.  get_sequences_from_bed
FAKE_CHR1_SEQ = "A" * 500 + "C" * 500   # 1000 bp total
FAKE_CHR2_SEQ = "G" * 300 + "T" * 300   # 600 bp total


@pytest.fixture
def fake_genome(tmp_path):
    """Write a minimal two-chromosome FASTA to a temp file."""
    fa = tmp_path / "genome.fa"
    fa.write_text(
        f">chr1 Homo sapiens chromosome 1\n{FAKE_CHR1_SEQ}\n"
        f">chr2 Homo sapiens chromosome 2\n{FAKE_CHR2_SEQ}\n"
    )
    return fa


@pytest.fixture
def fake_genome_numeric_chrom(tmp_path):
    """FASTA with numeric chromosome names (Ensembl style)."""
    fa = tmp_path / "genome_ensembl.fa"
    fa.write_text(
        f">1 Homo sapiens chromosome 1\n{FAKE_CHR1_SEQ}\n"
    )
    return fa


class TestCreateFastaSequences:
    """Unit tests for create_fasta_sequences (no file I/O)."""

    def test_region_record_7col(self, gsfb, logger):
        """7-column BED line (exon region) → SeqRecord with correct sequence."""
        # chr1, start=1, end=10 → first 10 'A's
        coords_list = ["chr1\t1\t10\tGENE1\tENSE001\texon\tchr1"]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert len(records) == 1
        assert str(records[0].seq) == "A" * 10
        assert records[0].id == "ENSE001"
        assert records[0].name == "GENE1"
        assert "Region" in records[0].description

    def test_snp_record_7col(self, gsfb, logger):
        """7-column SNP BED line → SeqRecord with 'snp' name."""
        # chr1, start=496, end=506 → 5 A's + 1 boundary + 5 C's (flank=5)
        # ref allele = 'A', middle position = 500 (0-based index 4 in 10-bp window)
        # sequence[5] = 'A' (1-based pos 501 is 'C', but let's use a safe position)
        # Use start=1, end=11, ref='A', middle index=5 → seq[5]='A'
        coords_list = ["chr1\t1\t11\tA\trs123\tsnp\tchr1"]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert len(records) == 1
        assert records[0].name == "snp"
        assert records[0].id == "rs123"
        assert "SNP" in records[0].description

    def test_snp_ref_allele_mismatch_skipped(self, gsfb, logger):
        """SNP with wrong ref allele → record skipped, warning emitted."""
        # Position 501-510 is all 'C', but ref='A' → mismatch
        coords_list = ["chr1\t501\t511\tA\trs_bad\tsnp\tchr1"]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert len(records) == 0

    def test_snp_ref_allele_match_kept(self, gsfb, logger):
        """SNP with correct ref allele → record kept."""
        # Position 1-11 is all 'A', ref='A' → match
        coords_list = ["chr1\t1\t11\tA\trs_good\tsnp\tchr1"]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert len(records) == 1

    def test_del_allele_not_validated(self, gsfb, logger):
        """del/ins alleles bypass ref-allele check → always kept."""
        coords_list = ["chr1\t501\t511\tdelATG\trs_del\tsnp\tchr1"]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert len(records) == 1

    def test_6col_fallback(self, gsfb, logger):
        """6-column BED line (no original_chrom) still produces a record."""
        coords_list = ["chr1\t1\t10\tGENE1\tENSE001\texon"]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert len(records) == 1
        assert str(records[0].seq) == "A" * 10

    def test_short_coords_fallback_description(self, gsfb, logger):
        """Fewer than 6 columns → description is the raw coords list."""
        coords_list = ["chr1\t1\t10"]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert len(records) == 1
        assert records[0].id == "ID"
        assert records[0].name == "NAME"

    def test_multiple_records(self, gsfb, logger):
        coords_list = [
            "chr1\t1\t5\tGENE1\tENSE001\texon\tchr1",
            "chr1\t6\t10\tGENE1\tENSE002\texon\tchr1",
        ]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert len(records) == 2
        assert str(records[0].seq) == "A" * 5
        assert str(records[1].seq) == "A" * 5

    def test_sequence_slice_correctness(self, gsfb, logger):
        """Verify 1-based inclusive slicing: start=498, end=503 spans A→C boundary."""
        # FAKE_CHR1_SEQ: positions 1-500 = 'A', 501-1000 = 'C'
        # start=498, end=503 → seq = 'AAA' + 'CCC' (3 A's + 3 C's)
        coords_list = ["chr1\t498\t503\tGENE1\tENSE001\texon\tchr1"]
        records = gsfb.create_fasta_sequences(coords_list, FAKE_CHR1_SEQ, logger)
        assert str(records[0].seq) == "AAACCC"


class TestGetFastaSeqsFromBed:
    """Integration tests for get_fasta_seqs_from_bed (reads actual FASTA file)."""

    def test_basic_extraction(self, gsfb, fake_genome, logger, tmp_path):
        coord_list = ["chr1\t1\t10\tGENE1\tENSE001\texon\tchr1"]
        records = gsfb.get_fasta_seqs_from_bed(
            coord_list, fake_genome, "GRCh38", False, tmp_path, logger
        )
        assert len(records) == 1
        assert str(records[0].seq) == "A" * 10

    def test_last_chromosome_flushed(self, gsfb, fake_genome, logger, tmp_path):
        """chr2 is the last chrom in FASTA → must be flushed at EOF."""
        coord_list = ["chr2\t1\t10\tGENE2\tENSE002\texon\tchr2"]
        records = gsfb.get_fasta_seqs_from_bed(
            coord_list, fake_genome, "GRCh38", False, tmp_path, logger
        )
        assert len(records) == 1
        assert str(records[0].seq) == "G" * 10

    def test_two_chromosomes(self, gsfb, fake_genome, logger, tmp_path):
        coord_list = [
            "chr1\t1\t5\tGENE1\tENSE001\texon\tchr1",
            "chr2\t1\t5\tGENE2\tENSE002\texon\tchr2",
        ]
        records = gsfb.get_fasta_seqs_from_bed(
            coord_list, fake_genome, "GRCh38", False, tmp_path, logger
        )
        assert len(records) == 2
        seqs = {r.id: str(r.seq) for r in records}
        assert seqs["ENSE001"] == "A" * 5
        assert seqs["ENSE002"] == "G" * 5

    def test_numeric_chrom_in_fasta_normalised(self, gsfb, fake_genome_numeric_chrom, logger, tmp_path):
        """FASTA with '1' header → normalised to 'chr1' for GRCh38."""
        coord_list = ["chr1\t1\t5\tGENE1\tENSE001\texon\tchr1"]
        records = gsfb.get_fasta_seqs_from_bed(
            coord_list, fake_genome_numeric_chrom, "GRCh38", False, tmp_path, logger
        )
        assert len(records) == 1
        assert str(records[0].seq) == "A" * 5

    def test_write_intermediate_creates_fasta(self, gsfb, fake_genome, logger, tmp_path):
        coord_list = ["chr1\t1\t10\tGENE1\tENSE001\texon\tchr1"]
        gsfb.get_fasta_seqs_from_bed(
            coord_list, fake_genome, "GRCh38", True, tmp_path, logger
        )
        assert (tmp_path / "INPUT_regions_sequences.fasta").exists()

    def test_chrom_not_in_fasta_returns_empty(self, gsfb, fake_genome, logger, tmp_path):
        """Requesting a chromosome absent from FASTA → empty result."""
        coord_list = ["chr99\t1\t10\tGENE1\tENSE001\texon\tchr99"]
        records = gsfb.get_fasta_seqs_from_bed(
            coord_list, fake_genome, "GRCh38", False, tmp_path, logger
        )
        assert len(records) == 0

    def test_snp_ref_mismatch_skipped_in_integration(self, gsfb, fake_genome, logger, tmp_path):
        """SNP with wrong ref allele is skipped end-to-end."""
        # chr1 pos 501-511 is all 'C'; ref='A' → mismatch → skip
        coord_list = ["chr1\t501\t511\tA\trs_bad\tsnp\tchr1"]
        records = gsfb.get_fasta_seqs_from_bed(
            coord_list, fake_genome, "GRCh38", False, tmp_path, logger
        )
        assert len(records) == 0

    def test_multiple_regions_same_chrom(self, gsfb, fake_genome, logger, tmp_path):
        coord_list = [
            "chr1\t1\t5\tGENE1\tENSE001\texon\tchr1",
            "chr1\t496\t500\tGENE1\tENSE002\texon\tchr1",
            "chr1\t501\t505\tGENE1\tENSE003\texon\tchr1",
        ]
        records = gsfb.get_fasta_seqs_from_bed(
            coord_list, fake_genome, "GRCh38", False, tmp_path, logger
        )
        assert len(records) == 3
        seqs = {r.id: str(r.seq) for r in records}
        assert seqs["ENSE001"] == "AAAAA"
        assert seqs["ENSE002"] == "AAAAA"
        assert seqs["ENSE003"] == "CCCCC"
