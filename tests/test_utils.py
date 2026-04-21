"""
Unit tests for utils.py — GTF annotation parser for TS capture probe design.

Coverage
--------
Helpers
  _parse_attributes       – quoted values, multi-value tag keys, empty string
  _get_attr               – key present / absent / fallback chain / default
  _has_tag                – tag present / absent / no tag key at all
  _detect_source          – Ensembl (bare chrom), GENCODE (chr-prefix), mixed,
                            empty annotation
  _biotype_key            – all four combinations of source × level
  _transcript_id_versioned – Ensembl (bare + version attr), Ensembl (already
                             dotted), GENCODE (version embedded), missing version
  _is_utr_feature         – all three UTR feature names + non-UTR names

prepare_annot_for_muts
  – returns None for genomic-only mutations (no c./n.)
  – returns None for empty mutation list
  – parses Ensembl GTF: versioned ID assembled, exons + start_codon collected
  – parses GENCODE GTF: versioned ID from embedded version, chr-prefix preserved
  – transcript without start_codon is excluded
  – multiple transcripts in one annotation
  – n. notation triggers parsing
  – transcript with exons but no start_codon → excluded
  – malformed / short lines are skipped gracefully

prepare_annot_for_genes
  – Ensembl GTF: gene_id stripped of version, gene_name_dict populated
  – GENCODE GTF: gene_id stripped of version, chr-prefix on chrom
  – protein_coding exons collected; non-protein_coding exons ignored
  – exon deduplication (same exon_id appearing in multiple transcripts)
  – UTR deduplication (same coords from multiple transcripts)
  – Ensembl UTR feature names (five_prime_utr / three_prime_utr)
  – GENCODE UTR feature name (UTR)
  – stop_codon trimming on + strand
  – stop_codon trimming on − strand
  – MANE_Select tag collected
  – MANE_Plus_Clinical tag collected
  – transcript with both MANE tags simultaneously
  – transcript with neither MANE tag
  – multiple MANE_Select transcripts per gene (edge case)
  – gene with no protein-coding exons → empty exon_list
  – gene with no UTRs → empty utr_list
  – gene with no stop_codon → no trimming applied
  – two consecutive genes parsed correctly (index-jump logic)
  – malformed / short lines are skipped gracefully
"""

import sys
import os
import logging
from datetime import datetime
import pytest

# Allow running from repo root or tests/ directory
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def create_logger(timestamp: datetime):
    """
    Create the main-process logger.
    Returns the Logger object; its handlers are also passed to the listener thread.
    """
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(fmt)

    file_handler = logging.FileHandler(
        f'log_oligodesigner_{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}-{timestamp.second}.txt',
        mode='w'
    )
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(fmt)

    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.handlers = []
    root.addHandler(stream_handler)
    root.addHandler(file_handler)

    return root, [stream_handler, file_handler]

timestamp = datetime.now()
logger, log_handlers = create_logger(timestamp)


from utils import (
    _parse_attributes,
    _get_attr,
    _has_tag,
    _detect_source,
    _biotype_key,
    _transcript_id_versioned,
    _is_utr_feature,
    prepare_annot_for_muts,
    prepare_annot_for_genes,
)


# ===========================================================================
# GTF line factories
# ===========================================================================
# These helpers build realistic, tab-separated GTF lines so that test fixtures
# read like documentation rather than raw string soup.

def _gtf_line(chrom, source, feature, start, end, strand, attrs: str) -> str:
    return f"{chrom}\t{source}\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"


def _attrs(**kwargs) -> str:
    """Build a GTF attribute string from keyword arguments.
    Pass tag as a list to emit multiple  tag "...";  entries."""
    parts = []
    for k, v in kwargs.items():
        if k == "tag" and isinstance(v, list):
            for t in v:
                parts.append(f'tag "{t}";')
        else:
            parts.append(f'{k} "{v}";')
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Minimal Ensembl gene block (+ strand, protein-coding)
# ---------------------------------------------------------------------------
#
#  ENSG00000000001.3  /  BRCA1
#    transcript ENST00000000001  version 5  → MANE_Select
#      exon  100–300  exon_id ENSE00000000001
#      exon  500–700  exon_id ENSE00000000002
#      start_codon 150
#      stop_codon  end=650
#      five_prime_utr  100–149
#      three_prime_utr 651–700
#    transcript ENST00000000002  version 2  → MANE_Plus_Clinical
#      exon  100–300  exon_id ENSE00000000001  (duplicate of above)
#      exon  800–900  exon_id ENSE00000000003
#      stop_codon  end=880

ENSEMBL_GENE1_LINES = [
    _gtf_line("1", "ensembl_havana", "gene", "100", "900", "+",
              _attrs(gene_id="ENSG00000000001.3", gene_version="3",
                     gene_name="BRCA1", gene_biotype="protein_coding")),

    # transcript 1 — MANE_Select
    _gtf_line("1", "ensembl_havana", "transcript", "100", "700", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000001",
                     transcript_version="5", transcript_biotype="protein_coding",
                     tag=["MANE_Select", "basic"])),

    _gtf_line("1", "ensembl_havana", "exon", "100", "300", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000001",
                     transcript_version="5", transcript_biotype="protein_coding",
                     exon_id="ENSE00000000001")),

    _gtf_line("1", "ensembl_havana", "exon", "500", "700", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000001",
                     transcript_version="5", transcript_biotype="protein_coding",
                     exon_id="ENSE00000000002")),

    _gtf_line("1", "ensembl_havana", "start_codon", "150", "152", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000001",
                     transcript_version="5", transcript_biotype="protein_coding")),

    _gtf_line("1", "ensembl_havana", "stop_codon", "648", "650", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000001",
                     transcript_version="5", transcript_biotype="protein_coding")),

    _gtf_line("1", "ensembl_havana", "five_prime_utr", "100", "149", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000001",
                     transcript_version="5", transcript_biotype="protein_coding")),

    _gtf_line("1", "ensembl_havana", "three_prime_utr", "651", "700", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000001",
                     transcript_version="5", transcript_biotype="protein_coding")),

    # transcript 2 — MANE_Plus_Clinical
    _gtf_line("1", "ensembl_havana", "transcript", "100", "900", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000002",
                     transcript_version="2", transcript_biotype="protein_coding",
                     tag=["MANE_Plus_Clinical", "basic"])),

    _gtf_line("1", "ensembl_havana", "exon", "100", "300", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000002",
                     transcript_version="2", transcript_biotype="protein_coding",
                     exon_id="ENSE00000000001")),   # duplicate exon_id

    _gtf_line("1", "ensembl_havana", "exon", "800", "900", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000002",
                     transcript_version="2", transcript_biotype="protein_coding",
                     exon_id="ENSE00000000003")),

    _gtf_line("1", "ensembl_havana", "stop_codon", "878", "880", "+",
              _attrs(gene_id="ENSG00000000001.3", transcript_id="ENST00000000002",
                     transcript_version="2", transcript_biotype="protein_coding")),
]

# ---------------------------------------------------------------------------
# Minimal Ensembl gene block (− strand)
# ---------------------------------------------------------------------------
#
#  ENSG00000000002.1  /  TP53
#    transcript ENST00000000003  version 1
#      exon  200–400  exon_id ENSE00000000004
#      exon  600–800  exon_id ENSE00000000005
#      stop_codon  start=210, end=212  (− strand: trims exon[0])

ENSEMBL_GENE2_LINES = [
    _gtf_line("17", "ensembl_havana", "gene", "200", "800", "-",
              _attrs(gene_id="ENSG00000000002.1", gene_version="1",
                     gene_name="TP53", gene_biotype="protein_coding")),

    _gtf_line("17", "ensembl_havana", "transcript", "200", "800", "-",
              _attrs(gene_id="ENSG00000000002.1", transcript_id="ENST00000000003",
                     transcript_version="1", transcript_biotype="protein_coding",
                     tag=["basic"])),

    _gtf_line("17", "ensembl_havana", "exon", "200", "400", "-",
              _attrs(gene_id="ENSG00000000002.1", transcript_id="ENST00000000003",
                     transcript_version="1", transcript_biotype="protein_coding",
                     exon_id="ENSE00000000004")),

    _gtf_line("17", "ensembl_havana", "exon", "600", "800", "-",
              _attrs(gene_id="ENSG00000000002.1", transcript_id="ENST00000000003",
                     transcript_version="1", transcript_biotype="protein_coding",
                     exon_id="ENSE00000000005")),

    # stop_codon on − strand: end=212 → exon [200,400] → exon[0] trimmed to 212
    _gtf_line("17", "ensembl_havana", "stop_codon", "210", "212", "-",
              _attrs(gene_id="ENSG00000000002.1", transcript_id="ENST00000000003",
                     transcript_version="1", transcript_biotype="protein_coding")),
]

# ---------------------------------------------------------------------------
# Minimal GENCODE gene block (+ strand)
# ---------------------------------------------------------------------------
#
#  ENSG00000000003.5  /  EGFR
#    transcript ENST00000000004.3  → MANE_Select + MANE_Plus_Clinical (both)
#      exon  1000–1200  exon_id ENSE00000000006.1
#      exon  1400–1600  exon_id ENSE00000000007.2
#      start_codon 1050
#      stop_codon  end=1580
#      UTR  1000–1049   (GENCODE unified UTR feature)
#      UTR  1581–1600
#    transcript ENST00000000005.1  → no MANE tag, non-protein_coding exon

GENCODE_GENE1_LINES = [
    _gtf_line("chr7", "HAVANA", "gene", "1000", "1600", "+",
              _attrs(gene_id="ENSG00000000003.5", gene_type="protein_coding",
                     gene_name="EGFR")),

    # transcript 1 — both MANE tags (unusual but valid edge case)
    _gtf_line("chr7", "HAVANA", "transcript", "1000", "1600", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000004.3",
                     transcript_type="protein_coding",
                     transcript_name="EGFR-201",
                     tag=["MANE_Select", "MANE_Plus_Clinical", "basic",
                          "appris_principal_1", "CCDS"])),

    _gtf_line("chr7", "HAVANA", "exon", "1000", "1200", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000004.3",
                     transcript_type="protein_coding",
                     exon_id="ENSE00000000006.1")),

    _gtf_line("chr7", "HAVANA", "exon", "1400", "1600", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000004.3",
                     transcript_type="protein_coding",
                     exon_id="ENSE00000000007.2")),

    _gtf_line("chr7", "HAVANA", "start_codon", "1050", "1052", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000004.3",
                     transcript_type="protein_coding")),

    _gtf_line("chr7", "HAVANA", "stop_codon", "1578", "1580", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000004.3",
                     transcript_type="protein_coding")),

    _gtf_line("chr7", "HAVANA", "UTR", "1000", "1049", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000004.3",
                     transcript_type="protein_coding")),

    _gtf_line("chr7", "HAVANA", "UTR", "1581", "1600", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000004.3",
                     transcript_type="protein_coding")),

    # transcript 2 — no MANE tag, lncRNA (exon must be ignored)
    _gtf_line("chr7", "ENSEMBL", "transcript", "1000", "1600", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000005.1",
                     transcript_type="lncRNA",
                     transcript_name="EGFR-202")),

    _gtf_line("chr7", "ENSEMBL", "exon", "1000", "1200", "+",
              _attrs(gene_id="ENSG00000000003.5",
                     transcript_id="ENST00000000005.1",
                     transcript_type="lncRNA",
                     exon_id="ENSE00000000008.1")),
]

# ---------------------------------------------------------------------------
# Transcript-only annotation for prepare_annot_for_muts
# ---------------------------------------------------------------------------

def _make_muts_ensembl_annot():
    """Ensembl-style annotation for prepare_annot_for_muts tests."""
    return [
        _gtf_line("2", "ensembl_havana", "transcript", "1000", "2000", "+",
                  _attrs(gene_id="ENSG00000000010.1", transcript_id="ENST00000000010",
                         transcript_version="3", transcript_biotype="protein_coding")),
        _gtf_line("2", "ensembl_havana", "exon", "1000", "1200", "+",
                  _attrs(gene_id="ENSG00000000010.1", transcript_id="ENST00000000010",
                         transcript_version="3", transcript_biotype="protein_coding",
                         exon_id="ENSE00000000010")),
        _gtf_line("2", "ensembl_havana", "exon", "1500", "2000", "+",
                  _attrs(gene_id="ENSG00000000010.1", transcript_id="ENST00000000010",
                         transcript_version="3", transcript_biotype="protein_coding",
                         exon_id="ENSE00000000011")),
        _gtf_line("2", "ensembl_havana", "start_codon", "1050", "1052", "+",
                  _attrs(gene_id="ENSG00000000010.1", transcript_id="ENST00000000010",
                         transcript_version="3", transcript_biotype="protein_coding")),
    ]


def _make_muts_gencode_annot():
    """GENCODE-style annotation for prepare_annot_for_muts tests."""
    return [
        _gtf_line("chr2", "HAVANA", "transcript", "1000", "2000", "-",
                  _attrs(gene_id="ENSG00000000010.1",
                         transcript_id="ENST00000000010.3",
                         transcript_type="protein_coding")),
        _gtf_line("chr2", "HAVANA", "exon", "1000", "1200", "-",
                  _attrs(gene_id="ENSG00000000010.1",
                         transcript_id="ENST00000000010.3",
                         transcript_type="protein_coding",
                         exon_id="ENSE00000000010.1")),
        _gtf_line("chr2", "HAVANA", "start_codon", "1180", "1182", "-",
                  _attrs(gene_id="ENSG00000000010.1",
                         transcript_id="ENST00000000010.3",
                         transcript_type="protein_coding")),
    ]


# ===========================================================================
# Tests: internal helpers
# ===========================================================================

class TestParseAttributes:
    def test_single_quoted_value(self):
        attrs = _parse_attributes('gene_id "ENSG00000000001.3";')
        assert attrs["gene_id"] == ["ENSG00000000001.3"]

    def test_multiple_different_keys(self):
        attrs = _parse_attributes(
            'gene_id "ENSG1"; transcript_id "ENST1"; gene_name "BRCA1";'
        )
        assert attrs["gene_id"] == ["ENSG1"]
        assert attrs["transcript_id"] == ["ENST1"]
        assert attrs["gene_name"] == ["BRCA1"]

    def test_multi_value_tag_key(self):
        attrs = _parse_attributes(
            'tag "MANE_Select"; tag "basic"; tag "CCDS";'
        )
        assert attrs["tag"] == ["MANE_Select", "basic", "CCDS"]

    def test_all_mane_tags_present(self):
        attrs = _parse_attributes(
            'tag "MANE_Select"; tag "MANE_Plus_Clinical"; tag "basic";'
        )
        assert "MANE_Select" in attrs["tag"]
        assert "MANE_Plus_Clinical" in attrs["tag"]

    def test_empty_string(self):
        assert _parse_attributes("") == {}

    def test_no_tag_key(self):
        attrs = _parse_attributes('gene_id "ENSG1";')
        assert "tag" not in attrs


class TestGetAttr:
    def setup_method(self):
        self.attrs = {"gene_id": ["ENSG1"], "gene_name": ["BRCA1"]}

    def test_key_present(self):
        assert _get_attr(self.attrs, "gene_id") == "ENSG1"

    def test_key_absent_returns_default(self):
        assert _get_attr(self.attrs, "missing") == ""

    def test_custom_default(self):
        assert _get_attr(self.attrs, "missing", default="X") == "X"

    def test_fallback_chain_first_wins(self):
        assert _get_attr(self.attrs, "gene_id", "gene_name") == "ENSG1"

    def test_fallback_chain_second_wins(self):
        assert _get_attr(self.attrs, "missing", "gene_name") == "BRCA1"

    def test_fallback_chain_all_missing(self):
        assert _get_attr(self.attrs, "a", "b", "c", default="Z") == "Z"


class TestHasTag:
    def test_tag_present(self):
        attrs = {"tag": ["MANE_Select", "basic"]}
        assert _has_tag(attrs, "MANE_Select") is True

    def test_tag_absent(self):
        attrs = {"tag": ["basic"]}
        assert _has_tag(attrs, "MANE_Select") is False

    def test_no_tag_key_at_all(self):
        assert _has_tag({}, "MANE_Select") is False

    def test_mane_plus_clinical(self):
        attrs = {"tag": ["MANE_Plus_Clinical"]}
        assert _has_tag(attrs, "MANE_Plus_Clinical") is True
        assert _has_tag(attrs, "MANE_Select") is False


class TestDetectSource:
    def test_ensembl_bare_chrom(self):
        annot = [_gtf_line("1", "ensembl", "gene", "1", "100", "+", 'gene_id "G1";')]
        assert _detect_source(annot) == "ensembl"

    def test_gencode_chr_prefix(self):
        annot = [_gtf_line("chr1", "HAVANA", "gene", "1", "100", "+", 'gene_id "G1";')]
        assert _detect_source(annot) == "gencode"

    def test_mixed_first_chr_line_wins(self):
        # First data line is chr-prefixed → gencode
        annot = [
            _gtf_line("chr1", "HAVANA", "gene", "1", "100", "+", 'gene_id "G1";'),
            _gtf_line("1", "ensembl", "gene", "1", "100", "+", 'gene_id "G2";'),
        ]
        assert _detect_source(annot) == "gencode"

    def test_empty_annotation_defaults_to_ensembl(self):
        assert _detect_source([]) == "ensembl"

    def test_chrM_detected_as_gencode(self):
        annot = [_gtf_line("chrM", "HAVANA", "gene", "1", "100", "+", 'gene_id "G1";')]
        assert _detect_source(annot) == "gencode"

    def test_chrX_detected_as_gencode(self):
        annot = [_gtf_line("chrX", "HAVANA", "gene", "1", "100", "+", 'gene_id "G1";')]
        assert _detect_source(annot) == "gencode"


class TestBiotypeKey:
    def test_ensembl_transcript(self):
        assert _biotype_key("ensembl", "transcript") == "transcript_biotype"

    def test_ensembl_gene(self):
        assert _biotype_key("ensembl", "gene") == "gene_biotype"

    def test_gencode_transcript(self):
        assert _biotype_key("gencode", "transcript") == "transcript_type"

    def test_gencode_gene(self):
        assert _biotype_key("gencode", "gene") == "gene_type"


class TestTranscriptIdVersioned:
    def test_ensembl_bare_id_with_version_attr(self):
        attrs = {"transcript_id": ["ENST00000000001"], "transcript_version": ["5"]}
        assert _transcript_id_versioned(attrs, "ensembl") == "ENST00000000001.5"

    def test_ensembl_already_dotted_id_unchanged(self):
        # If Ensembl ID somehow already has a dot, don't double-append
        attrs = {"transcript_id": ["ENST00000000001.5"]}
        assert _transcript_id_versioned(attrs, "ensembl") == "ENST00000000001.5"

    def test_ensembl_missing_version_attr(self):
        attrs = {"transcript_id": ["ENST00000000001"]}
        assert _transcript_id_versioned(attrs, "ensembl") == "ENST00000000001"

    def test_gencode_version_embedded(self):
        attrs = {"transcript_id": ["ENST00000000004.3"]}
        assert _transcript_id_versioned(attrs, "gencode") == "ENST00000000004.3"

    def test_gencode_no_version_suffix(self):
        # Unusual but should not crash
        attrs = {"transcript_id": ["ENST00000000004"]}
        assert _transcript_id_versioned(attrs, "gencode") == "ENST00000000004"


class TestIsUtrFeature:
    def test_five_prime_utr(self):
        assert _is_utr_feature("five_prime_utr") is True

    def test_three_prime_utr(self):
        assert _is_utr_feature("three_prime_utr") is True

    def test_gencode_utr(self):
        assert _is_utr_feature("UTR") is True

    def test_exon_is_not_utr(self):
        assert _is_utr_feature("exon") is False

    def test_cds_is_not_utr(self):
        assert _is_utr_feature("CDS") is False

    def test_gene_is_not_utr(self):
        assert _is_utr_feature("gene") is False


# ===========================================================================
# Tests: prepare_annot_for_muts
# ===========================================================================

class TestPrepareAnnotForMuts:

    # --- triggering / skipping logic ----------------------------------------

    def test_returns_none_for_genomic_only_mutations(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["BRCA1:g.12345A>G"], annot, logger)
        assert result is None

    def test_returns_none_for_empty_mutation_list(self):
        annot = _make_muts_ensembl_annot()
        assert prepare_annot_for_muts([], annot, logger) is None

    def test_returns_none_for_mutations_without_colon(self):
        annot = _make_muts_ensembl_annot()
        assert prepare_annot_for_muts(["BRCA1"], annot, logger) is None

    def test_c_notation_triggers_parsing(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["GENE:c.123G>C"], annot, logger)
        assert result is not None

    def test_n_notation_triggers_parsing(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["GENE:n.456A>T"], annot, logger)
        assert result is not None

    def test_mixed_list_one_hgvs_triggers_parsing(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["BRCA1:g.1A>G", "TP53:c.100G>A"], annot, logger)
        assert result is not None

    # --- Ensembl parsing ----------------------------------------------------

    def test_ensembl_versioned_id_assembled(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        assert "ENST00000000010.3" in result

    def test_ensembl_chrom_preserved(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        chrom, *_ = result["ENST00000000010.3"]
        assert chrom == "2"

    def test_ensembl_exon_list_correct(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        _, exon_list, _, _ = result["ENST00000000010.3"]
        assert ["1000", "1200"] in exon_list
        assert ["1500", "2000"] in exon_list

    def test_ensembl_start_codon_coord(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        _, _, start_codon, _ = result["ENST00000000010.3"]
        assert start_codon == "1050"

    def test_ensembl_strand_correct(self):
        annot = _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        _, _, _, strand = result["ENST00000000010.3"]
        assert strand == "+"

    # --- GENCODE parsing ----------------------------------------------------

    def test_gencode_versioned_id_from_embedded_version(self):
        annot = _make_muts_gencode_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        assert "ENST00000000010.3" in result

    def test_gencode_chr_prefix_preserved(self):
        annot = _make_muts_gencode_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        chrom, *_ = result["ENST00000000010.3"]
        assert chrom == "chr2"

    def test_gencode_minus_strand(self):
        annot = _make_muts_gencode_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        _, _, _, strand = result["ENST00000000010.3"]
        assert strand == "-"

    # --- exclusion of transcripts without start_codon -----------------------

    def test_transcript_without_start_codon_excluded(self):
        annot = [
            _gtf_line("1", "ensembl", "transcript", "1", "500", "+",
                      _attrs(gene_id="ENSG1", transcript_id="ENST_NO_START",
                             transcript_version="1",
                             transcript_biotype="protein_coding")),
            _gtf_line("1", "ensembl", "exon", "1", "500", "+",
                      _attrs(gene_id="ENSG1", transcript_id="ENST_NO_START",
                             transcript_version="1",
                             transcript_biotype="protein_coding",
                             exon_id="ENSE1")),
            # no start_codon line
        ]
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        assert result == {}

    # --- multiple transcripts -----------------------------------------------

    def test_multiple_transcripts_all_collected(self):
        annot = _make_muts_ensembl_annot() + [
            _gtf_line("2", "ensembl", "transcript", "3000", "4000", "-",
                      _attrs(gene_id="ENSG2", transcript_id="ENST00000000020",
                             transcript_version="1",
                             transcript_biotype="protein_coding")),
            _gtf_line("2", "ensembl", "exon", "3000", "3500", "-",
                      _attrs(gene_id="ENSG2", transcript_id="ENST00000000020",
                             transcript_version="1",
                             transcript_biotype="protein_coding",
                             exon_id="ENSE20")),
            _gtf_line("2", "ensembl", "start_codon", "3450", "3452", "-",
                      _attrs(gene_id="ENSG2", transcript_id="ENST00000000020",
                             transcript_version="1",
                             transcript_biotype="protein_coding")),
        ]
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        assert "ENST00000000010.3" in result
        assert "ENST00000000020.1" in result

    # --- malformed lines ----------------------------------------------------

    def test_malformed_short_lines_skipped(self):
        annot = ["this is not a valid gtf line"] + _make_muts_ensembl_annot()
        result = prepare_annot_for_muts(["GENE:c.1A>T"], annot, logger)
        assert "ENST00000000010.3" in result


# ===========================================================================
# Tests: prepare_annot_for_genes
# ===========================================================================

class TestPrepareAnnotForGenes:

    # --- Ensembl: basic structure -------------------------------------------

    def test_ensembl_gene_id_stripped_of_version(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        assert "ENSG00000000001" in gene_dict
        assert "ENSG00000000001.3" not in gene_dict

    def test_ensembl_gene_name_dict_populated(self):
        _, gene_name_dict = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        assert gene_name_dict["BRCA1"] == "ENSG00000000001"

    def test_ensembl_chrom_correct(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        chrom = gene_dict["ENSG00000000001"][0]
        assert chrom == "1"

    def test_ensembl_gene_coords(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        coords = gene_dict["ENSG00000000001"][1]
        assert coords == ["100", "900"]

    # --- Ensembl: exon collection and deduplication -------------------------

    def test_ensembl_exon_deduplication(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        exon_id_list = gene_dict["ENSG00000000001"][4]
        assert len(exon_id_list) == len(set(exon_id_list)), "Duplicate exon IDs found"

    def test_ensembl_exon_list_sorted(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        exon_list = gene_dict["ENSG00000000001"][2]
        starts = [e[0] for e in exon_list]
        assert starts == sorted(starts)

    # --- Ensembl: stop-codon trimming on + strand ---------------------------

    def test_ensembl_plus_strand_stop_codon_trims_exon_end(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        exon_list = gene_dict["ENSG00000000001"][2]
        # stop_codon end=650 is inside exon [500,700] → exon[1] trimmed to 650
        exon_500 = next(e for e in exon_list if e[0] == 500)
        assert exon_500[1] == 650

    # --- Ensembl: stop-codon trimming on − strand ---------------------------

    def test_ensembl_minus_strand_stop_codon_trims_exon_start(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE2_LINES, logger)
        exon_list = gene_dict["ENSG00000000002"][2]
        # stop_codon end=212 is inside exon [200,400] → exon[0] trimmed to 212
        exon_200 = next(e for e in exon_list if e[1] == 400)
        assert exon_200[0] == 212

    def test_ensembl_minus_strand_unaffected_exon_unchanged(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE2_LINES, logger)
        exon_list = gene_dict["ENSG00000000002"][2]
        exon_600 = next(e for e in exon_list if e[0] == 600)
        assert exon_600 == [600, 800]

    # --- Ensembl: UTR collection (five_prime_utr / three_prime_utr) ---------

    def test_ensembl_utr_collected(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        utr_list = gene_dict["ENSG00000000001"][3]
        assert ["100", "149"] in utr_list
        assert ["651", "700"] in utr_list

    def test_ensembl_utr_deduplication(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        utr_list = gene_dict["ENSG00000000001"][3]
        utr_tuples = [tuple(u) for u in utr_list]
        assert len(utr_tuples) == len(set(utr_tuples))

    # --- Ensembl: MANE tags -------------------------------------------------

    def test_ensembl_mane_select_collected(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        mane_sel = gene_dict["ENSG00000000001"][5]
        assert "ENST00000000001.5" in mane_sel

    def test_ensembl_mane_plus_clinical_collected(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE1_LINES, logger)
        mane_plus = gene_dict["ENSG00000000001"][6]
        assert "ENST00000000002.2" in mane_plus

    def test_ensembl_non_mane_transcript_not_in_mane_lists(self):
        gene_dict, _ = prepare_annot_for_genes(ENSEMBL_GENE2_LINES, logger)
        mane_sel = gene_dict["ENSG00000000002"][5]
        mane_plus = gene_dict["ENSG00000000002"][6]
        assert mane_sel == []
        assert mane_plus == []

    # --- GENCODE: basic structure -------------------------------------------

    def test_gencode_gene_id_stripped_of_version(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        assert "ENSG00000000003" in gene_dict
        assert "ENSG00000000003.5" not in gene_dict

    def test_gencode_chr_prefix_preserved(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        chrom = gene_dict["ENSG00000000003"][0]
        assert chrom == "chr7"

    def test_gencode_gene_name_dict(self):
        _, gene_name_dict = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        assert gene_name_dict["EGFR"] == "ENSG00000000003"

    # --- GENCODE: biotype filtering (transcript_type) -----------------------

    def test_gencode_non_protein_coding_exon_excluded(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        exon_id_list = gene_dict["ENSG00000000003"][4]
        # ENSE00000000008.1 belongs to lncRNA transcript → must be absent
        assert "ENSE00000000008.1" not in exon_id_list

    def test_gencode_protein_coding_exons_collected(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        exon_list = gene_dict["ENSG00000000003"][2]
        assert len(exon_list) == 2

    # --- GENCODE: UTR (unified "UTR" feature name) --------------------------

    def test_gencode_utr_feature_collected(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        utr_list = gene_dict["ENSG00000000003"][3]
        assert ["1000", "1049"] in utr_list
        assert ["1581", "1600"] in utr_list

    # --- GENCODE: stop-codon trimming on + strand ---------------------------

    def test_gencode_plus_strand_stop_codon_trims_exon_end(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        exon_list = gene_dict["ENSG00000000003"][2]
        # stop_codon end=1580 inside exon [1400,1600] → exon[1] trimmed to 1580
        exon_1400 = next(e for e in exon_list if e[0] == 1400)
        assert exon_1400[1] == 1580

    # --- GENCODE: MANE tags -------------------------------------------------

    def test_gencode_mane_select_collected(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        mane_sel = gene_dict["ENSG00000000003"][5]
        assert "ENST00000000004.3" in mane_sel

    def test_gencode_mane_plus_clinical_collected(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        mane_plus = gene_dict["ENSG00000000003"][6]
        assert "ENST00000000004.3" in mane_plus

    def test_gencode_transcript_with_both_mane_tags(self):
        """A single transcript tagged with both MANE_Select and MANE_Plus_Clinical
        must appear in both lists."""
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        mane_sel = gene_dict["ENSG00000000003"][5]
        mane_plus = gene_dict["ENSG00000000003"][6]
        assert "ENST00000000004.3" in mane_sel
        assert "ENST00000000004.3" in mane_plus

    def test_gencode_lncrna_transcript_not_in_mane_lists(self):
        gene_dict, _ = prepare_annot_for_genes(GENCODE_GENE1_LINES, logger)
        mane_sel = gene_dict["ENSG00000000003"][5]
        mane_plus = gene_dict["ENSG00000000003"][6]
        assert "ENST00000000005.1" not in mane_sel
        assert "ENST00000000005.1" not in mane_plus

    # --- edge cases ---------------------------------------------------------

    def test_gene_with_no_protein_coding_exons(self):
        annot = [
            _gtf_line("chr1", "HAVANA", "gene", "1", "500", "+",
                      _attrs(gene_id="ENSG99.1", gene_type="lncRNA",
                             gene_name="LINC001")),
            _gtf_line("chr1", "HAVANA", "transcript", "1", "500", "+",
                      _attrs(gene_id="ENSG99.1", transcript_id="ENST99.1",
                             transcript_type="lncRNA")),
            _gtf_line("chr1", "HAVANA", "exon", "1", "500", "+",
                      _attrs(gene_id="ENSG99.1", transcript_id="ENST99.1",
                             transcript_type="lncRNA", exon_id="ENSE99.1")),
        ]
        gene_dict, _ = prepare_annot_for_genes(annot, logger)
        assert gene_dict["ENSG99"][2] == []

    def test_gene_with_no_utrs(self):
        annot = [
            _gtf_line("1", "ensembl", "gene", "1", "500", "+",
                      _attrs(gene_id="ENSG88.1", gene_biotype="protein_coding",
                             gene_name="GENE88")),
            _gtf_line("1", "ensembl", "transcript", "1", "500", "+",
                      _attrs(gene_id="ENSG88.1", transcript_id="ENST88",
                             transcript_version="1",
                             transcript_biotype="protein_coding")),
            _gtf_line("1", "ensembl", "exon", "1", "500", "+",
                      _attrs(gene_id="ENSG88.1", transcript_id="ENST88",
                             transcript_version="1",
                             transcript_biotype="protein_coding",
                             exon_id="ENSE88")),
        ]
        gene_dict, _ = prepare_annot_for_genes(annot, logger)
        assert gene_dict["ENSG88"][3] == []

    def test_gene_with_no_stop_codon_no_trimming(self):
        annot = [
            _gtf_line("1", "ensembl", "gene", "1", "500", "+",
                      _attrs(gene_id="ENSG77.1", gene_biotype="protein_coding",
                             gene_name="GENE77")),
            _gtf_line("1", "ensembl", "transcript", "1", "500", "+",
                      _attrs(gene_id="ENSG77.1", transcript_id="ENST77",
                             transcript_version="1",
                             transcript_biotype="protein_coding")),
            _gtf_line("1", "ensembl", "exon", "100", "400", "+",
                      _attrs(gene_id="ENSG77.1", transcript_id="ENST77",
                             transcript_version="1",
                             transcript_biotype="protein_coding",
                             exon_id="ENSE77")),
        ]
        gene_dict, _ = prepare_annot_for_genes(annot, logger)
        exon_list = gene_dict["ENSG77"][2]
        assert exon_list == [[100, 400]]  # untrimmed

    def test_two_consecutive_genes_both_parsed(self):
        annot = ENSEMBL_GENE1_LINES + ENSEMBL_GENE2_LINES
        gene_dict, gene_name_dict = prepare_annot_for_genes(annot, logger)
        assert "ENSG00000000001" in gene_dict
        assert "ENSG00000000002" in gene_dict
        assert gene_name_dict["BRCA1"] == "ENSG00000000001"
        assert gene_name_dict["TP53"] == "ENSG00000000002"

    def test_malformed_short_lines_skipped(self):
        annot = ["bad line\twithout\tnine\tcolumns"] + ENSEMBL_GENE1_LINES
        gene_dict, _ = prepare_annot_for_genes(annot, logger)
        assert "ENSG00000000001" in gene_dict

    def test_multiple_mane_select_transcripts_per_gene(self):
        """Edge case: two transcripts both tagged MANE_Select."""
        annot = [
            _gtf_line("1", "ensembl", "gene", "1", "1000", "+",
                      _attrs(gene_id="ENSG55.1", gene_biotype="protein_coding",
                             gene_name="GENE55")),
            _gtf_line("1", "ensembl", "transcript", "1", "500", "+",
                      _attrs(gene_id="ENSG55.1", transcript_id="ENST55A",
                             transcript_version="1",
                             transcript_biotype="protein_coding",
                             tag=["MANE_Select"])),
            _gtf_line("1", "ensembl", "transcript", "500", "1000", "+",
                      _attrs(gene_id="ENSG55.1", transcript_id="ENST55B",
                             transcript_version="2",
                             transcript_biotype="protein_coding",
                             tag=["MANE_Select"])),
        ]
        gene_dict, _ = prepare_annot_for_genes(annot, logger)
        mane_sel = gene_dict["ENSG55"][5]
        assert "ENST55A.1" in mane_sel
        assert "ENST55B.2" in mane_sel
