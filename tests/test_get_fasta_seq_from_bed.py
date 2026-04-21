import pytest
from contextlib import nullcontext as does_not_raise
from Bio.SeqRecord import SeqRecord

import logging

from OligoDesigner.utils import *
from OligoDesigner.get_regions_for_mutations import *
from OligoDesigner.get_regions_for_genes import *
from OligoDesigner.get_sequences_from_bed import *
from OligoDesigner.constants import genome_versions_chromosome_dict

@pytest.fixture
def logger():
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    return logger

@pytest.mark.parametrize(
    "coords_list, chrom, answ, clas",
   [
        (
            ['chr1\t5\t10\ttest_g\texo3\texon\t1'],
            'ATGCTGATCGTAGCATCGATCGATCGATCGATCGATGCATCG',
            ['>Region chrom:chr1 start:5 end:10 gene:test_g exon:exo3', 'TGATCG'],
            list,
        ),
        (
            ['chr1\t5\t9\tG\tchr1:7A>T\tsnp\t1'],
            'ATGCTGATCGTAGCATCGATCGATCGATCGATCGATGCATCG',
            ['>SNP chrom:chr1 start:5 end:9 ref_allele:G snp_id:chr1:7A>T', 'TGATC'],
            list,
        ),
        (
            ['chr1\t5\t9\tA\tchr1:7A>T\tsnp\t1'],
            'ATGCTGATCGTAGCATCGATCGATCGATCGATCGATGCATCG',
            [],
            list,
        ),
        (
            ['chr1\t5\t9\tdelA\tchr1:7A>T\tsnp\t1'],
            'ATGCTGATCGTAGCATCGATCGATCGATCGATCGATGCATCG',
            ['>SNP chrom:chr1 start:5 end:9 ref_allele:delA snp_id:chr1:7A>T', 'TGATC'],
            list,
        ),
    ]
)
def test_create_fasta_sequences(coords_list, chrom, answ, clas, logger):
    result = create_fasta_sequences(coords_list, chrom, logger)
    assert isinstance(result, clas)

@pytest.mark.parametrize(
    "coords_list, genome, answ",
    [
        (['chr1\t5\t10\ttest_g\texo3\texon\t1'], './tests/FileTest_genome.fasta', ['>Region chrom:chr1 start:5 end:10 gene:test_g exon:exo3', 'TCGATC']),
        (['chr1\t4\t15\ttest_g\texo3\texon\t1', 'chr3\t10\t20\ttest_g\texo3\texon\t1'], './tests/FileTest_genome.fasta', ['>Region chrom:chr1 start:4 end:15 gene:test_g exon:exo3', 'ATCGATCACGAC', '>Region chrom:chr3 start:10 end:20 gene:test_g exon:exo3', 'CTACGATCGAC'])
    ]
)
def test_get_fasta_seqs_from_bed(coords_list, genome, answ, logger):
    result = get_fasta_seqs_from_bed(coords_list, genome, 'GRCh38.p14', False, '', logger)
    assert isinstance(result[0], SeqRecord)