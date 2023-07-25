import pytest
from contextlib import nullcontext as does_not_raise

import logging

from OligoDesigner.utils import *
from OligoDesigner.get_regions_for_mutations import *
from OligoDesigner.get_regions_for_genes import *
from OligoDesigner.constants import genome_versions_chromosome_dict

@pytest.fixture
def annot():
    test_gtf = './tests/FileTest_annot.gtf'
    with open(test_gtf) as infile:
        annot = [line.strip() for line in infile.readlines() if not line.startswith('#')]
    return annot

@pytest.fixture
def logger():
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    return logger

def test_prepare_annot_for_muts(annot):
    annot_test = prepare_annot_for_muts(['Chr2:c.123G>C'], annot, '')
    answer = {'ENST00000673477.1': ['1', [['1471765', '1472089'], ['1477274', '1477350'], ['1478644', '1478745'], ['1479049', '1479108'], ['1480867', '1480936'], ['1482138', '1482303'], ['1482545', '1482614'], ['1485016', '1485171'], ['1485782', '1485838'], ['1486110', '1486235'], ['1486544', '1486668'], ['1487863', '1487914'], ['1489204', '1489274'], ['1490257', '1490424'], ['1490563', '1490671'], ['1495485', '1497848']]], 'test.7': ['1', [['10', '90'], ['100', '200'], ['400', '500']]], 'ENST00000612315.1': ['KI270713.1', [['31698', '32528']]], 'ENST00000619109.1': ['KI270713.1', [['21861', '22024']]]}

    assert annot_test == answer

def test_prepare_annot_for_muts_full(annot, logger):
    annot_test = prepare_annot_for_muts(['Chr2:c.123G>C'], annot, '')
    muts_list = ['test.7:c.47-2_50-2insATGC']
    answer = get_regions_for_muts(muts_list, annot_test, 10, False, "", logger, 'GRCh38')

    assert answer == ['chr1\t45\t65\tinsATGC\ttest.7:c.47-2_50-2insATGC\tsnp\t1']

def test_prepare_annot_for_genes(annot, logger):
    gene_dict, gene_name_dict = prepare_annot_for_genes(annot, logger)
    gene_dict_answ = {'ENSG00000160072': ['1', ['1471765', '1497848'], [['1471765', '1472089'], ['1477274', '1477350'], ['1478644', '1478745'], ['1479049', '1479108'], ['1480867', '1480936'], ['1482138', '1482303'], ['1482545', '1482614'], ['1485016', '1485171'], ['1485782', '1485838'], ['1486110', '1486235'], ['1486544', '1486668'], ['1487863', '1487914'], ['1489204', '1489274'], ['1490257', '1490424'], ['1490563', '1490671'], ['1495485', '1497848']], [['1471765', '1471884'], ['1495818', '1497848']], ['ENSE00003889014', 'ENSE00003467707', 'ENSE00003569130', 'ENSE00003608502', 'ENSE00003474888', 'ENSE00003654064', 'ENSE00003510521', 'ENSE00003459370', 'ENSE00003655926', 'ENSE00003594545', 'ENSE00003892109', 'ENSE00003689846', 'ENSE00003670332', 'ENSE00003505365', 'ENSE00003497242', 'ENSE00003892597']], 'test_g': ['1', ['10', '500'], [['10', '90'], ['100', '200'], ['400', '500']], [['1', '9'], ['501', '600']], ['exo1', 'exo2', 'exo3']], 'ENSG00000277475': ['KI270713.1', ['31968', '32528'], [['31698', '32528']], [['32373', '32528'], ['31698', '31841']], ['ENSE00003739404']], 'ENSG00000275405': ['KI270713.1', ['21861', '22024'], [['21861', '22024']], [], ['ENSE00003747303']]}
    gene_dict_name_answ = {'ATAD3B': 'ENSG00000160072', 'test_name': 'test_g', 'GENENAME': 'ENSG00000277475', 'U1': 'ENSG00000275405'}

    assert gene_dict == gene_dict_answ and gene_name_dict == gene_dict_name_answ

def test_prepare_annot_for_genes_full(annot, logger):
    gene_dict, gene_name_dict = prepare_annot_for_genes(annot, logger)
    input_genes_list = ["ATAD3B\t1-2,4-5\t100","test_g\texon\t0"]

    gene_coords_regions = get_regions_for_genes(input_genes_list, gene_dict, gene_name_dict,
                            'exon', 50, False, '', logger, 'GRCh38')

    answ = ['chr1\t1471665\t1472189\tENSG00000160072\tENSE00003889014\t1-2\t1', 'chr1\t1477174\t1477450\tENSG00000160072\tENSE00003467707\t1-2\t1',\
            'chr1\t1478949\t1479208\tENSG00000160072\tENSE00003608502\t4-5\t1', 'chr1\t1480767\t1481036\tENSG00000160072\tENSE00003474888\t4-5\t1',\
            'chr1\t10\t90\ttest_g\texo1\texon\t1', 'chr1\t100\t200\ttest_g\texo2\texon\t1',\
            'chr1\t400\t500\ttest_g\texo3\texon\t1']

    assert gene_coords_regions == answ