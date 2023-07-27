import pytest
from contextlib import nullcontext as does_not_raise

import logging

from OligoDesigner.utils import *
from OligoDesigner.get_regions_for_mutations import *
from OligoDesigner.get_regions_for_genes import *
from OligoDesigner.constants import genome_versions_chromosome_dict

@pytest.fixture
def logger():
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    return logger

@pytest.mark.parametrize(
    "input_genes_list, answer, ext",
    [
        (["ATAD3B\t1-2,4-5\t100"], ['chr1\t1471665\t1472189\tENSG00000160072\tENSE00003889014\t1-2\t1', 'chr1\t1477174\t1477450\tENSG00000160072\tENSE00003467707\t1-2\t1',\
            'chr1\t1478949\t1479208\tENSG00000160072\tENSE00003608502\t4-5\t1', 'chr1\t1480767\t1481036\tENSG00000160072\tENSE00003474888\t4-5\t1'], does_not_raise()),
        (["ENSG00000160072\tfull\t0"], ['chr1\t1471765\t1497848\tENSG00000160072\tfull\tfull\t1'], does_not_raise()),
        (["ENSG00000160072\tENSE00003889014,\t0"], ['chr1\t1471765\t1472089\tENSG00000160072\tENSE00003889014\texon\t1'], does_not_raise()),
        #(['1','.','1I1','A'], 100, 'ERROR!!!!!!!!!!!!!!!!!!!!!', pytest.raises(ValueError))
    ]
)
def test_get_regions_for_genes(input_genes_list, answer, ext, logger):
    annot_for_genes = {'ENSG00000160072': ['1', ['1471765', '1497848'], [['1471765', '1472089'], ['1477274', '1477350'], ['1478644', '1478745'], ['1479049', '1479108'], ['1480867', '1480936'], ['1482138', '1482303'], ['1482545', '1482614'], ['1485016', '1485171'], ['1485782', '1485838'], ['1486110', '1486235'], ['1486544', '1486668'], ['1487863', '1487914'], ['1489204', '1489274'], ['1490257', '1490424'], ['1490563', '1490671'], ['1495485', '1497848']], [['1471765', '1471884'], ['1495818', '1497848']], ['ENSE00003889014', 'ENSE00003467707', 'ENSE00003569130', 'ENSE00003608502', 'ENSE00003474888', 'ENSE00003654064', 'ENSE00003510521', 'ENSE00003459370', 'ENSE00003655926', 'ENSE00003594545', 'ENSE00003892109', 'ENSE00003689846', 'ENSE00003670332', 'ENSE00003505365', 'ENSE00003497242', 'ENSE00003892597']], 'test_g': ['1', ['10', '500'], [['10', '90'], ['100', '200'], ['400', '500']], [['1', '9'], ['501', '600']], ['exo1', 'exo2', 'exo3']], 'ENSG00000277475': ['KI270713.1', ['31968', '32528'], [['31698', '32528']], [['32373', '32528'], ['31698', '31841']], ['ENSE00003739404']], 'ENSG00000275405': ['KI270713.1', ['21861', '22024'], [['21861', '22024']], [], ['ENSE00003747303']]}
    gene_name_dict = {'ATAD3B': 'ENSG00000160072', 'test_name': 'test_g', 'GENENAME': 'ENSG00000277475', 'U1': 'ENSG00000275405'}
    
    with ext:
        gene_coords = get_regions_for_genes(input_genes_list, annot_for_genes, gene_name_dict,
                    'exon', 50, False, '', logger, 'GRCh38')
        assert gene_coords == answer