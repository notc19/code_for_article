import pytest
from contextlib import nullcontext as does_not_raise

import logging

from OligoDesigner.get_regions_for_mutations import *
from OligoDesigner.utils import *
from OligoDesigner.constants import genome_versions_chromosome_dict

def test_prepare_annot_for_muts():
    test_gtf = './tests/test_annot.gtf'

    annot = prepare_annot_for_muts(['Chr2:c.123G>C'], test_gtf, '')

    answer = {'ENST00000673477.1': ['1', [['1471765', '1472089'], ['1477274', '1477350'], ['1478644', '1478745'], ['1479049', '1479108'], ['1480867', '1480936'], ['1482138', '1482303'], ['1482545', '1482614'], ['1485016', '1485171'], ['1485782', '1485838'], ['1486110', '1486235'], ['1486544', '1486668'], ['1487863', '1487914'], ['1489204', '1489274'], ['1490257', '1490424'], ['1490563', '1490671'], ['1495485', '1497848']]], 'test.7': ['1', [['10', '90'], ['100', '200'], ['400', '500']]], 'ENST00000612315.1': ['KI270713.1', [['31698', '32528']]], 'ENST00000619109.1': ['KI270713.1', [['21861', '22024']]]}

    assert annot == answer

def test_prepare_annot_for_muts_full():

    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)

    test_gtf = './tests/test_annot.gtf'

    annot = prepare_annot_for_muts(['Chr2:c.123G>C'], test_gtf, '')
    muts_list = ['test.7:c.47-2_50-2insATGC']

    answer = get_regions_for_muts(muts_list, annot, 10, False, "", logger, 'GRCh38.p14')

    assert answer == ['chr1\t45\t65\tinsATGC\ttest.7:c.47-2_50-2insATGC\t1']
