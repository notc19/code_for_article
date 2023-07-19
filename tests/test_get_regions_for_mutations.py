import pytest
from contextlib import nullcontext as does_not_raise

import logging

from OligoDesigner.get_regions_for_mutations import *
from OligoDesigner.constants import genome_versions_chromosome_dict

@pytest.mark.parametrize(
    "vcf, flank, answer, exp",
    [
        (['1','.','111','A'], 100, '1\t11\t211\tA\t.', does_not_raise()),
        (['1','.','111','A', '', 'PASS','AC=5'], 100, '1\t11\t211\tA\t.', does_not_raise()),
        (['1','.','1I1','A'], 100, 'ERROR!!!!!!!!!!!!!!!!!!!!!', pytest.raises(ValueError))
    ]
)
def test_get_coords_for_vcf(vcf, flank, answer, exp):
    with exp:
        assert get_coords_for_vcf(vcf,flank) == answer

@pytest.mark.parametrize(
    "hgvsg, flank, answer, exp",
    [
        (['chr1','123215A>T'],   100, 'chr1\t123115\t123315\tA\tchr1:123215A>T',   does_not_raise()),
        (['chr1','g.123215A>T'], 100, 'chr1\t123115\t123315\tA\tchr1:g.123215A>T', does_not_raise()),
        (['chr1','g.123+215A>T'],100, 'ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!', pytest.raises(ValueError))
    ]
)
def test_get_coords_for_hgvs_g(hgvsg, flank, answer, exp):
    with exp:
        assert get_coords_for_hgvs_g(hgvsg, flank) == answer

@pytest.mark.parametrize(
        "hgvsc, flank, annot, logger, answer",
    [
        (['test','c.47A>T'], 10, {}, "", '1\t47\t67\tA\ttest:c.47A>T'),
        (['test','c.47-5A>T'], 10, {}, "", '1\t42\t62\tA\ttest:c.47-5A>T'),
        (['test','c.47_48delA'], 10, {}, "", '1\t47\t67\tdelA\ttest:c.47_48delA'),
        (['test','c.47_50insATGC'], 10, {}, "", '1\t47\t67\tinsATGC\ttest:c.47_50insATGC'),
        (['test','c.47-2_50-2insATGC'], 10, {}, "", '1\t45\t65\tinsATGC\ttest:c.47-2_50-2insATGC'),
        (['test','c.-2A>G'], 5, {}, "", '1\t3\t13\tA\ttest:c.-2A>G'),
        (['test','c.-2_3insAT'], 5, {}, "", '1\t3\t13\tinsAT\ttest:c.-2_3insAT'),
        (['test','c.-2_3delAT'], 5, {}, "", '1\t3\t13\tdelAT\ttest:c.-2_3delAT'),
        (['test','c.-2-1A>G'], 5, {}, "", '1\t2\t12\tA\ttest:c.-2-1A>G'),
        (['test','c.-2+1A>G'], 5, {}, "", '1\t4\t14\tA\ttest:c.-2+1A>G'),
        (['test','c.-2-1_-2-2insG'], 5, {}, "", '1\t2\t12\tinsG\ttest:c.-2-1_-2-2insG'),
        (['test','c.*10A>G'], 50, {}, "", '1\t460\t560\tA\ttest:c.*10A>G'),
        (['test','c.*10_11insAT'], 50, {}, "", '1\t460\t560\tinsAT\ttest:c.*10_11insAT'),
        (['test','c.*10-2A>G'], 50, {}, "", '1\t458\t558\tA\ttest:c.*10-2A>G'),
        (['test','c.*10+2A>G'], 50, {}, "", '1\t462\t562\tA\ttest:c.*10+2A>G'),
        (['test','c.*10-2_10-1delAT'], 50, {}, "", '1\t458\t558\tdelAT\ttest:c.*10-2_10-1delAT'),
        (['test','c.?10-2_10-1delAT'], 50, {}, "", None),
        (['NotFound','c.47A>T'], 10, {}, "", None)
    ]
)
def test_get_coords_for_hgvs_c(hgvsc, flank, annot, logger, answer):
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)
    annot = {'test':['1', [[10,90],[100,200],[400,500]]]}
    
    assert get_coords_for_hgvs_c(hgvsc, flank, annot, logger) == answer

@pytest.mark.parametrize(
    "coord_list, genome_version, logger, answer",
    [
        (['chr1\t111\t211\tA\t.\tchr1:111A>T'], 'GRCh38.p14', '', ['chr1\t111\t211\tA\t.\tchr1']),
        (['1\t111\t211\tA\t.\tchr1:111A>T'],    'GRCh38.p14', '', ['chr1\t111\t211\tA\t.\t1']),
        (['NC_000001.11\t111\t211\tA\t.\tchr1:111A>T'], 'GRCh38.p14', '', ['chr1\t111\t211\tA\t.\tNC_000001.11']),
        (['NC_000001.10\t111\t211\tA\t.\tchr1:111A>T'], 'GRCh37.p13', '', ['chr1\t111\t211\tA\t.\tNC_000001.10']),
        (['NC_000001.9\t111\t211\tA\t.\tchr1:111A>T'], 'GRCh37.p13', '', []),
        (['Z\t111\t211\tA\t.\tchr1:111A>T'],    'GRCh38.p14', '', []),
        (['chr1\t111\t211\tA'], 'GRC14', '', []),
    ]
)
def test_fix_chromosome_name_and_check_coords(coord_list, genome_version, logger, answer):
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)

    assert fix_chromosome_name_and_check_coords(coord_list, genome_version, logger) == answer

@pytest.mark.parametrize(
    "muts_list, annot, flank_length, write_intermediate, tmpdir, logger, genome_version, answer, exp",
    [
        (['chr1\t.\t211\tA\t.\tchr1:111A>T'], {}, 10, False, "", "", 'GRCh38.p14', ['chr1\t201\t221\tA\t.\tchr1'], does_not_raise()),
        (['chr1:211A>T'], {}, 10, False, "", "", 'GRCh38.p14', ['chr1\t201\t221\tA\tchr1:211A>T\tchr1'], does_not_raise()),
        (['chr1:g.211A>T'], {}, 10, False, "", "", 'GRCh38.p14', ['chr1\t201\t221\tA\tchr1:g.211A>T\tchr1'], does_not_raise()),
        (['test:c.47-2_50-2insATGC'], {}, 10, False, "", "", 'GRCh38.p14', ['chr1\t45\t65\tinsATGC\ttest:c.47-2_50-2insATGC\t1'], does_not_raise()),
        (['test:p.47-2_50-2insATGC'], {}, 10, False, "", "", 'GRCh38.p14', [], pytest.raises(SystemExit)),
        (['Z:211A>T'], {}, 10, False, "", "", 'GRCh38.p14', [], pytest.raises(SystemExit)),
        (['chr1:g.211A>T'], {}, 1000, False, "", "", 'GRCh38.p14', [], pytest.raises(SystemExit))
    ]
)
def test_get_regions_for_muts(muts_list, annot, flank_length, write_intermediate, tmpdir, logger, genome_version, answer, exp):
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)

    annot = {'test':['1', [[10,90],[100,200],[400,500]]]}

    with exp:
        assert get_regions_for_muts(muts_list, annot, flank_length, write_intermediate, tmpdir, logger, genome_version) \
            == answer