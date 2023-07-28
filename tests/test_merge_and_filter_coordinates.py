import logging

from OligoDesigner.merge_and_filter_coordinates import *

def test_merge_and_filter_coordinates():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    input_bed_list = ['chr1\t100\t200\tENSG00000160072\tENSE00003889014\fghfgh\t1',
                      'chr1\t100\t200\tENSG00000160072\tENSE00003889014\texon\t1',
                      'chr1\t100\t500\tENSG00000160072\tENSE00003889014\texon\t1',
                      'chr1\t100\t250\tENSG00000160072\tENSE00003889014\texon\t1',
                      'chr1\t120\t300\tENSG00000160072\tENSE00003889014\texon\t1',
                      'chr1\t200\t800\tENSG00000160072\tfull\tfull\t1',
                      'chr1\t400\t450\tENSG00000160072\tfull\tfull\t1',
                      'chr1\t700\t900\tENSG00000160072\tfull\tfull\t1',
                      'chr1\t7100\t9400\tENSG00000160072\tfull\tfull\t1',
                      'chr1\t900\t1000\tENSG00000160072\tfull\tfull\t1']
    
    answer = ['chr1\t100\t1000\tENSG00000160072\tfull,full,full,ENSE00003889014\tfull,full,full,exon\t1', 
              'chr1\t7100\t9400\tENSG00000160072\tfull\tfull\t1']
    
    assert process_coordinates(input_bed_list, False, '', logger) == answer