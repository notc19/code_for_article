import logging
from pathlib import Path
import os
import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def process_sequences(input_fasta_list:[SeqRecord], tmpdir:Path, logger:logging, vsearch_binary:str, max_id:float) -> list:
    """
    Эта функция выполняет кластеризацию поданных на вход последовательностей при помощи VSEARCH cluster_fast.
    Затем кластеризованные последовательности считываются при помощи SeqIO, записываются в список и этот список возвращается
    :param input_fasta_list: Список с последовательностями в формате fasta
    :param tmpdir: Объект Path с путем до промежуточной директории
    :param logger: Объект класса logging
    :return: Список с fasta-последовательсотями (объект SeqRecord)
    """
    
    if not os.path.exists(tmpdir/'INPUT_regions_sequences.fasta'):
        with open(tmpdir/'INPUT_regions_sequences.fasta', 'w') as outfile:
             SeqIO.write(input_fasta_list, outfile, "fasta")
    vsearch_binary = vsearch_binary if vsearch_binary else 'vsearch'
    vsearch_options = [vsearch_binary,
                        '--cluster_fast', f'{tmpdir}/INPUT_regions_sequences.fasta',
                        '--centroids', f'{tmpdir}/FIXED_clustered_sequences.fasta',
                        '--clusterout_id',
                        '--id', f'{max_id}',
                        '--iddef', '0',
                        '--threads', '1'
                       ]
    
    vsearch_command = ' '.join(vsearch_options)
    logger.info(f'Running VSEARCH with the following command: {vsearch_command}')
    logger.warning(f'WARNING! Sequences with length >50000 bp will be removed! This is vsearch internal filter')
    output = subprocess.run(vsearch_options, capture_output=True)
    stderr = str(output.stderr)
    if output.returncode:
        logger.error('ERROR in VSEARCH! Error message is:')
        logger.error(stderr)
        exit(1)
    else:
        logger.info(f'VSEARCH clustering is done!')
        logger.info(stderr)
    
    fixed_fasta_list = list(SeqIO.parse(f'{tmpdir}/FIXED_clustered_sequences.fasta', "fasta"))

    logger.info(f'There were {len(input_fasta_list)} sequence in input list, but after clustering {len(fixed_fasta_list)} was remained!')

    return fixed_fasta_list

if __name__ == '__main__':
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    input_fasta_list = []
    tmpdir=Path('./tests/test_reuslts')
    process_sequences(input_fasta_list, tmpdir, logger)