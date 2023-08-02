from pathlib import Path
import logging
import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from run_oligo_design import Oligonucleotide

def cluster_oligos(tmpdir:Path, outdir:Path, logger:logging, vsearch_binary:str, max_id:float) -> [Oligonucleotide]:
    """
    Эта функция выполняет кластеризацию поданных на вход олигов.
    :param tmpdir: директория для промежуточных результатов
    :param outdir: Директория, в которую будут сохраняться результаты
    :param logger: Объект класса logging
    :param vsearch_binary: Путь до исполняемого файла vsearch
    :param max_id: Максимально допустимая доля совпадающих нуклеотидов в последовательности
    :return: Список олигонуклеотидов (объект SeqRecords)
    """
    vsearch_binary = vsearch_binary if vsearch_binary else 'vsearch'
    vsearch_options = [vsearch_binary,
                        '--cluster_fast', f'{tmpdir}/OUTPUT_Oligonucleotides_RAW.fasta',
                        '--centroids', f'{outdir}/OUTPUT_Oligonucleotides_FIXED.fasta',
                        '--clusterout_id',
                        '--id', f'{max_id}',
                        '--iddef', '0',
                        '--threads', '1',
                        '--maxseqlength', '500000',
                        '--relabel_keep'
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
    
    fixed_fasta_list = list(SeqIO.parse(f'{outdir}/OUTPUT_Oligonucleotides_FIXED.fasta', "fasta"))
    
    return fixed_fasta_list

def trim_oligos(olig_list:list) -> [Oligonucleotide]:
    """
    Эта функция выполняет тримминг олигонуклеотидов. С конца олига будут удаляться нуклеотиды до достижения
    минимально допустимой длины, если это удаление будет приводить к увеличению качества олигонуклеотида
    :param olig_list: список олигонуклеотидов (объект Oligonucleotide)
    """
    
    return olig_list

def postprocess_oligos(tmpdir: Path, outdir:Path, vsearch_binary:str, max_id:float, logger:logging) -> [Oligonucleotide]:
    """
    Эта функция выполняет построцессинг олигонуклеотидов (фильтрацию по кластеризации и тримминг)
    :param outdir: Директория с результатами
    :param vsearch_binary: Путь до исполняемого файла vsearch
    :param max_id: Максимально допустимая доля совпадающих нуклеотидов в последовательности
    :param logger: Объект класса logging
    :return: Список олигонуклеотидов (объект Oligonucleotide)
    """

    # First step - cluster oligonucleotides using vsearch
    clustered_oligs = cluster_oligos(tmpdir, outdir, logger, vsearch_binary, max_id)
    seq_records_dict = {record.id:record.description \
                        for record in list(SeqIO.parse(f'{tmpdir}/OUTPUT_Oligonucleotides_RAW.fasta', "fasta"))}
    for record in clustered_oligs:
        record.description = seq_records_dict[record.id.split(';')[0]]
    olig_list = [Oligonucleotide(o.description.split(' ')[1], str(o.seq), *o.description.split(' ')[2:7]) for o in clustered_oligs]

    # Second step - trim oligonucleotides
    olig_list = trim_oligos(olig_list)

    return olig_list


if __name__ == '__main__':
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    input_fasta_list = []
    tmpdir=Path('./tests/test_reuslts')
    #cluster_oligos(tmpdir, logger, '/home/a_vikhorev/soft/miniconda3/envs/mrbait/bin/vsearch', 0.80)

    postprocess_oligos(input_fasta_list, tmpdir, '/home/a_vikhorev/soft/miniconda3/envs/mrbait/bin/vsearch', 0.80)