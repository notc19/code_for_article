import multiprocessing
from itertools import product
import logging
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sequence_tools as st


class Oligonucleotide():
    """
    Этот класс содержит объект олигонуклеотида
    """
    def __init__(self, id, seq, gc, entr, mt, repeat, masked):
        self.id = id
        self.seq = seq
        self.gc = gc
        self.entr = entr
        self.tm = mt
        self.repeat = repeat
        self.masked = masked
        self.blastresult = 1

        try:
            self.qual = round((entr - repeat - masked) / (abs(0.65-gc) + abs(70-mt)/100), 3)
        except ZeroDivisionError:
            self.qual = round(entr / 0.01, 3)

        if gc < 0.2 or mt < 40 or repeat > 0.20 or masked > 0.25:
            self.status = 'bad'
        elif 0.50 <= gc <= 0.60 and 65 <= mt <=80 and entr > 1.9:
            self.status = 'perfect'
        else:
            self.status = 'best'

    def __str__(self):
        return f'{self.id} {self.gc} {self.entr} {self.tm} {self.repeat} {self.masked} {self.blastresult} {self.qual} {self.status} {self.seq}'

    def __repr__(self):
        return f'{self.id} {self.gc} {self.entr} {self.tm} {self.repeat} {self.masked} {self.blastresult} {self.qual} {self.status} {self.seq[:]}...\n' #


def design_oligos(sequence: SeqRecord, length: int) -> Oligonucleotide:
    pass

def design_oligos_snp(sequence: SeqRecord, length: int) -> Oligonucleotide:
    """
    Эта функция выполняет дизайн олигонуклеотида для заданного региона со снипом. Перебираются все возможные олигонуклеотиды
    длиной length от первого до последнего в последовательности, и для каджого считаются параметры. Затем выбирается один 
    (или первый попавшийся идеальный, или тот, у которого максимальный суммарный параметр qual).
    Подразумевается, что снип находится в центре последовательности и даже крайние олигонуклеотиды будут включать этот снип
    (length > len(sequence/2))
    :param sequence: Объект SeqRecord с последовательностью региона
    :param length: длина олигонуклеотида
    :return: Объект Oligonucleotide
    """
    best_qual = 0
    for i in range(0,len(sequence)-length):
        subseq = sequence.seq[i:i+length]

        gc = st.gc_content(subseq)
        entropy = st.estimate_shannon_entropy(subseq)
        mt = st.calculate_melting_temp(subseq)
        longest_repeat = st.calculate_longest_repeat_percent(subseq)
        longest_masked = st.calculate_repeat_masked_percent(subseq)

        olig = Oligonucleotide(sequence.id, subseq, gc, entropy, mt, longest_repeat, longest_masked)

        if blast_strategy == 'filter':
            if olig.status == 'perfect' or olig.qual > best_qual:
                blast_result = st.blastn(blast_binary, olig.seq, blast_db, 10, 1)
                blast_numhits = st.read_blast_results(blast_result, length, max_id, max_cov, max_eval)
                olig.blastresult = blast_numhits
        if olig.status == 'perfect':
            if olig.blastresult == 1:
                best_olig = olig
                break
        elif olig.qual >= best_qual:
            best_olig = olig

    if blast_strategy == 'analyze':
        blast_result = st.blastn(blast_binary, best_olig.seq, blast_db, 10, 1)
        blast_numhits = st.read_blast_results(blast_result, length, max_id, max_cov, max_eval)
        best_olig.blastresult = blast_numhits

    return best_olig


def run_oligo_design(seq_list:[SeqRecord], threads:int, logger:logging, params):
    """
    
    """
    global length, blast_binary, blast_db, blast_strategy, max_id, max_cov, max_eval
    length, blast_binary, blast_db, blast_strategy = params.bait_length, params.blast_binary, params.blast_db, params.blast_strategy
    max_id, max_cov, max_eval = params.blast_filtering_params.split(',')
    logger.info('Start designing oligonucleotides!')

    snp_seq_list = [seq for seq in seq_list if seq.name == 'snp' or 'snp' in seq.description]
    other_seq_list = [seq for seq in seq_list if not seq.name == 'snp' and not 'snp' in seq.description]
    
    logger.info(f'There were found {len(snp_seq_list)} sequences with SNP and {len(other_seq_list)} other regions')

    output_oligos = []
    lock = multiprocessing.Lock()
    if snp_seq_list:
        logger.info(f'Start designin oligonucleotides for SNPs')
        try:
            with multiprocessing.Pool(threads, initargs=(lock,),) as pool:
                snp_oligos = pool.starmap(design_oligos_snp, product(snp_seq_list, [length]))
            output_oligos.extend(snp_oligos)
            logger.info(f'Oligonucleotides for SNPs were designed!')
        except Exception as e:
            logger.error(e)

    if other_seq_list:
        logger.info(f'Start designin oligonucleotides for all regions')
        try:
            with multiprocessing.Pool(threads, initargs=(lock,),) as pool:
                other_oligos = pool.starmap(design_oligos_snp, product(other_seq_list, [length]))
            output_oligos.extend(other_oligos)
            logger.info(f'Oligonucleotides were designed!')
        except Exception as e:
            logger.error(e)

    return output_oligos

if __name__ == '__main__':
    inpt = '/home/a_vikhorev/OligoDesigner/tests/test_reuslts/FIXED_clustered_sequences.fasta'
    seq_list = list(SeqIO.parse(inpt, "fasta"))

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('-l','--bait_length', help='Length of the designed olignucleotide bait',
                             type=int, required=False, default=100)
    args_parser.add_argument('--threads', help='Number of threads',
                             type=int, required=False, default=4)
    args_parser.add_argument('--blast_strategy', help='Для каких целей использовать blast. Варианты: filter/analyze',
                             type=str, required=False, default='analyze')
    args_parser.add_argument('--blast_binary', help='Путь до исполняемого файла blastn',
                             type=str, required=False, default='/home/a_vikhorev/soft/miniconda3/envs/mrbait/bin/blastn')
    args_parser.add_argument('--blast_db', help='Путь до базы данны blast',
                             type=str, required=False, default='/home/a_vikhorev/Data/DataBases/hg38_blast/hg38')
    args_parser.add_argument('--blast_filtering_params', help='Параметры для фильтрации результатов бласта (identity, coverage, степень e-value)',
                             type=str, required=False, default='0.85,0.85,5')
    params = args_parser.parse_args()

    binary='/home/a_vikhorev/soft/miniconda3/envs/mrbait/bin/blastn'
    db='/home/a_vikhorev/Data/DataBases/hg38_blast/hg38'

    oligos = run_oligo_design(seq_list, 14, logger, params)

    print(oligos)