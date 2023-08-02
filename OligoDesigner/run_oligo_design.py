import multiprocessing
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
    def __init__(self, idi, seq, gc, entr, mt, repeat, masked):
        self.idi = idi
        self.seq = str(seq)
        self.gc = float(gc)
        self.entr = float(entr)
        self.mt = float(mt)
        self.repeat = float(repeat)
        self.masked = float(masked)
        self.blastresult = 1

        try:
            self.qual = round((self.entr - self.repeat - self.masked) / (abs(0.65-self.gc) + abs(70-self.mt)/100), 3)
        except ZeroDivisionError:
            self.qual = round(self.entr / 0.01, 3)

        if self.gc < 0.2 or self.mt < 40 or self.repeat > 0.20 or self.masked > 0.25:
            self.status = 'bad'
        elif 0.50 <= self.gc <= 0.60 and 65 <= self.mt <=80 and self.entr > 1.9:
            self.status = 'perfect'
        else:
            self.status = 'best'

    def __str__(self):
        return f'{self.idi} {self.gc} {self.entr} {self.mt} {self.repeat} {self.masked} {self.blastresult} {self.qual} {self.status} {self.seq}'

    def __repr__(self):
        return f'{self.idi} {self.gc} {self.entr} {self.mt} {self.repeat} {self.masked} {self.blastresult} {self.qual} {self.status} {self.seq[:]}...\n' #


def design_oligos_snp(sequence: SeqRecord) -> Oligonucleotide:
    """
    Эта функция выполняет дизайн олигонуклеотида для заданного региона со снипом. Перебираются все возможные олигонуклеотиды
    длиной length от первого до последнего в последовательности, и для каджого считаются параметры. Затем выбирается один 
    (или первый попавшийся идеальный, или тот, у которого максимальный суммарный параметр qual).
    Подразумевается, что снип находится в центре последовательности и даже крайние олигонуклеотиды будут включать этот снип
    (length > len(sequence/2))
    :param sequence: Объект SeqRecord с последовательностью региона
    :return: Объект Oligonucleotide
    """
    best_qual = 0
    max_id, max_cov, max_eval = params.blast_filtering_params.split(',')[0], params.blast_filtering_params.split(',')[1], \
                                params.blast_filtering_params.split(',')[2]

    for i in range(0,len(sequence)-params.length+1):
        subseq = sequence.seq[i:i+params.length]

        olig = Oligonucleotide(sequence.id, subseq, *st.get_olig_stats(subseq))

        if params.blast_strategy == 'filter':
            if olig.status == 'perfect' or olig.qual > best_qual:
                blast_result = st.blastn(params.blast_binary, olig.seq, params.blast_db, 10, 1)
                blast_numhits = st.read_blast_results(blast_result, params.length, max_id, max_cov, max_eval)
                olig.blastresult = blast_numhits
        if olig.status == 'perfect':
            if olig.blastresult == 1:
                best_olig = olig
                break
        elif olig.qual >= best_qual:
            best_olig = olig

    if params.blast_strategy == 'analyze':
        blast_result = st.blastn(params.blast_binary, best_olig.seq, params.blast_db, 10, 1)
        blast_numhits = st.read_blast_results(blast_result, params.length, max_id, max_cov, max_eval)
        best_olig.blastresult = blast_numhits

    return best_olig

def design_oligos(sequence: SeqRecord) -> Oligonucleotide:
    """
    Эта функция выполняет дизайн олигонуклеотида для заданного региона.
    :param sequence: последовательость (объект SeqRecord)
    :return: Олигонуклеотид (объект Oligonucleotide)
    """
    olig_list = []

    # str with regions sequence from SeqRecord object
    seq = str(sequence.seq)

    # create first olig - just first {length} nucleotides from start
    subseq = seq[:params.length]
    olig = Oligonucleotide(sequence.id, subseq, *st.get_olig_stats(subseq))
    olig_list.append(olig)

    while True:
        previous_end = seq.index(olig.seq) + len(olig.seq)-1
        start = previous_end-params.max_i
        end = previous_end-params.min_i if previous_end+params.min_i < len(seq) else len(seq)
        new_subseq = seq[start:end+params.length]

        if len(new_subseq) >= params.length:
            olig = design_oligos_snp(SeqRecord(Seq(new_subseq), id = sequence.id))
            olig_list.append(olig)
        else:
            olig = Oligonucleotide(sequence.id, seq[-params.length-1:], *st.get_olig_stats(seq[-params.length-1:]))
            olig_list.append(olig)
            break

        if end == len(seq):
            break

    return olig_list

def run_oligo_design(seq_list:[SeqRecord], threads:int, logger:logging, param:argparse) -> [Oligonucleotide]:
    """
    Эта функция запускает дизайн олигов для последовательностей снипов и прочих регионов.
    :param seq_list: Список последовательностей (SeqRecord)
    :param threads: Количество ядер
    :param logger: Объект класса logging
    :param param: Объект argparse с параметрами
    :return: Список олигонуклеотидов (объект Oligonucleotide)
    """
    logger.info('Start designing oligonucleotides!')
    global params
    params = param

    snp_seq_list = [seq for seq in seq_list if seq.name == 'snp' or 'snp' in seq.description]
    other_seq_list = [seq for seq in seq_list if not seq.name == 'snp' and not 'snp' in seq.description]
    
    logger.info(f'There were found {len(snp_seq_list)} sequences with SNP and {len(other_seq_list)} other regions')

    output_oligos = []
    lock = multiprocessing.Lock()
    if snp_seq_list:
        logger.info(f'Start designin oligonucleotides for SNPs')
        try:
            with multiprocessing.Pool(threads, initargs=(lock,),) as pool:
                snp_oligos = pool.map(design_oligos_snp, snp_seq_list)
            output_oligos.extend(snp_oligos)
            logger.info(f'Oligonucleotides for SNPs were designed!')
        except Exception as e:
            logger.error(e)

    if other_seq_list:
        logger.info(f'Start designin oligonucleotides for all regions')
        try:
            other_oligos = []
            threads = len(other_seq_list) if len(other_seq_list) < threads else threads
            with multiprocessing.Pool(threads, initargs=(lock,),) as pool:
                other_oligos_tmp = pool.map(design_oligos, other_seq_list)
            for oligos in other_oligos_tmp:
                other_oligos.extend(oligos)
            output_oligos.extend(other_oligos)
            logger.info(f'Oligonucleotides were designed!')
        except Exception as e:
            logger.error(e)

    return output_oligos

if __name__ == '__main__':
    inpt = '/home/a_vikhorev/OligoDesigner/tests/test_reuslts/FIXED_clustered_sequences.fasta'
    #inpt = '/home/a_vikhorev/OligoDesigner/tests/test_reuslts/TEST_sequences.fasta'
    seq_list = list(SeqIO.parse(inpt, "fasta"))

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    global params
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('--max_i', help='Maximum intersection of oligos in region',
                             type=int, required=False, default=40)
    args_parser.add_argument('--min_i', help='Minimum intersection of oligos in region',
                             type=int, required=False, default=-20)
    args_parser.add_argument('-l','--length', help='Length of the designed olignucleotide bait',
                             type=int, required=False, default=100)
    args_parser.add_argument('--threads', help='Number of threads',
                             type=int, required=False, default=4)
    args_parser.add_argument('--blast_strategy', help='Для каких целей использовать blast. Варианты: filter/analyze',
                             type=str, required=False, default='analyzee')
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

    print(len(oligos))
    print(oligos[-5:])