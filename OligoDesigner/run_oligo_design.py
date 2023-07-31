import multiprocessing
from itertools import product
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sequence_tools as st


class Oligonucleotide():
    """
    Этот класс содержит объект олигонуклеотида
    """
    def __init__(self, id, seq, gc, entr, mt):
        self.id = id
        self.seq = seq
        self.gc = gc
        self.entr = entr
        self.tm = mt

        try:
            self.qual = round(entr / (abs(0.65-gc) + abs(70-mt)/100), 3)
        except ZeroDivisionError:
            self.qual = round(entr / 0.01, 3)
        self.status = 'perfect' if 0.55 <= gc <= 0.65 and 65 <= mt <=75 and entr > 1.9 else ''

    def __str__(self):
        return f'{self.id} {self.gc} {self.entr} {self.tm} {self.qual} {self.status} {self.seq}'

    def __repr__(self):
        return f'{self.id} {self.gc} {self.entr} {self.tm} {self.qual} {self.status} {self.seq[:10]}...\n'


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
    #print(sequence)
    best_qual = 0
    for i in range(0,len(sequence)-length):
        subseq = sequence.seq[i:i+length]
        gc = st.gc_content(subseq)
        entropy = st.estimate_shannon_entropy(subseq)
        mt = st.calculate_melting_temp(subseq)
        olig = Oligonucleotide(sequence.id, subseq, gc, entropy, mt)
        if olig.status == 'perfect':
            return olig
        elif olig.qual >= best_qual:
            best_olig = olig
    best_olig.status = 'best'

    return best_olig


def run_oligo_design(seq_list, threads, length, logger):
    """
    
    """
    logger.info('Start designing oligonucleotides!')

    snp_seq_list = [seq for seq in seq_list if seq.name == 'snp' or 'snp' in seq.description]
    other_seq_list = [seq for seq in seq_list if not seq.name == 'snp' and not 'snp' in seq.description]
    
    logger.info(f'There were found {len(snp_seq_list)} regions with SNP and {len(other_seq_list)} other regions')

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

    oligos = run_oligo_design(seq_list, 4, 100, logger)

    print([o.split(' ') for o in list(map(str,oligos))])