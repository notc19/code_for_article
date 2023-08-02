import math
import collections
import subprocess
from subprocess import Popen, PIPE, CalledProcessError
import os

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast.Applications import NcbiblastpCommandline

def gc_content(seq:str) -> float:
    """
    Эта функция возвращает ГЦ-состав последовательности
    :param seq: последовательность олигонуклеотида
    :return: ГЦ-состав (доля нуклеотидов Г и Ц) олигонуклеотида
    """
    return (seq.count('G') + seq.count('C'))/len(seq)

def estimate_shannon_entropy(seq:str) -> float:
    """
    Эта функция вычисляет энтропию Шеннона для последовательности.
    Чем выше значение энтропии, тем более разнообразна структура последовательности. В идеале должна быть как можно ближе к 2.
    :param seq: последовательность ДНК
    :return: Значение энтропии (число от 0 до 2)
    """
    m = len(seq)
    seq = seq.upper()
    bases = collections.Counter([tmp_base for tmp_base in seq])

    shannon_entropy_value = 0
    for base in bases:
        n_i = bases[base]
        p_i = n_i / float(m)
        entropy_i = p_i * (math.log(p_i, 2))
        shannon_entropy_value += entropy_i
 
    return round(shannon_entropy_value * -1, 3)

def calculate_melting_temp(seq:str) -> float:
    """
    Эта функция вычисляет температуру плавления олигонуклеотида
    :param seq: последовательность олигонуклеотида
    :return: температура плавления олигонуклеотида
    """
    return round(mt.Tm_GC(seq, strict=False), 3)

def calculate_longest_repeat_percent(seq:str) -> float:
    """
    Эта функция считает отношение длины самого длинного однонуклеотидного повтора в последовательности к длине последовательности.
    :param seq: Последовательность олигонуклеотида (строка)
    :return: Отношение длины самого длинного повтора к длине последовательности
    """
    i=-1
    counter = 0
    while i < len(seq)-1:
        k = 1
        while (i+k < len(seq) - 1) and (seq[i] == seq[i+k]):
            k += 1
            k += 1 if (i+k == len(seq) - 1) and (seq[i] == seq[i+k]) else 0
        counter = k if k > counter else counter
        i+=k

    return round(counter / len(seq), 3)

def calculate_repeat_masked_percent(seq:str) -> float:
    """
    Эта функция считает процент последовательности, который был маскирован RepeatMasker-ом (маленькие буквы (sm) или N (rm)).
    :param seq: Последовательность олигонуклеотида
    :return: Процент маскированого повтора в последовательности (или 0.05, если повторов нет)
    """
    i=-1
    counter = 0
    while i < len(seq)-1:
        k = 0
        if seq[i+k].islower():
            while (i+k < len(seq) - 1) and (seq[i+k].islower()):
                k += 1
                k += 1 if (i+k == len(seq) - 1) and (seq[i+k].islower()) else 0
        elif seq[i+k] == 'N':
            while (i+k < len(seq) - 1) and (seq[i+k] == 'N'):
                k += 1
                k += 1 if (i+k == len(seq) - 1) and (seq[i+k] == 'N') else 0

        counter = k if k > counter else counter

        i+=k if k > 0 else 1

    return round(counter / len(seq), 3) if counter else 0.05

def blastn(binary:str, seq:str, blastdb:str, maxhits:int, threads:int):
    """
    Эта функция выполняет бласт с заданными параметрами для поданной на вход последовательности и возвращает
    список результатов бласта в формате outfmt6
    :param binary: путь до binary файла бласта
    :param seq: последовательность олигонуклеотида (строка)
    :param blastdb: путь до базы данных blast
    :param maxhits: сколько максимум хитов для одной последовательности выписывать
    :param threads: сколько ядер использовать
    :return: список результатов бласта в формате outfmt 6
    """
    binary = binary if binary else 'blastn'
    blastn_cli = [
        binary,
		"-db", str(blastdb),
		"-evalue", str(1),
		"-outfmt", '6',
		"-num_threads", str(threads),
		"-max_target_seqs", str(maxhits),
		"-max_hsps", str(1)
       ]

    proc = Popen(blastn_cli, stdout=PIPE, stdin=PIPE, env={'PATH': os.getenv('PATH')})
    t = proc.communicate(str.encode(str(seq)))#[0]
    blastresult = t[0].decode().split('\n')
    
    return blastresult

def read_blast_results(blastresult:list, length:int, max_id:float, max_cov:float, max_eval:int) -> int:
    """
    Эта функция парсит результаты бласта и считает, сколько нашлось последовательностей с заданными
    значениями pident, length, evalue
    :param blastresult: Список, содержащий результаты бласта в формате outfmt 6
    :return: Число, означающие количество хитов с заданными параметрами
    """
    mismatch_counter = 0
    for result in blastresult:
        if result:
            pident, length, evalue = float(result.split()[2]), int(result.split()[3]), result.split()[-2]
            if pident > float(max_id) and length > float(max_cov)*length and '-' in evalue and int(evalue.split('-')[1]) > int(max_eval):
                mismatch_counter += 1
    
    return mismatch_counter

def get_olig_stats(subseq):
    """
    
    """
    gc = gc_content(subseq)
    entropy = estimate_shannon_entropy(subseq)
    mt = calculate_melting_temp(subseq)
    longest_repeat = calculate_longest_repeat_percent(subseq)
    longest_masked = calculate_repeat_masked_percent(subseq)

    return gc, entropy, mt, longest_repeat, longest_masked

if __name__ == '__main__':
    #seq = 'TACGATCGATTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCGATCGTAGCATCGATCGATCGAAAAAAAACATCGATCGATCGATCGATCGACTTTTT'
    #seq = 'TTGTTAACATCAAAGCTAATGAATGACATTGCAGATATTTGTAAAAGTTTAGTAAGTATGCTTCCTGTTTTGCTATCATATTTTGATTCTAATAGGCATA'
    seq = 'CCACCGGCCTGGAGCGCGCCGCCAAGGCGGCGCGCGAGCTGGAGCACTCGCGTGAGTgcggcggggcggggcggggcgggcgggcgggcgggacgggCCG'
    #entropy = estimate_shannon_entropy(seq)
    #print(entropy)
    #print(calculate_repeat_masked_percent(seq))

    binary='/home/a_vikhorev/soft/miniconda3/envs/mrbait/bin/blastn'
    db='/home/a_vikhorev/Data/DataBases/hg38_blast/hg38'

    blastresult = blastn(binary, seq, db, 10, 'out.outfmt6', 4)

    print(read_blast_results(blastresult))