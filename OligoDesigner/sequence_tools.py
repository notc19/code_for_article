import math
import collections

from Bio.SeqUtils import MeltingTemp as mt

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
    return round(mt.Tm_GC(seq), 3)

if __name__ == '__main__':
    seq = 'TACGATCGATTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCGATCGTAGCATCGATCGATCGAAAAAAAACATCGATCGATCGATCGATCGACT'
    entropy = estimate_shannon_entropy(seq)
    #print(entropy)

