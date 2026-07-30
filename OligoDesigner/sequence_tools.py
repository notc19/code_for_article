import math
import re
import collections
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils.MeltingTemp import DNA_NN4
import DNAFMIndex
import g4_hunter as g4

def get_absolute_seq_coordinates(sequence: SeqRecord) -> tuple[str, int, int]:
    pattern = r"chrom:(\S+)\s+start:(\d+)\s+end:(\d+)"
    match = re.search(pattern, sequence.description)

    chrom, start, end = match.groups()
    start, end = int(start), int(end)
    if start > end:
        minus_strand_start = end
        end = start
        start = minus_strand_start
    return chrom, start+1, end+1

def gc_content(seq:str) -> float:
    """
    Эта функция возвращает ГЦ-состав последовательности
    :param seq: последовательность олигонуклеотида
    :return: ГЦ-состав (доля нуклеотидов Г и Ц) олигонуклеотида
    """
    gc = (seq.count('G') + seq.count('C'))
    return round(gc / len(seq), 3)

def estimate_shannon_entropy(seq:str) -> float:
    """
    Эта функция вычисляет энтропию Шеннона для последовательности.
    Чем выше значение энтропии, тем более разнообразна структура последовательности. В идеале должна быть как можно ближе к 2.
    :param seq: последовательность ДНК
    :return: Значение энтропии (число от 0 до 2)
    """
    seq = seq.upper()
    bases = collections.Counter([tmp_base for tmp_base in seq])

    shannon_entropy_value = 0
    for base in bases:
        n_i = bases[base]
        p_i = n_i / float(len(seq))
        shannon_entropy_value -= p_i * math.log2(p_i)
 
    return round(shannon_entropy_value, 3)

def shannon_kmer_entropy(seq:str, k=2) -> float:
    """
    Эта функция вычисляет 2-mer энтропию Шеннона для последовательности.
    Чем выше значение энтропии, тем более разнообразна структура последовательности.
    :param seq: последовательность ДНК
    :return: Значение энтропии (число от 3 до 4)
    """
    seq = seq.upper()
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    counts = collections.Counter(kmers)
    counts_sum = sum(counts.values())
    
    kmer_entropy_value = 0
    for c in counts.values():
        p_i = c / counts_sum
        kmer_entropy_value -= p_i * math.log2(p_i)
    
    return round(kmer_entropy_value, 3)

def calculate_melting_temp(seq:str) -> float:
    """
    Эта функция вычисляет температуру плавления олигонуклеотида
    :param seq: последовательность олигонуклеотида
    :return: температура плавления олигонуклеотида
    """
    return round(mt.Tm_NN(seq, nn_table=DNA_NN4, Na=800, saltcorr=6), 3)

def calculate_repeats(seq: str) -> tuple[float]:
    """
    Эта функция считает отношение длины самого длинного моно-/ди-/три-нуклеотидного повтора к длине последовательности,
    а также суммарные доли моно-, ди-, три-нуклеотидных повторов относительно длины последовательности.
    :param seq: Последовательность олигонуклеотида (строка)
    :return: tuple с долями длин повторов к длине последовательности
    """
    seq = seq.upper()
    n = len(seq)

    max_lens = {1: 0, 2: 0, 3: 0}
    sum_lens = {1: 0, 2: 0, 3: 0}

    for k in (1, 2, 3):
        i = 0
        while i <= n - k:
            motif = seq[i:i+k]
            if any(c not in "ACGT" for c in motif):
                i += 1
                continue

            # проверка, что есть хотя бы 2 motif подряд
            if seq[i+k:i+2*k] != motif:
                i += 1
                continue

            # считаем длину повтора
            j = i + k
            repeat_len = k
            while j <= n - k and seq[j:j+k] == motif:
                j += k
                repeat_len += k

            # проверка минимального порога
            valid = False
            if k == 1:
                # мононукл. повторы начинаются от 5 нуклеотидов
                valid = repeat_len >= 5
            else:
                # ди-/три-нукл. повторы начинаются от трехкратного повторения
                valid = repeat_len >= 3 * k

            if valid:
                if repeat_len > max_lens[k]:
                    max_lens[k] = repeat_len
                sum_lens[k] += repeat_len
                # перескакиваем шаг, равный длине повтора
                i = j
            else:
                i += 1

    mono_max = round(max_lens[1] / n, 3)
    di_max   = round(max_lens[2] / n, 3)
    tri_max  = round(max_lens[3] / n, 3)

    mono_sum = round(sum_lens[1] / n, 3)
    di_sum   = round(sum_lens[2] / n, 3)
    tri_sum  = round(sum_lens[3] / n, 3)

    total_sum = round((sum_lens[1] + sum_lens[2] + sum_lens[3]) / n, 3)

    return mono_max, di_max, tri_max, mono_sum, di_sum, tri_sum, total_sum

def calculate_repeat_masked_percent(seq:str) -> float:
    """
    Эта функция считает процент последовательности, который был маскирован RepeatMasker-ом (маленькие буквы (sm) или N (rm)).
    :param seq: Последовательность олигонуклеотида
    :return: Процент маскированого повтора в последовательности (или 0.05, если повторов нет)
    """
    masked = sum(1 for base in seq if base.islower() or base == 'N')
    return round(masked / len(seq), 3)

def get_prefixes(seq: str, kmer_len: int, olig_len: int) -> list[str]:
    """
    Эта функция составляет список префиксов последовательности
    :param seq: Последовательность, минимальная длина k-mer, максимальная длина k-mer (равная длине олигонуклеотида) 
    :return: Список префиксов
    """
    L = len(seq)
    prefixes = []

    for j in range(L - kmer_len + 1):
        if j+olig_len <= L:
            subseq = seq[j:j+olig_len]
            prefixes += [subseq[:i] for i in range(kmer_len, olig_len+1)]
        else:
            subseq = seq[j:L]
            prefixes += [subseq[:i] for i in range(kmer_len, L-j+1)]

    return prefixes

def kmer_counter(prefixes: list[str], fmindex: DNAFMIndex.Index, threads: int) -> int:
    """
    Эта функция рассчитывает количество 100%-ных совпадений префиксов последовательности в геноме
    :param seq: Последовательность
    :param prefixes: Список префиксов
    :param index: геномный индекс
    :param threads: количество ядер
    :return: Количество 100% off-target для префиксов
    """
    kmer_search_list = DNAFMIndex.KmerSearchList(len(prefixes)*2) # KmerSearchList capacity
    kmer_search_list.fill(prefixes)
    kmer_search_list.parallel_search_locate(fmindex, threads)
    
    offtarget_sum = 0
    for i in range(kmer_search_list.count):
        if kmer_search_list.kmer_search_data[i].count > 1: # count only off-target
            prefix_count = kmer_search_list.kmer_search_data[i].count
            offtarget_sum += prefix_count

    return offtarget_sum

def g4_hunt(subseq: str):
    """
    Эта функция рассчитывает количество и score G4-квадруплексов для участка ДНК
    :param subseq: Последовательность участка ДНК
    :return: Кортеж из score, количества G4-квадруплексов
    """
    model = g4.G4Hunter(window_size=25, score_threshold=0.5)
    
    g4_score = float(model.run(subseq)['sum_score'])
    g4_num = len(model.run(subseq)['g4_clusters'])
    return g4_score, g4_num

def get_locus_stats(subseq: str):
    """
    Эта функция рассчитывает характеристики участка ДНК для модели классификации
    :param subseq: Последовательность участка ДНК
    :return: Кортеж из характеристик последовательности
    """
    length = len(subseq)
    gc = gc_content(subseq)
    entropy = estimate_shannon_entropy(subseq)
    kmer_entropy = shannon_kmer_entropy(subseq)
    mt = calculate_melting_temp(subseq)
    g4_score, g4_num = g4_hunt(subseq)
    mono_max, di_max, tri_max, mono_sum, di_sum, tri_sum, total_sum = calculate_repeats(subseq)
    
    return length, gc, entropy, kmer_entropy, g4_score, g4_num, mt, mono_max, di_max, tri_max, mono_sum, di_sum, tri_sum, total_sum

def get_olig_stats(subseq, fmindex):
    """
    Эта функция рассчитывает характеристики олигонуклеотидной последовательности.
    :param subseq: Последовательность олигонуклеотида
    :return: Кортеж из характеристик олигонуклеотида
    """
    gc = gc_content(subseq)
    entropy = estimate_shannon_entropy(subseq)
    mt = calculate_melting_temp(subseq)
    masked = calculate_repeat_masked_percent(subseq)
    
    repeats_info = calculate_repeats(subseq)
    mono, di, tri = repeats_info[3], repeats_info[4], repeats_info[5]
    
    prefixes = get_prefixes(subseq, 30, 100)
    offtarget = kmer_counter(prefixes, fmindex, 1)

    return gc, entropy, mt, masked, mono, di, tri, offtarget

if __name__ == '__main__':
    seq1 = 'TACGATCGATTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCGATCGTAGCATCGATCGATCGAAAAAAAACATCGATCGATCGATCGATCGACTTTTT'
    seq2 = 'TTGTTAACATCAAAGCTAATGAATGACATTGCAGATATTTGTAAAAGTTTAGTAAGTATGCTTCCTGTTTTGCTATCATATTTTGATTCTAATAGGCATA'
    seq3 = 'CCACCGGCCTGGAGCGCGCCGCCAAGGCGGCGCGCGAGCTGGAGCACTCGCGTGAGTgcggcggggcggggcggggcgggcgggcgggcgggacgggCCG'
    
    path_index = "../Data/GRCh38_index.awfmi"
    fmindex = DNAFMIndex.read_index_from_file(path_index, False)

    print(get_olig_stats(seq1, fmindex))
    print(get_olig_stats(seq2, fmindex))
    print(get_olig_stats(seq3, fmindex))
    
    print(get_locus_stats(seq1))
    print(get_locus_stats(seq2))
    print(get_locus_stats(seq3))