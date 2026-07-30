import DNAFMIndex
import numpy as np
from scipy.signal import find_peaks, peak_widths
from Bio import SeqIO

def create_fm_index_from_genome(path_to_genome: str, path_to_index: str | None):
    if path_to_index == None:
        path_to_index = path_to_genome + '.awfmi'
    config = DNAFMIndex.IndexConfiguration(8, 12, 2, True, False)
    fmindex = DNAFMIndex.Index(config, path_to_index, fasta_path=path_to_genome)
    return fmindex

def get_locus_pos(locus: str, index: DNAFMIndex.Index) -> int:
    locus_search = DNAFMIndex.KmerSearchList(2)
    locus_search.fill([locus])
    locus_search.parallel_search_locate(index, 1)

    locus_start = locus_search.kmer_search_data[0].position_list[0]

    locus_end = locus_start + len(locus)
    return locus_start, locus_end

def merge_intervals(left, right):
    intervals = list(zip(left, right))
    intervals.sort()
    
    merged = []
    current_start, current_end = intervals[0]
    
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.extend([current_start, current_end])
            current_start, current_end = start, end
    
    merged.extend([current_start, current_end])
    
    return np.array(merged)

def get_prefixes(locus: str, kmer_len: int, olig_len: int) -> list[str]:
    L = len(locus)
    prefixes = []

    for j in range(L - kmer_len + 1):
        if j+olig_len <= L:
            subseq = locus[j:j+olig_len]
            prefixes += [subseq[:i] for i in range(kmer_len, olig_len+1)]
        else:
            subseq = locus[j:L]
            prefixes += [subseq[:i] for i in range(kmer_len, L-j+1)]

    return prefixes

def kmer_counter(locus: str, index: DNAFMIndex.Index, kmer_len: int, olig_len: int, threads: int):

    prefixes = get_prefixes(locus, kmer_len, olig_len)

    kmer_search_list = DNAFMIndex.KmerSearchList(len(locus)*100) #capacity
    kmer_search_list.fill(prefixes)
    kmer_search_list.parallel_search_locate(index, threads)
    
    locus_data = {}
    for i in range(kmer_search_list.count):
        prefix =  kmer_search_list.kmer_search_data[i].kmer_string.decode()
        count = kmer_search_list.kmer_search_data[i].count
        locus_data[prefix] = count

    return locus_data

def build_offtarget_alignment(locus: str, locus_data: dict, kmer_len: int, olig_len: int):

    L = len(locus)
    max_prefix_len = olig_len

    alignment = np.zeros((L, max_prefix_len - kmer_len + 1), dtype=np.int32)
    for i in range(L - kmer_len + 1):
        max_len = min(olig_len, L - i)
        
        for prefix_len in range(kmer_len, max_len + 1):
            prefix = locus[i:i+prefix_len]
            count = locus_data[prefix]
            col_idx = prefix_len - kmer_len
            alignment[i, col_idx] = 0 if count in (0, 1) else count - 1

    return alignment

def analyze_peaks(alignment: np.ndarray, min_height, min_width):
    result = {
        "peaks": np.array([]),
        "heights": np.array([]),
        "widths": np.array([]),
        "borders": [],
        "local_minimums": [],
        "local_min_heights": [],
    }

    peaks, properties = find_peaks(alignment, height=min_height, width=min_width, rel_height=0.9)
    if len(peaks) == 0:
        return result
    
    result["peaks"] = peaks
    result["heights"] = properties.get("peak_heights", np.array([]))
    result["widths"] = properties.get("widths", np.array([]))
    
    try:
        widths, _, left_bound, right_bound = peak_widths(alignment, peaks, rel_height=0.9)
        borders = merge_intervals(left_bound, right_bound)
        result["borders"] = borders
    except Exception:
        result["borders"] = []

    if len(peaks) >= 2:
        local_minimums = []
        local_min_heights = []

        for i in range(len(peaks) - 1):
            left, right = peaks[i], peaks[i + 1]
            segment = alignment[left:right]

            if len(segment) < 3:
                continue

            min_peaks, _ = find_peaks(-segment)
            if len(min_peaks) == 0:
                continue

            min_idx = np.argmin(alignment[left + min_peaks])
            absolute_idx = min_peaks[min_idx] + left

            local_minimums.append(absolute_idx)
            local_min_heights.append(alignment[absolute_idx])

        result["local_minimums"] = local_minimums
        result["local_min_heights"] = local_min_heights

    return result
    
    
if __name__ == '__main__':
    import pandas as pd
    import numpy as np
    import DNAFMIndex
    #path_to_genome = "../Data/GRCh38.primary_assembly.genome.fa"
    #fmindex = create_fm_index_from_genome(path_to_genome, None)
    genome_index = '/home/a_makarova/Data/GRCh38_index.awfmi'
    fmindex = DNAFMIndex.read_index_from_file(genome_index, False)
    
    inpt = "/home/a_makarova/Data/oligs/United_second_panel_oligs.fasta"
    seq_list = list(SeqIO.parse(inpt, "fasta"))
    sum_dict = {}
    for i in range(len(seq_list)):
        locus = str(seq_list[i].seq)
        id = str(seq_list[i].id)
        locus_data = kmer_counter(locus, fmindex, 30, 100, 14)
        
        alignment = build_offtarget_alignment(locus, locus_data, 30, 100)
        aligned_sum = np.sum(alignment, axis=1)
        sum_dict[id] = aligned_sum

    df = pd.DataFrame.from_dict(sum_dict, orient='index')
    df.to_csv('../tests/test_results/Second_panel_DNAFM.tsv')

