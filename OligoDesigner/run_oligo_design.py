import multiprocessing
import logging
import logging.handlers
import threading
from logging import Logger
import traceback
from datetime import datetime
import argparse
import numpy as np
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import DNAFMIndex

import sequence_tools as st
import models

params = None
logger = None
fmindex = None

def _log_listener_thread(queue: multiprocessing.Queue, handlers: list):
    """
    Runs in a background thread in the main process.
    Pulls LogRecords from the queue and passes them to the real handlers.
    Stops when it receives a None sentinel.
    """
    while True:
        record = queue.get()
        if record is None:
            break
        for handler in handlers:
            try:
                handler.handle(record)
            except Exception:
                pass

def start_log_listener(queue: multiprocessing.Queue, handlers: list) -> threading.Thread:
    """Start the listener thread and return it so the caller can join() it later."""
    t = threading.Thread(target=_log_listener_thread, args=(queue, handlers), daemon=True)
    t.start()
    return t

def stop_log_listener(queue: multiprocessing.Queue, thread: threading.Thread):
    """Send the sentinel and wait for the listener thread to finish."""
    queue.put(None)
    thread.join()

def create_logger(timestamp: datetime) -> Logger:
    """
    Create the main-process logger.
    Returns the Logger object; its handlers are also passed to the listener thread.
    """
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(fmt)

    file_handler = logging.FileHandler(
        f'log_oligodesigner_{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}-{timestamp.second}.txt',
        mode='w'
    )
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(fmt)

    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.handlers = []
    root.addHandler(stream_handler)
    root.addHandler(file_handler)

    return root, [stream_handler, file_handler]

def worker_init(parameters_main, log_queue: multiprocessing.Queue):
    """
    Pool initializer - called once per spawned worker process.
    Prevents logging and threading conflicts between multiple processes. 
    Initializes models with n_jobs=1.
    :param parameters_main: parameters taken from def main()
    :log_queue: queue of logger objects
    """
    global params, logger, fmindex
    params = parameters_main

    root = logging.getLogger()
    for handler in root.handlers[:]:
        try:
            handler.close()
        except Exception:
            pass
        root.removeHandler(handler)

    for name in list(logging.Logger.manager.loggerDict.keys()):
        log = logging.getLogger(name)
        for handler in log.handlers[:]:
            try:
                handler.close()
            except Exception:
                pass
            log.removeHandler(handler)

    queue_handler = logging.handlers.QueueHandler(log_queue)
    root.addHandler(queue_handler)
    root.setLevel(logging.INFO)

    logger = logging.getLogger(multiprocessing.current_process().name)
    fmindex = DNAFMIndex.read_index_from_file(params.genome_index, False)
    models.init_models()    

class Oligonucleotide():
    def __init__(self, idi, seq, length, chr, start, end, gc, entr, mt, masked, monorepeat, direpeat, trirepeat, offtarget):
        self.idi = idi
        self.seq = str(seq)
        self.length = length
        self.chr = chr
        self.start = start
        self.end = end

        self.gc = float(gc)
        self.entr = float(entr)
        self.mt = float(mt)
        self.masked = float(masked)
        self.monorepeat = float(monorepeat)
        self.direpeat = float(direpeat)
        self.trirepeat = float(trirepeat)
        self.offtarget = float(offtarget)
        self.intersection = None
        self.qual = None

    def __str__(self):
        return f'{self.idi} {self.gc} {self.entr} {self.mt} {self.masked} {self.monorepeat} {self.direpeat} {self.trirepeat} {self.offtarget} {self.intersection} {self.length} {self.qual} {self.seq} {self.chr} {self.start} {self.end}'
    
    def __repr__(self):
        return f'{self.idi} {self.gc} {self.entr} {self.mt} {self.masked} {self.monorepeat} {self.direpeat} {self.trirepeat} {self.offtarget} {self.intersection} {self.length} {self.qual} {self.seq} {self.chr} {self.start} {self.end}'

def olig_intersection_strategy(coverage_result: str) -> int:
    """
    Эта функция возвращает кратность перекрытия зондов в зависимости от класса сложности 
    участка дизайна, рассчитывается перекрытие зондов,
    необходимое для эффективного покрытия участка ДНК.
    :param coverage_result: класс участка - результат
                            предсказания моделью классификации
    :return: кратность перекрытия зондов
    """    
    if coverage_result == 'high_coverage':
        return 1
    elif coverage_result == 'moderate_coverage':
        return 2
    return 4

def choose_k_intersecting_oligos(scores: list[float], k: int, lambda_dist: int=1000) -> list[int]:
    """
    Function chooses k best-scored oligos in window.
    Strategy is to choose best oligos densely placed within window.
    :param scores: list of olig.qual scores
    :param k: кратность перекрытия зондов
    :param lambda_dist: coefficient for olig distance sum
    :return: список индексов лучших олигонуклеотидов
    """     
    l = len(scores)
    k = min(k, l)

    if k==1:
        return [int(np.argmax(scores))]

    best = None
    best_val = -np.inf

    for combination in combinations(range(l), k):
        score_sum = sum(scores[i] for i in combination)

        dist_sum = 0
        for i in range(k):
            for j in range(i + 1, k):
                dist_sum += abs(combination[i] - combination[j])

        val = score_sum + lambda_dist * dist_sum

        if val > best_val:
            best_val = val
            best = combination

    return list(best)
    
def select_best_in_window(sequence: SeqRecord, chr: str, abs_start: int) -> list[Oligonucleotide]:
    """
    Эта функция выбирает лучшие олигонуклеотиды в окне дизайна.
    В зависимости от сложности участка дизайна, выбирается такое количество лучших зондов,
    которое необходимо для эффективного покрытия участка ДНК.

    :param sequence: последовательость (объект SeqRecord)
    :param chr: chromosome name from SeqRecord description
    :param abs_start: sequence start from SeqRecord description
    :return: список олигонуклеотидов (объектов класса Oligonucleotide)
    """
    oligos = []
    scores = []

    window_characteristics = st.get_locus_stats(sequence.seq)
    window_coverage = models.predict_coverage_class(window_characteristics)
    k = olig_intersection_strategy(window_coverage[0])

    for i in range(len(sequence) - params.length + 1):
        subseq = str(sequence.seq[i:i+params.length])

        start = abs_start+i-1
        olig = Oligonucleotide(sequence.id, subseq, len(subseq), chr, start, start+params.length, *st.get_olig_stats(subseq, fmindex))
        olig.intersection = (k-1)*100 + params.max_i
        olig_characteristics = np.array([olig.gc, olig.entr, olig.mt, olig.masked, olig.monorepeat, 
                                         olig.direpeat, olig.trirepeat, olig.offtarget, olig.intersection]).reshape(1,-1)

        predicted_nbases_array = models.predict_ontarget(olig_characteristics)
        predicted_nbases = float(predicted_nbases_array[0])

        olig.qual = round(predicted_nbases)
        oligos.append(olig)
        scores.append(predicted_nbases)

    best_indices = choose_k_intersecting_oligos(scores, k)
    best_oligos = [oligos[i] for i in best_indices]
    return best_oligos

def design_oligos(sequence: SeqRecord) -> list[Oligonucleotide]:
    """
    Эта функция выполняет дизайн олигонуклеотидов для заданного региона.
    :param sequence: последовательость (объект SeqRecord)
    :return: список олигонуклеотидов (объектов класса Oligonucleotide)
    """
    chrom, abs_start, abs_end = st.get_absolute_seq_coordinates(sequence)
    seq = str(sequence.seq)
    olig_list = []

    subseq = seq[:params.length + abs(params.min_i)]

    # create first olig - just first {length} nucleotides from start
    try:
        oligos = select_best_in_window(SeqRecord(Seq(subseq), id = sequence.id), chrom, abs_start)
        olig_list.extend(oligos)

    except Exception as e:
        raise e
    
    olig = olig_list[-1]
    while True:
        previous_end = seq.index(olig.seq) + len(olig.seq) - 1
        start = previous_end-params.max_i
        end = previous_end-params.min_i if previous_end+params.min_i < len(seq) else len(seq)
        seq = str(sequence.seq)
        new_subseq = seq[start:end + params.length]

        if len(new_subseq) >= params.length:
            try:
                oligos = select_best_in_window(SeqRecord(Seq(new_subseq), id = sequence.id), chrom, abs_start+start)
                olig_list.extend(oligos)
                olig = olig_list[-1]
            except Exception as e:
                raise e

        elif len(new_subseq) >= params.length * 0.7:
            olig = Oligonucleotide(sequence.id, seq[-params.length:], len(seq[-params.length:]),
                                   chrom, abs_end-params.length, abs_end, *st.get_olig_stats(seq[-params.length:], fmindex))
            olig.intersection = 30
            olig_characteristics = np.array([olig.gc, olig.entr, olig.mt, olig.masked, olig.monorepeat, 
                                    olig.direpeat, olig.trirepeat, olig.offtarget, olig.intersection]).reshape(1,-1)
            predicted_nbases_array = models.predict_ontarget(olig_characteristics)
            olig.qual = round(float(predicted_nbases_array[0]))
            olig_list.append(olig)
            break
        else:
            break

        if end >= len(seq):
            break
    return olig_list

def design_oligos_snp(sequence: SeqRecord) -> list[Oligonucleotide]:
    """
    Эта функция выполняет дизайн олигонуклеотида для заданного региона с SNP.
    :param sequence: последовательость (объект SeqRecord)
    :return: список олигонуклеотидов (объектов класса Oligonucleotide)
    """
    chrom, abs_start, _ = st.get_absolute_seq_coordinates(sequence)
    oligs = select_best_in_window(sequence, chrom, abs_start)

    return oligs

def safe_design_oligos(sequence: SeqRecord):
    try:
        return design_oligos(sequence)
    except Exception as e:
        raise RuntimeError(f"Worker failed on {sequence.id}:\n{traceback.format_exc()}") from e

def safe_design_oligos_snp(sequence: SeqRecord):
    try:
        return design_oligos_snp(sequence)
    except Exception as e:
        raise RuntimeError(f"Worker failed on {sequence.id}:\n{traceback.format_exc()}") from e

def run_oligo_design(seq_list: list[SeqRecord], threads:int, parameters_main, logger: Logger, log_handlers: list) -> list[Oligonucleotide]:
    """
    Эта функция запускает дизайн олигонуклеотидов для списка последовательностей SeqRecord.
    :param seq_list: список последовательностей SeqRecord
    :param threads: сколько ядер использовать
    :param logger: объект класса Logger
    :param parameters_main: parameters from def main()
    :param log_handlers: хэндлеры логгера
    :return: список олигонуклеотидов (объектов класса Oligonucleotide)
    """
    logger.info('Start designing oligonucleotides!')

    global params
    params = parameters_main

    snp_seq_list, other_seq_list = [], []
    for seq in seq_list:
        if seq.name == "snp" or "snp" in seq.description:
            snp_seq_list.append(seq)
        else:
            other_seq_list.append(seq)

    logger.info(
        f'Found {len(snp_seq_list)} SNP sequences and {len(other_seq_list)} other regions'
    )

    multiproc_context = multiprocessing.get_context('spawn')

    log_queue = multiproc_context.Queue()
    listener_thread = start_log_listener(log_queue, log_handlers)

    output_oligos = []

    try:
        if snp_seq_list:
            logger.info('Designing oligos for SNPs...')
            with multiproc_context.Pool(
                threads,
                initializer=worker_init,
                initargs=(parameters_main, log_queue)
            ) as pool:
                result = pool.map_async(safe_design_oligos_snp, snp_seq_list)
                snp_oligos_tmp = result.get(timeout=3600)

            for oligos in snp_oligos_tmp:
                output_oligos.extend(oligos)
            logger.info('SNP oligos done.')

        if other_seq_list:
            logger.info('Designing oligos for other regions...')
            n_workers = min(len(other_seq_list), threads)
            with multiproc_context.Pool(
                n_workers,
                initializer=worker_init,
                initargs=(parameters_main, log_queue)
            ) as pool:
                result = pool.map_async(safe_design_oligos, other_seq_list)
                other_oligos_tmp = result.get(timeout=3600)

            for oligos in other_oligos_tmp:
                output_oligos.extend(oligos)
            logger.info('Other region oligos done.')

    finally:
        stop_log_listener(log_queue, listener_thread)

    return output_oligos

if __name__ == '__main__':
    timestamp = datetime.now()
    logger, log_handlers = create_logger(timestamp)
    logger.info(f'######### Start OligoDesigner at {timestamp.day}.{timestamp.month}.{timestamp.year} {timestamp.hour}:{timestamp.minute}:{timestamp.second} #########')
    inpt = '../tests/test_results/INPUT_regions_sequences.fasta'
    seq_list = list(SeqIO.parse(inpt, "fasta"))

    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('--max_i', help='Maximum intersection of oligos in region',
                             type=int, required=False, default=10)
    args_parser.add_argument('--min_i', help='Minimum intersection of oligos in region',
                             type=int, required=False, default=-10)
    args_parser.add_argument('-l','--length', help='Length of the designed olignucleotide bait',
                             type=int, required=False, default=100)
    args_parser.add_argument('--threads', help='Number of threads',
                             type=int, required=False, default=1)
    args_parser.add_argument('--blast_strategy', help='Для каких целей использовать blast. Варианты: filter/analyze',
                             type=str, required=False, default='analyze')
    args_parser.add_argument('--blast_binary', help='Путь до исполняемого файла blastn',
                             type=str, required=False, default='')
    args_parser.add_argument('--blast_db', help='Путь до базы данныx blast',
                             type=str, required=False, default='')
    args_parser.add_argument('--blast_filtering_params', help='Параметры для фильтрации результатов бласта (identity, coverage, степень e-value)',
                             type=str, required=False, default='0.10, 0.10, 2')
                             
    args_parser.add_argument('--genome_index', help='path to file with genome index in awfmi format (build by DNAFMIndex)',
                             type=str, required=False, default='../Data/GRCh38_index.awfmi')
    params = args_parser.parse_args()

    oligos = run_oligo_design(seq_list, params.threads, params, logger, log_handlers)
    for olig in oligos:
        print(olig,'\n')
