import os
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
import xgboost as xgb
import joblib

import sequence_tools as st

rf_clf = None
model_ontarget_offtarget = None
model_blast = None

def init_models():
    """
    Load models with single-threaded settings. 
    Called once in each spawning process (in run_oligo_design.py)
    """
    global rf_clf, model_ontarget_offtarget, model_blast

    rf_clf = joblib.load("../lib/rf_dna_classifier.json")
    rf_clf.n_jobs = 1

    model_ontarget_offtarget = xgb.XGBRegressor()
    model_ontarget_offtarget.load_model("../lib/ontarget_offtarget_reg_model.json")
    model_ontarget_offtarget._Booster.set_param('n_jobs', 1)

    model_blast = xgb.XGBRegressor()
    model_blast.load_model("../lib/offt_xgboost_model.json")
    model_blast._Booster.set_param('n_jobs', 1)


def predict_coverage_class(locus_characteristics):
    """
    Функция передает массив со значениями 14 характеристик таргетного участка ДНК
    модели Random Forest, которая предсказывает класс покрытия таргета
    :param locus_characteristics: numpy-массив характеристик тарегтного участка
    :return: предсказанный класс покрытия ('high', 'moderate' or 'low' coverage)
    """
    locus_characteristics_array = np.array(locus_characteristics).reshape(1, -1)
    coverage_class = rf_clf.predict(locus_characteristics_array)
    return coverage_class

def predict_ontarget(olig_characteristics) -> float:
    """
    Эта функция принимает numpy-массив из характеристик олигонуклеотида и передает модели XGBoost
    :param olig_characteristics: numpy-массив характеристик олигонуклеотида
    :return: Количество on-target Nbases, предсказанное регрессионной моделью
    """
    quality = (model_ontarget_offtarget.predict(olig_characteristics))
    return quality

def blastn(seq:str, blast_binary:str, blast_db:str,  threads:int):
    """
    Эта функция выполняет BLAST с заданными параметрами для поданной на вход последовательности и возвращает
    список результатов BLAST в формате outfmt 6
    :param seq: последовательность олигонуклеотида (строка)
    :param blast_binary: путь до binary файла BLAST
    :param blast_db: путь до базы данных BLAST
    :param threads: сколько ядер использовать
    :return: список результатов BLAST в формате outfmt 6
    """
    binary = blast_binary
    blastn_cli = [
        binary,
        "-db", str(blast_db),
        "-evalue", str(1e-02),
        "-outfmt", '6 qseqid sseqid pident length sstart send qlen',
        "-num_threads", str(threads),
        "-perc_identity", str(10),
        "-qcov_hsp_perc", str(10)
       ]

    proc = Popen(blastn_cli, stdout=PIPE, stdin=PIPE, env={'PATH': os.getenv('PATH')})
    t = proc.communicate(str.encode(str(seq)))
    blastresult = t[0].decode().split('\n')
    return blastresult

def predict_offtarget(blastresult: list) -> tuple[int, int]:
    """
    Эта функция читает результат BLAST, преобразует его в dataframe,
    рассчитывает % ЭИ BLAST-хитов и количество BLAST-хитов в интервалах % ЭИ,
    а затем передает массив модели XGBoost
    :param blastresult: Список, содержащий результаты BLAST в формате outfmt 6
    :return: Количество off-target BLAST-hits и Nbases, предсказанных регрессионной моделью 
    """
    columns=['Olig_id', 'Chromosome', 'pident', 'length', 'Start', 'End', 'qlen']
    df = pd.DataFrame([line.split('\t') for line in blastresult[:-1]], columns=columns).astype({
        'Olig_id': str,
        'Chromosome': str,
        'pident': float,
        'length': int,
        'Start': int,
        'End': int,
        'qlen': int})
    
    if (len(df)==1) and df.pident.values[0]==100:
        # no off-target - no XGBoost prediction
        blast_offt_hits = 0
        offtarget_nbases = 0

    else:
        df = df.assign(eff_ident = (df.pident * df.length) / df.qlen)
        df = df.loc[:,['Olig_id','eff_ident']]

        if len(df[df['eff_ident']==100]) > 1:
        # 100% off-target - нежелательный олиг, offtarget_nbases "максимальный"
            blast_offt_hits = len(df) - 1
            offtarget_nbases = 1000000000
        else:
            blast_offt_hits = len(df) - 1
            df = df.assign(percentage_tail = (df.eff_ident//10)*10)
            df = pd.pivot_table(df, index='Olig_id', values='eff_ident', columns='percentage_tail', aggfunc='count')
            if df.shape[1] == 9:
                X = np.array(df.iloc[:,0:8])
      
            elif df.shape[1] < 9:
                X = [df[i].iloc[0] if i in df.columns else 0 for i in range(20, 100, 10)]
                X = np.array([X])
                
            offtarget_nbases = abs(int(model_blast.predict(X)[0]))
    return blast_offt_hits, offtarget_nbases

def get_blast_offtarget(seq: str, blast_binary: str, blast_db: str, threads: int):
    """
    Краткая реализация предсказания off-target Nbases BLAST для последовательности олигонуклеотида
    :param seq: последовательность олигонуклеотида
    :param blast_binary: путь до binary файла BLAST
    :param blast_db: путь до базы данных BLAST
    :param threads: сколько ядер использовать
    :return: список результатов BLAST в формате outfmt 6
    """
    blast_result = blastn(seq, blast_binary, blast_db, threads)
    blast_offt_hits, offtarget_nbases = predict_offtarget(blast_result)

    return blast_offt_hits, offtarget_nbases

if __name__ == '__main__':
    import pandas as pd
    import DNAFMIndex
    from Bio import SeqIO
    import sequence_tools as st

    init_models()
    genome_index = '/home/a_makarova/Data/GRCh38_index.awfmi'
    fmindex = DNAFMIndex.read_index_from_file(genome_index, False)
    binary = '/home/a_makarova/miniconda3/envs/od/bin/blastn'
    db = '/home/a_makarova/Data/blastdb_GCA/hg38'    

    seq = 'CTCCCCAGTAGCTAGGATTGCAGGCACATGCTACCACACCTGGTTAATTTTTGTATTTTTAGTAGAGACGGGGTTTCACGATGTTGGCCAGGCTGGCTGG'
    seq_characteristics = list(st.get_olig_stats(seq, fmindex))
    seq_characteristics.append(10) #example intersection percent
    seq_characteristics = np.array(seq_characteristics).reshape(1, -1)
    regr_result = predict_ontarget(seq_characteristics)
    print(regr_result)

    locus_chars = st.get_locus_stats(seq)
    class_result = predict_coverage_class(locus_chars)

    fasta_path = '/home/a_makarova/Data/oligs/United_second_panel_oligs.fasta' #'../tests/test_results/OUTPUT_Oligonucleotides_RAW.fasta'
    d = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(record.seq)
        blast_offt_hits, offt_nbases = get_blast_offtarget(seq_str, binary, db, 20)
        d[record.id] = [blast_offt_hits, offt_nbases]
    df = pd.DataFrame.from_dict(d, orient='index', columns=['offt_BLAST_hits', 'offt_Nbases'])
    df.to_csv('../tests/test_results/Second_panel_oligs_BLAST_offt', sep='\t')
