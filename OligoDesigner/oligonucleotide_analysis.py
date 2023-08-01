import logging
import pandas as pd

def analyze_oligos(oligo_df:pd.DataFrame, logger:logging) -> None:
    """
    Эта функция выполняет анализ разработанной панели олигонуклеотидов
    :param oligo_df: Датафрейм с олигонуклеотидами
    :param logger: Объект класса logging
    :return: Ничего не возвращает
    """
    oligo_df[["GC", "Entropy", "Tm", "Repeats", "Masked", "Blast_Hits"]] = \
        oligo_df[["GC", "Entropy", "Tm", "Repeats", "Masked", "Blast_Hits"]].apply(pd.to_numeric)

    logger.info(f'The mean GC content in designed oligos is {round(oligo_df["GC"].mean(), 3)}')
    logger.info(f'The mean Shannon Entropy in designed oligos is {round(oligo_df["Entropy"].mean(), 3)}')
    logger.info(f'The mean Melting Temperature in designed oligos is {round(oligo_df["Tm"].mean(), 3)}')
    logger.info(f'The mean Repeat Length in designed oligos is {round(oligo_df["Repeats"].mean(), 3)}')
    logger.info(f'The mean Masked Length in designed oligos is {round(oligo_df["Masked"].mean(), 3)}')
    logger.info(f'The mean Blast Hits in designed oligos is {round(oligo_df["Blast_Hits"].mean(), 3)}')