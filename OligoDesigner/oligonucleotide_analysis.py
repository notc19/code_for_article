from logging import Logger
import pandas as pd
from oligo_plot import Designer

def norm_qual(qual_data: dict) -> dict:
    """
    Нормализует значения качества в qual_data на глобальный максимум,
    приводя все значения к диапазону [0, 1].
    Метод изменяет qual_data на месте. Если qual_data пуст или не является
    словарём, либо глобальный максимум равен нулю, метод завершается без изменений.
    :param qual_data: словарь с метриками качества, рассчитанными во всех 'окнах' таргетного региона
    :return: словарь с нормализованными метриками качества данного региона
    """
    if not qual_data or not isinstance(qual_data, dict):
        return
    normal_lim = 0
    for i in qual_data.values():
        for j in i:
            normal_lim = max(normal_lim, j)
    if normal_lim == 0:
        return
    for i in qual_data.keys():
        for j in range(len(qual_data[i])):
            qual_data[i][j] = qual_data[i][j] / normal_lim
    return qual_data

def analyze_oligos_draw_plots(oligo_df:pd.DataFrame, qual_dict:dict, tmpdir, logger:Logger, outdir) -> None:
    """
    Эта функция выполняет анализ разработанной панели олигонуклеотидов
    :param oligo_df: Датафрейм с олигонуклеотидами
    :param tmpdir: Путь к директории с временными файлами
    :param logger: Объект класса Logger
    :param outdir: Путь к выходной директории
    :return: Ничего не возвращает
    """
    logger.info(f'The mean GC content of the designed panel is {round(oligo_df["GC"].mean(), 3)}')
    logger.info(f'The mean Shannon Entropy of the designed panel is {round(oligo_df["Entropy"].mean(), 3)}')
    logger.info(f'The mean Melting Temperature of the designed panel is {round(oligo_df["Tm"].mean(), 3)}')
    logger.info(f'The mean Repeat Length of the designed panel is {round(oligo_df["%1-nt_repeats"].mean(), 3)}')
    logger.info(f'The mean Masked Length of the designed panel is {round(oligo_df["Masked"].mean(), 3)}')
    logger.info(f'The mean Off-Target Hits of the designed panel is {round(oligo_df["Off-target_100%"].mean(), 3)}')

    normalized_qual_dict = norm_qual(qual_dict)

    des = Designer(f'{tmpdir}/INPUT_coords_FIXED.bed', f'{tmpdir}/OUTPUT_Oligonucleotides.bed', 
                   quality_dict = normalized_qual_dict, result_path = outdir)    
    des.save()
