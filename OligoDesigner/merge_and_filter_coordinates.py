import logging
from pathlib import Path

def drop_nested(coords_sorted:list, logger:logging) -> list:
    """
    Эта функция удаляет из списка вложенные списки координат. Координаты считаются вложенными в двух случаях.
    1. Если начальная координата региона совпадает с начальной координатой следующего региона, а координата конца меньше или равна координате конца следующего региона
    2. Если начальная координата региона больше, чем начальная координата следующего региона, но координата конца меньше, чем координата конца следующего региона
    :param coords_sorted: Список регионов, в котором под индексом 1 координата начала, а под индексом 2 координата конца
    :param logger: Объект класса logging
    :return: Список, из которого удалены вложенные координаты
    """
    clear_list = []
    for idx in range(len(coords_sorted)):
        if idx < len(coords_sorted)-1 and \
            (int(coords_sorted[idx][1]) == int(coords_sorted[idx+1][1]) and \
            int(coords_sorted[idx][2]) <= int(coords_sorted[idx+1][2])):
            var = ' '.join(map(str,coords_sorted[idx]))
            logger.debug(f'Coordinates {coords_sorted[idx][1]}:{coords_sorted[idx][2]} and {coords_sorted[idx+1][1]}:{coords_sorted[idx+1][2]} are nested. Skipping {var} variant!')
            if (int(coords_sorted[idx][2])-int(coords_sorted[idx][1])) > \
                    (int(coords_sorted[idx+1][2])-int(coords_sorted[idx+1][1])):
                logger.error('BUG! Coordinates are not really nested!')
            continue
        if idx > 0 and \
            int(coords_sorted[idx][1]) >  int(coords_sorted[idx-1][1]) and \
            int(coords_sorted[idx][2]) <= int(coords_sorted[idx-1][2]):
            var = ' '.join(map(str,coords_sorted[idx]))
            logger.debug(f'Coordinates {coords_sorted[idx][1]}:{coords_sorted[idx][2]} and {coords_sorted[idx-1][1]}:{coords_sorted[idx-1][2]} are nested. Skipping {var} variant!')
            if (int(coords_sorted[idx][2])-int(coords_sorted[idx][1])) > \
                    (int(coords_sorted[idx-1][2])-int(coords_sorted[idx-1][1])):
                logger.error('BUG! Coordinates are not really nested!')
            continue
        clear_list.append(coords_sorted[idx])
        
    return clear_list

def merge_regions(coords_sorted:list, logger:logging, counter:int) -> [list, int]:
    """
    Эта функция сливает координаты пересекающихся регионов. Регионы считаются пересекающимися, если координата начала региона лежит
    между координатами начала и конца предыдущего региона, а координата конца лежит за пределами предыдущего региона.
    :param coords_sorted: Список регионов, в котором под индексом 1 координата начала, а под индексом 2 координата конца
    :param logger: Объект класса logging
    :param counter: Число, счетчик пересекающихся регионов
    :return: Список, в котором удалены пересекающиеся координаты, и счетчик пересекающихся координат
    """
    clear_list = []
    for idx in range(len(coords_sorted)):
        if idx > 0 and \
        int(coords_sorted[idx][1]) <= int(coords_sorted[idx-1][2]) and \
        int(coords_sorted[idx][1]) >  int(coords_sorted[idx-1][1]) and \
        int(coords_sorted[idx][2]) >  int(coords_sorted[idx-1][2]):
            logger.debug(f'Coordinates {coords_sorted[idx][1]}:{coords_sorted[idx][2]} and {coords_sorted[idx-1][1]}:{coords_sorted[idx-1][2]} are intersecting')
            coords_sorted[idx][1] = coords_sorted[idx-1][1]
            if len(coords_sorted[idx]) == 7:
                coords_sorted[idx][4] = coords_sorted[idx][4] + ',' + coords_sorted[idx-1][4]
                coords_sorted[idx][5] = coords_sorted[idx][5] + ',' + coords_sorted[idx-1][5]
            clear_list.append(coords_sorted[idx])
            counter +=1
        elif idx < len(coords_sorted)-1 and \
        int(coords_sorted[idx][2]) >= int(coords_sorted[idx+1][1]) and \
        int(coords_sorted[idx][1]) <  int(coords_sorted[idx+1][1]) and \
        int(coords_sorted[idx][2]) <  int(coords_sorted[idx+1][2]):
            continue
        else:
            clear_list.append(coords_sorted[idx])
    return clear_list, counter

def process_coordinates(input_bed_list:list, write_intermediate:bool, tmpdir:Path, logger:logging) -> list:
    """
    Эта функция принимает на вход список координат, и удаляет из него вложенные координаты и объединяет пересекающиеся.
    :param input_bed_list: Список регионов, в котором под индексом 1 координата начала, а под индексом 2 координата конца
    :param write_intermediate: Сохранять или нет промежуточный файл с регионами в формате bed
    :param tmpdir: Объект Path, содержащий путь до директории с промежуточными файлами
    :param logger: Объект класса logging
    :return: Список с исправленными координатами
    """
    fixed_bed_list = []
    logger.info(f'Start fixing coordinates! There were found {len(input_bed_list)} coordinates in input')
    coord_chrom_dict = {}
    for coord in input_bed_list:
        coord_list = coord.split()
        coord_chrom_dict.setdefault(coord_list[0], []).append(coord_list)

    nested_counter, merged_counter, iter_counter = 0, 0, 0
    for chrom in coord_chrom_dict:
        coords_sorted = sorted(coord_chrom_dict[chrom], key=lambda x: (int(x[1]),int(x[2])))
        while True:
            coords_sorted_wo_nested = drop_nested(coords_sorted, logger)
            nested_counter += len(coords_sorted) - len(coords_sorted_wo_nested)
            coords_sorted_merged, merged_counter = merge_regions(coords_sorted_wo_nested, logger, merged_counter) \
                if len(coords_sorted_wo_nested) > 1 else (coords_sorted_wo_nested, merged_counter)
            if len(coords_sorted) != len(coords_sorted_merged):
                coords_sorted = sorted(coords_sorted_merged, key=lambda x: (int(x[1]),int(x[2])))
                iter_counter += 1
                logger.info(f'Iteration {iter_counter} was done!')
            else:
                break
        fixed_bed_list.extend(coords_sorted_merged)
        
    logger.info(f'Coordinates were fixed! There were {nested_counter} nested coordinates and {merged_counter} intersecting coordinates. {len(fixed_bed_list)} coordinates were saved')

    fixed_bed_list = list(map(('\t').join,fixed_bed_list))

    if write_intermediate:
        with open(tmpdir/'FIXED_coords.bed', 'w') as outfile:
            logger.info(f"Coordinates for input genes were written in {tmpdir/'FIXED_coords.bed'} file")
            #print('\n'.join(fixed_bed_list))
            outfile.write('\n'.join(fixed_bed_list))

    return fixed_bed_list


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    input_bed_list = ['chr1\t100\t200\tENSG00000160072\tENSE00003889014\fghfgh\t1',
                      'chr1\t100\t200\tENSG00000160072\tENSE00003889014\texon\t1',
                      'chr1\t100\t500\tENSG00000160072\tENSE00003889014\texon\t1',
                      'chr1\t100\t250\tENSG00000160072\tENSE00003889014\texon\t1',
                      'chr1\t120\t300\tENSG00000160072\tENSE00003889014\texon\t1',
                      'chr1\t200\t800\tENSG00000160072\tfull\tfull\t1',
                      'chr1\t400\t450\tENSG00000160072\tfull\tfull\t1',
                      'chr1\t700\t900\tENSG00000160072\tfull\tfull\t1',
                      'chr1\t7100\t9400\tENSG00000160072\tfull\tfull\t1',
                      'chr1\t900\t1000\tENSG00000160072\tfull\tfull\t1']
    
    print(process_coordinates(input_bed_list, False, '', logger))