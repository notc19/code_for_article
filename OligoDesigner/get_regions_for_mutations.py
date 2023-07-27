import logging
from pathlib import Path

from constants import genome_versions_chromosome_dict, rev_dict

def get_coords_for_vcf(vcf:list, flank_length:int) -> str:
    """
    Эта функция принимает мутацию в формате vcf и выписывает координаты региона, фланкирующего мутацию, в формате bed
    :param vcf: список, сделанный из мутации в формате vcf (при помощи split())
    :param flank_length: длина участка, фланкирующего мутацию
    :return: строка с регионом, фланкирующем мутацию, в формате bed
    """
    chrom = vcf[0]
    rsid = vcf[1]
    pos = int(vcf[2])
    ref = vcf[3]

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref + '\t' + rsid + '\t' + 'snp'

def get_coords_for_hgvs_g(hgvsg:list, flank_length:int) -> str:
    """
    Эта функция принимает мутацию в формате hgvsg и выписывает координаты региона, фланкирующего мутацию, в формате bed
    :param hgvsg: список, сделанный из мутации в формате hgvsg (при помощи split(':'))
    :param flank_length: длина участка, фланкирующего мутацию
    :return: строка с регионом, фланкирующем мутацию, в формате bed
    """
    chrom = hgvsg[0]
    pos_and_nucl = hgvsg[1].split('.')[1] if '.' in hgvsg[1] else hgvsg[1]
    i = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
    pos = int(pos_and_nucl[:i])
    ref = pos_and_nucl[i:].split('>')[0]

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref + '\t' + ':'.join(hgvsg) + '\t' + 'snp'

def get_pos_c_pos_m_ref(pos_and_nucl:str, sign:str) -> [int, int, str]:
    """
    Эта функция извлекает позицию, смещение позиции и референсный нуклеотид из мутации в формате hgvsc (hgvsc.split()[1])
    :param pos_and_nucl: часть мутации в формате hgvsc, содержащая позиции и референсный нуклеотид
    :param sign: знак, определяющий позицию в интроне относительно экзона (+ или -)
    :return: Список с координатой в экзоне (pos_c), смещением координаты в интрон (pos_m), референсным нуклеотидом (ref)
    """
    pos_c = int(pos_and_nucl.split(sign)[0])
    j = len(['_' for c in range(len(pos_and_nucl.split(sign)[1])) if pos_and_nucl.split(sign)[1][c].isdigit()]) # index of the last digit +1
    ref = pos_and_nucl.split(sign)[1][j:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
        else 'del'+pos_and_nucl.split('del')[1] if 'ins' not in pos_and_nucl \
        else 'ins'+pos_and_nucl.split('ins')[1]
    pos_m = int(pos_and_nucl.split(sign)[1][:j])
    if '_' in str(pos_and_nucl):
        pos_c = int(str(pos_and_nucl).split('_')[0].split(sign)[0])
        pos_m = int(str(pos_and_nucl).split('_')[0].split(sign)[1][:j])
    
    return pos_c, pos_m, ref

def get_pos(strand:str, start_codon_coord:int, exons:list, pos_c:int, ref:str) -> [int, str]:
    """
    Эта функция находит координату в геноме для мутации в формате hgvsc
    :param strand: цепь, на которой находится мутация
    :param start_codon_coord: координата старт-кодона
    :param exons: список с экзонами гена
    :param pos_c: позиция в экзоне
    :param ref: референсный нуклеотид
    :return: Возвращает позицию в геноме и референсный нуклеотид
    """
    i = 0
    if strand == '+':
        while True:
            if int(exons[i][1]) < start_codon_coord:
                i += 1
                continue
            elif int(exons[i][0]) <= start_codon_coord <= int(exons[i][1]):
                exons[i][0] = start_codon_coord
            if pos_c > (int(exons[i][1]) - int(exons[i][0]) + 1):
                pos_c -= (int(exons[i][1]) - int(exons[i][0]) + 1)
                i += 1
            elif pos_c <= (int(exons[i][1]) - int(exons[i][0]) + 1):
                pos = int(exons[i][0]) + pos_c + 2
                break 
    elif strand == '-':
        exons = [[x[1],x[0]] for x in exons]
        ref = rev_dict[ref] if 'del' not in ref and 'ins' not in ref else ref
        while True:
            if int(exons[i][1]) > start_codon_coord:
                i += 1
                continue
            elif int(exons[i][1]) <= start_codon_coord <= int(exons[i][0]):
                exons[i][0] = start_codon_coord
            if pos_c > (int(exons[i][0]) - int(exons[i][1]) + 1):
                pos_c -= (int(exons[i][0]) - int(exons[i][1]) + 1)
                i += 1
            elif pos_c <= (int(exons[i][0]) - int(exons[i][1]) + 1):
                pos = int(exons[i][0]) - pos_c
                break
    return pos, ref

def get_coords_for_hgvs_c(hgvsc:str, flank_length:int, annot:dict, logger:logging) -> str:
    """
    Эта функция принимает мутацию в формате hgvsc и выписывает координаты региона, фланкирующего мутацию, в формате bed
    :param hgvsc: список, сделанный из мутации в формате hgvsc (при помощи split(':'))
    :param flank_length: длина участка, фланкирующего мутацию
    :param annot: словарь с аннотацией генома {транскрипт: [хромосома, [экзоны], стар_кодон, цепь]}
    :param logger: Объект класса logging
    :return: строка с регионом, фланкирующем мутацию, в формате bed
    """
    if hgvsc[0] in annot:
        trans = hgvsc[0]
        chrom = annot[trans][0]
        exons = annot[trans][1]
        start_codon_coord = int(annot[trans][2])
        strand = annot[trans][3]
        pos_and_nucl = hgvsc[1].split('.')[1]
        if pos_and_nucl[0].isdigit():
            if not '-' in pos_and_nucl and not '+' in pos_and_nucl:
                j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
                ref = pos_and_nucl[j:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
                    else 'del'+pos_and_nucl[1:].split('del')[1] if 'ins' not in pos_and_nucl \
                    else 'ins'+pos_and_nucl.split('ins')[1]
                pos_c = int(pos_and_nucl[:j]) - 3
                if '_' in str(pos_and_nucl[:j]):
                    pos_c = int(str(pos_and_nucl[:j]).split('_')[0]) - 3
                pos, ref = get_pos(strand, start_codon_coord, exons, pos_c, ref)

            elif '-' in pos_and_nucl:
                pos_c, pos_m, ref = get_pos_c_pos_m_ref(pos_and_nucl, '-')
                pos, ref = get_pos(strand, start_codon_coord, exons, pos_c, ref)
                pos -= pos_m + 2 #TODO: Пофиксить координаты (пока что точность +- несколько нуклеотидов)

            elif '+' in pos_and_nucl:
                pos_c, pos_m, ref = get_pos_c_pos_m_ref(pos_and_nucl, '+')
                pos, ref = get_pos(strand, start_codon_coord, exons, pos_c, ref)
                pos += pos_m #TODO: Пофиксить координаты (пока что точность +- несколько нуклеотидов)

        elif pos_and_nucl.startswith('-'):
            if not '-' in pos_and_nucl[1:] and not '+' in pos_and_nucl[1:]:
                j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
                ref = pos_and_nucl[j+1:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
                    else 'del'+pos_and_nucl[1:].split('del')[1] if 'ins' not in pos_and_nucl \
                    else 'ins'+pos_and_nucl.split('ins')[1]
                if '_' in str(pos_and_nucl):
                    pos_c = int(str(pos_and_nucl[1:]).split('_')[0])
                else:
                    pos_c = int(pos_and_nucl[1:j+1])
                pos = int(exons[0][0]) - pos_c if strand == '+' else exons[0][1] + pos_c
            elif '-' in pos_and_nucl[1:]:
                pos_c, pos_m, ref = get_pos_c_pos_m_ref(pos_and_nucl[1:], '-')
                pos = int(exons[0][0]) - pos_c - pos_m if strand == '+' else exons[0][1] + pos_c - pos_m
            elif '+' in pos_and_nucl[1:]:
                pos_c, pos_m, ref = get_pos_c_pos_m_ref(pos_and_nucl[1:], '+')
                pos = int(exons[0][0]) - pos_c + pos_m if strand == '+' else exons[0][1] + pos_c - pos_m

        elif pos_and_nucl.startswith('*'):
            if not '-' in pos_and_nucl[1:] and not '+' in pos_and_nucl[1:]:
                j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
                ref = pos_and_nucl[j+1:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
                    else 'del'+pos_and_nucl[1:].split('del')[1] if 'ins' not in pos_and_nucl \
                    else 'ins'+pos_and_nucl.split('ins')[1]
                if '_' in str(pos_and_nucl):
                    pos_c = int(str(pos_and_nucl[1:]).split('_')[0])
                else:
                    pos_c = int(pos_and_nucl[1:j+1])
                pos = int(exons[-1][1]) + pos_c if strand == '+' else exons[-1][0] - pos_c
            elif '-' in pos_and_nucl[1:]:
                pos_c, pos_m, ref = get_pos_c_pos_m_ref(pos_and_nucl[1:], '-')
                pos = int(exons[-1][1]) + pos_c - pos_m if strand == '+' else exons[-1][0] - pos_c - pos_m
            elif '+' in pos_and_nucl[1:]:
                pos_c, pos_m, ref = get_pos_c_pos_m_ref(pos_and_nucl[1:], '+')
                pos = int(exons[-1][1]) + pos_c + pos_m if strand == '+' else exons[-1][0] - pos_c + pos_m
        else:
            logger.warning(f'This format of HGVS mutation ({trans +":"+hgvsc[1]}) is not supproted yet. Skipping this variant')
            return None
    else:
        logger.warning(f'Transcript {hgvsc[0]} not present in provided annotation!!! Skipping {hgvsc[0] +":"+hgvsc[1]} variant!')
        return None

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref + '\t' + ':'.join(hgvsc) + '\t' + 'snp'

def fix_chromosome_name_and_check_coords(coord_list:list, genome_version:str, logger:logging) -> list:
    """
    Эта функция проверяет список регионов и фиксит при необходимости названия хромосом дял человека
    (приводит все хромосомы к названию типа chr1)
    :param coord_list: Список координат регионов, фланкирующих мутации
    :param genome_version: Версия генома (GRCh38/37)
    :param logger: Объект класса logging
    :return: Исправленный и проверенный список с координатами в формате bed
    """
    fixed_coord_list = []
    for coord in coord_list:
        coord = coord.split()
        chr_in_input = coord[0]
        if coord[0][:3] != 'chr':
            if coord[0].isdigit() or coord[0] == 'X' or coord[0] == 'Y':
                fixed_chr = 'chr' + str(coord[0])
            elif coord[0] in genome_versions_chromosome_dict[genome_version].keys():
                fixed_chr = genome_versions_chromosome_dict[genome_version][coord[0]]
            else:
                logger.warning(f"Entry {coord} has invalid chromosome identifier '{coord[0]}' !!!!")
                logger.warning(f"For genome {genome_version} only the following identifier formats are available : 1, chr1, NC_000001.11")
                logger.warning(f"Skipping this variant")
                continue
        else:
            fixed_chr = coord[0]
        
        if int(coord[1]) < 0 or int(coord[2]) < 0:
            continue

        try:
            fixed_coord_list.append('\t'.join([fixed_chr, coord[1], coord[2], coord[3], coord[4], coord[5], chr_in_input]))
        except:
            logger.warning(f'Something wrong with {coord} variant! It has no ref allele! Skipping this variant!')
    
    return fixed_coord_list

def get_regions_for_muts(muts_list:list, annot:dict, flank_length:int, write_intermediate:bool, tmpdir:Path,
                         logger:logging, genome_version:str) -> list: 
    """
    Эта функция выполняет поиск регионов, на которые будут подбираться олигонуклеотиды, для поданных на вход мутаций.
    :param muts_list: Список мутаций, которые были поданы на вход (vcf, hgvsg или hgvsc формат)
    :params annot: аннотация генома
    :param flank_length: длина фланка вокруг мутации, который будет включен в регион
    :param write_intermediate: Сохранять или нет промежуточный файл с регионами
    :param tmpdir: Путь к промежуточной директории (должна быть создана)
    :param logger: Объект класса logger, который будет сохранять логи
    :param genomve_version: Версия генома. Используется, чтобы фиксить номера хромосом при необходимости.
    :return: Список регионов для мутаций в формате bed
    """
    bed_coords_list = []
    vcf_counter, hgvs_counter = 0, 0
    for mut in muts_list:
        if mut.startswith('#'):
            continue
        if len(mut.split('\t')) > 1 and mut.split('\t')[2].isdigit():
            # If there is tab-separated values, this is vcf-like mutations
            vcf_counter += 1
            bed_coords = get_coords_for_vcf(mut.split('\t'), flank_length)
            bed_coords_list.append(bed_coords)
        elif (len(mut.split(':')) > 1) and ('g.' in mut.split(':')[1] or '.' not in mut.split(':')[1]):
            hgvs_counter += 1
            bed_coords = get_coords_for_hgvs_g(mut.split(':'), flank_length)
            bed_coords_list.append(bed_coords)
        elif len(mut.split(':')) > 1 and ('c.') in mut.split(':')[1]:
            hgvs_counter += 1
            bed_coords = get_coords_for_hgvs_c(mut.split(':'), flank_length, annot, logger)
            if bed_coords:
                bed_coords_list.append(bed_coords)
    
    if vcf_counter:
        logger.info(f'There were found {vcf_counter} vcf-formatted mutations in input')
    if hgvs_counter:
        logger.info(f'There were found {hgvs_counter} hgvs-formatted mutations in input')

    if not any([vcf_counter, hgvs_counter]):
        logger.error('ERROR! Mutations in appropriate format (vcf or hgvs) WERE NOT FOUND!')
        logger.error('ERROR! Check input files! Exiting...')
        exit(1)

    # fix chromosome names for human genome assemblies if needed
    if genome_version in ['GRCh38', 'GRCh37']:
        bed_coords_list = fix_chromosome_name_and_check_coords(bed_coords_list, genome_version, logger)

    if not bed_coords_list:
        logger.error('ERROR! No mutation were found! Check input files! Exiting...')
        exit(1)

    mut_found_percent = len(bed_coords_list)/(vcf_counter + hgvs_counter)
    logger.info(f'Genome coordinates were found for {len(bed_coords_list)} ({round(mut_found_percent*100)}%) mutations')
    if mut_found_percent < 0.8:
        logger.warning(f'Too little input mutations have coordinates in provided genome assembly! Please check if you use right genome assmbly')

    if write_intermediate:
        with open(tmpdir/'INPUT_coords_for_mutations.bed', 'w') as outfile:
            logger.info(f"Coordinates for input mutations were written in {tmpdir/'INPUT_coords_for_mutations.bed'} file")
            outfile.write('\n'.join(bed_coords_list))

    return bed_coords_list


if __name__ == '__main__':
    annot = {'test':['1', [[10,90],[100,200],[400,500]], 12, '+']}
    hgvsc = ['test','c.47A>T']
    flank_length = 10
   
    print(get_coords_for_hgvs_c(hgvsc, flank_length, annot, ""))