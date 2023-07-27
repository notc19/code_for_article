import logging
from pathlib import Path

from constants import genome_versions_chromosome_dict
from get_regions_for_mutations import fix_chromosome_name_and_check_coords

def get_regions_for_genes(input_genes_list:list, annot_for_genes:dict, gene_name_dict:dict,
                          gene_regions:str, intron_depth:int, write_intermediate:bool, tmpdir:Path, 
                          logger:logging, genome_version:str) -> list:
    """
    Эта функция находит геномные координаты для каждого гена в инпуте (учитывая необходимые регионы и фланки)
    :param input_genes_list: Список генов, для которых необходимо найти координаты
    :param annot_for_genes: Словарь с аннотацией генов
    :param gene_name_dict: Словарь с соответствием названий и идентификаторов генов
    :param gene_regions: Какой тип регионов по умолчанию выписывать для каждого гена (exons/exons_utr/full)
    :param intron_depth: Какую глубину интрона по умолчанию включать в регион для каждого гена
    :param write_intermediate: Сохранять или нет промежуточный файл с регионами в формате bed
    :param tmpdir: Объект Path с путем до папки с промежуточными результатами
    :param logger: Объект класса logging
    :param genome_version: Версия генома (нужна для того, чтобы фиксить названия хромосом для человека)
    :return: Возвращает список регионов в формате bed
    """
    gene_coord_list = []

    for gene in input_genes_list:
        genesplit = gene.split()
        gene_name = genesplit[0]
        if len(genesplit) >= 2:
            if genesplit[1]:
                gene_regions = genesplit[1]
        if len(genesplit) == 3:
            intron_depth = int(genesplit[2].strip())
        if len(genesplit) > 3:
            logger.warning(f'There is too many columns in input file for gene {gene_name}!!! Skipping extra columns {genesplit[3:]}')
        
        if gene_name in gene_name_dict:
            gene_name = gene_name_dict[gene_name]
        elif gene_name not in gene_name_dict.values():
            logger.warning(f'Gene {gene_name} not present in provided annotatin!!! Skipping this gene!')
            continue

        gene_annot = annot_for_genes[gene_name]
        gene_chr = gene_annot[0]

        if gene_regions == 'full':
            gene_start = int(gene_annot[1][0]) - intron_depth
            gene_end = int(gene_annot[1][1]) + intron_depth
            gene_coord_list.append(f'{gene_chr}\t{str(gene_start)}\t{str(gene_end)}\t{gene_name}\t{gene_regions}\t{gene_regions}')
        elif gene_regions == 'exon' or gene_regions == 'exon_utr':
            for i in range(len(gene_annot[2])):
                gene_start = int(gene_annot[2][i][0]) - intron_depth
                gene_end   = int(gene_annot[2][i][1]) + intron_depth
                gene_coord_list.append(f'{gene_chr}\t{str(gene_start)}\t{str(gene_end)}\t{gene_name}\t{gene_annot[4][i]}\t{gene_regions}') #
        elif gene_regions == 'exon_utr':
            for i in range(len(gene_annot[3])):
                gene_start = int(gene_annot[3][i][0]) - intron_depth
                gene_end   = int(gene_annot[3][i][1]) + intron_depth
                gene_coord_list.append(f'{gene_chr}\t{str(gene_start)}\t{str(gene_end)}\t{gene_name}\t{gene_regions[5:]}_{i+1}\t{gene_regions}')
        elif gene_regions.isdigit() or '-' in gene_regions or ',' in gene_regions:
            for c in gene_regions.split(','):
                if c.isdigit():
                    gene_start = int(gene_annot[2][int(c)][0]) - intron_depth
                    gene_end   = int(gene_annot[2][int(c)][1]) + intron_depth
                    gene_coord_list.append(f'{gene_chr}\t{str(gene_start)}\t{str(gene_end)}\t{gene_name}\t{gene_annot[4][int(c)-1]}\t{c}')
                elif '-' in c:
                    for j in range(int(c.split('-')[0])-1, int(c.split('-')[1])):
                        gene_start = int(gene_annot[2][j][0]) - intron_depth
                        gene_end   = int(gene_annot[2][j][1]) + intron_depth
                        gene_coord_list.append(f'{gene_chr}\t{str(gene_start)}\t{str(gene_end)}\t{gene_name}\t{gene_annot[4][j]}\t{c}')
                elif c.isalnum():
                    gene_start = int(gene_annot[2][gene_annot[4].index(c)][0]) - intron_depth
                    gene_end   = int(gene_annot[2][gene_annot[4].index(c)][1]) + intron_depth
                    gene_coord_list.append(f'{gene_chr}\t{str(gene_start)}\t{str(gene_end)}\t{gene_name}\t{c}\texon')
                else:
                    logger.warning(f'WARNING! Unknown region type {gene_regions}! Skipping this gene!')
                    continue
        else:
            logger.warning(f'WARNING! Unknown region type {gene_regions}! Skipping this gene!')
            continue
    
    if genome_version in ['GRCh38', 'GRCh37']:
        gene_coord_list = fix_chromosome_name_and_check_coords(gene_coord_list, genome_version, logger)

    if not gene_coord_list:
        logger.error('ERROR! No genes were found! Check input files! Exiting...')
        exit(1)

    ### TODO: Merge intersecting regions and drop nested regions!

    if write_intermediate:
        with open(tmpdir/'INPUT_coords_for_genes.bed', 'w') as outfile:
            logger.info(f"Coordinates for input genes were written in {tmpdir/'INPUT_coords_for_genes.bed'} file")
            outfile.write('\n'.join(gene_coord_list))

    logger.info(f'{len(gene_coord_list)} coordinates were found for {len(input_genes_list)} genes')

    return gene_coord_list

if __name__ == '__main__':
    pass