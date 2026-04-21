from logging import Logger
from pathlib import Path

from constants import genome_versions_chromosome_dict
from get_regions_for_mutations import fix_chromosome_name_and_check_coords
import merge_and_filter_coordinates as mafc

def get_regions_for_genes(input_genes_list:list, annot_for_genes:dict, gene_name_dict:dict,
                          gene_regions:str, intron_depth:int, write_intermediate:bool, tmpdir:Path, 
                          logger:Logger, genome_version:str, length:int) -> list:
    """
    Эта функция находит геномные координаты для каждого гена в input (учитывая необходимые регионы и фланки)
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
        current_regions = gene_regions #param 'gene_regions' was overwritten inside the loop! local copy will not mutate 'gene_regions'
        current_depth   = intron_depth
        genesplit = gene.split()
        gene_name = genesplit[0]
        if len(genesplit) == 2 and genesplit[1]:
            current_regions = genesplit[1]
        if len(genesplit) == 3:
            current_regions = genesplit[1]
            current_depth = int(genesplit[2].strip())
        if len(genesplit) > 3:
            current_regions = genesplit[1]
            current_depth = int(genesplit[2].strip())
            logger.warning(f'There are too many columns in input file for gene {gene_name}! Skipping extra columns {genesplit[3:]}')
        
        if gene_name in gene_name_dict:
            logger.info(f'Data for {gene_name} was found')
            gene_name = gene_name_dict[gene_name]
        elif gene_name not in gene_name_dict.values():
            logger.warning(f'Gene {gene_name} is not present in provided annotation! Skipping this gene!')
            continue

        gene_annot = annot_for_genes[gene_name]
        gene_chr = gene_annot[0]
        
        if current_regions == 'full':
            gene_start = int(gene_annot[1][0]) - current_depth
            gene_end   = int(gene_annot[1][1]) + current_depth
            gene_coord_list.append(f'{gene_chr}\t{gene_start}\t{gene_end}\t{gene_name}\t{gene_annot[8]}\t{current_regions}')
        elif current_regions in ('exon', 'exon_utr'):
            for i, (exon, exon_id) in enumerate(zip(gene_annot[2], gene_annot[4])):
                gene_start = exon[0] - current_depth
                gene_end   = exon[1] + current_depth
                exon_length = gene_end - gene_start
                if exon_length < length:
                    diff = length - exon_length
                    gene_start = gene_start - int(diff/2) - 1
                    gene_end = gene_end + int(diff/2) + 1
                gene_coord_list.append(f'{gene_chr}\t{gene_start}\t{gene_end}\t{gene_name}\t{exon_id}\t{current_regions}')
        elif current_regions == 'exon_utr':
            for i in range(len(gene_annot[3])):
                gene_start = int(gene_annot[3][i][0]) - current_depth
                gene_end   = int(gene_annot[3][i][1]) + current_depth
                gene_coord_list.append(f'{gene_chr}\t{gene_start}\t{gene_end}\t{gene_name}\t{gene_annot[3]}_{i+1}\t{current_regions}')
        elif current_regions.isdigit() or '-' in current_regions or ',' in current_regions:
            for c in current_regions.split(','):
                if c.isdigit():
                    idx = int(c)-1
                    gene_start = int(gene_annot[2][idx][0]) - current_depth
                    gene_end   = int(gene_annot[2][idx][1]) + current_depth
                    gene_coord_list.append(f'{gene_chr}\t{gene_start}\t{gene_end}\t{gene_name}\t{gene_annot[4][idx]}\t{c}')
                elif '-' in c:
                    for j in range(int(c.split('-')[0])-1, int(c.split('-')[1])):
                        gene_start = int(gene_annot[2][j][0]) - current_depth
                        gene_end   = int(gene_annot[2][j][1]) + current_depth
                        gene_coord_list.append(f'{gene_chr}\t{gene_start}\t{gene_end}\t{gene_name}\t{gene_annot[4][j]}\t{c}')
                elif c.isalnum():
                    gene_start = int(gene_annot[2][gene_annot[4].index(c)][0]) - current_depth
                    gene_end   = int(gene_annot[2][gene_annot[4].index(c)][1]) + current_depth
                    gene_coord_list.append(f'{gene_chr}\t{gene_start}\t{gene_end}\t{gene_name}\t{c}\texon')
                else:
                    logger.warning(f'WARNING! Unknown region type {gene}, {current_regions}! Skipping this gene!')
                    continue
        else:
            logger.warning(f'WARNING! Unable to parse region type for {gene}, {current_regions}! Skipping this gene!')
            continue
    
    if genome_version in ['GRCh38', 'GRCh37']:
        gene_coord_list = fix_chromosome_name_and_check_coords(gene_coord_list, genome_version, logger)
    else:
        logger.error('ERROR! Genome version is not supported. Exiting...')
        exit(1)

    if not gene_coord_list:
        logger.error('ERROR! No genes were found! Check input files! Exiting...')

    #mafc.process_coordinates is called in oligodesigner.py
    #gene_coord_list = mafc.process_coordinates(gene_coord_list, True, tmpdir, logger)

    if write_intermediate:
        with open(f'{tmpdir}/INPUT_coords_for_genes.bed', 'w') as outfile:
            logger.info(f"Coordinates for input genes were written in {tmpdir}/INPUT_coords_for_genes.bed file")
            outfile.write('\n'.join(gene_coord_list))

    logger.info(f'{len(gene_coord_list)} coordinates were found for {len(input_genes_list)} genes')

    return gene_coord_list

if __name__ == '__main__':
    import utils
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)

    input_genes_file = '../tests/FileTest_cardio_genes.txt'
    annot_file = '../tests/FileTest_cardio_gencode_annot.gtf'
    gene_regions = 'exon'
    intron_depth = 50
    length = 100
    write_intermediate = True
    tmpdir = '../tests/test_results'
    genome_version = 'GRCh38'

    with open(input_genes_file) as infile:
        logger.info(f'The File with genes: {input_genes_file}')
        input_genes_list = [line.strip() for line in infile.readlines()]

    with open(annot_file) as infile:
        annot = [line.strip() for line in infile.readlines() if not line.startswith('#')]
        if input_genes_list:
            annot_for_genes, gene_name_dict = utils.prepare_annot_for_genes(annot, logger)

    input_genes_bed_list = get_regions_for_genes(input_genes_list, annot_for_genes, gene_name_dict,
                                                gene_regions, intron_depth,
                                                write_intermediate, tmpdir, logger,
                                                genome_version, length)
    print(input_genes_bed_list)
