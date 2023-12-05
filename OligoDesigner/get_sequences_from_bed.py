import logging
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from constants import genome_versions_chromosome_dict

def create_fasta_sequences(coords_list:list, chromosome: str, logger:logging) -> list[SeqRecord]:
    """
    Эта функция создает объект SeqRecord с последовательностью региона.
    :param coords_list: список координат с данной хромосомы, для которых будет найдена последовательность
    :param chromosome: последовательность хромосомы в одну строку
    :param logger: объект класса logging
    :return: список объектов класса SeqRecord
    """
    sequence_list = []
    for coords in coords_list:
        coords = coords.split()
        start = int(coords[1])
        end = int(coords[2])
        if len(coords) == 7:
            if coords[5] == 'snp':
                description = f'SNP chrom:{coords[0]} start:{coords[1]} end:{coords[2]} ref_allele:{coords[3]} snp_id:{coords[4]}'
            else:
                description = f'Region chrom:{coords[0]} start:{start} end:{end} gene:{coords[3]} exon:{coords[4]}'
        else:
            description = f'{coords}'
        sequence = chromosome[start-1:end]
        if len(coords) > 5 and coords[5] == 'snp':
            if (coords[3].upper().count('A') + coords[3].upper().count('T') + coords[3].upper().count('G') \
                + coords[3].upper().count('C')) == len(coords[3]):
                if not (allele:=sequence[int(len(sequence)/2)]) == coords[3]:
                    logger.warning(f'For snp {coords[4]} ref allele in snp ({coords[3]}) not equal to allele in provided reference genome ({allele}) in position {int(start + (len(sequence)/2))}. SKIPPING this region!')
                    continue

        record_id = coords[4] if len(coords) > 4 else 'ID'
        record_name = 'snp' if len(coords) > 5 and coords[5] == 'snp' else coords[3] if len(coords) > 3 else 'NAME'
        record = SeqRecord(
            Seq(sequence),
            id=record_id,
            name=record_name,
            description=description
        )
        sequence_list.append(record)

    return sequence_list

def get_fasta_seqs_from_bed(coord_list:list, genome:Path, genome_version:str, write_intermediate:bool, 
                            tmpdir:Path, logger:logging) -> list[SeqRecord]:
    """
    Эта функция создает fasta последовательности для поданных на вход регионов.
    :param coord_list: Список координат, для которых нужно найти последовательности
    :param genome: путь к файлу с геномом
    :genome_version: версия генома. Нужна для только, чтобы фиксить названия хромосом в геноме человека
    :write intermediate: Сохранять ли файл с последовательностями
    :tmpdir: Директория, в которую будет сохранен файл с последовательностями
    :logger: Объект класса logging
    :return: Список с последовательностями в виде объекта класса SeqRecord
    """
    chromosomes = list(set([i.split()[0] for i in coord_list]))
    chrom_coord_dict = {}
    for chrom in chromosomes:
        chrom_coord_dict[chrom] = [i for i in coord_list if i.split()[0] == chrom]

    logger.info("Start extracting fasta sequence from genome")
    logger.info(f"Reading {genome} file")
    chr_str = ', '.join(sorted(chromosomes))
    logger.info(f"Regions were found on chromosomes {chr_str}")
    fasta_seqs_list = []
    with open(genome) as infile:
        flag = 0
        while True:
            line = infile.readline().strip()
            if not line:
                break
            if line.startswith('>'):
                if flag:
                    logger.info(f'Reading chromosome {chromosome}')
                    fasta_seqs_list.extend(create_fasta_sequences(chrom_coord_dict[chromosome], chrom_seq, logger))
                    chrom_seq = ''
                    flag = 0
                chromosome = line.split('>')[1].split(' ')[0]
                if genome_version in ['GRCh38', 'GRCh37']:
                    if chromosome.isdigit() or chromosome == 'X' or chromosome == 'Y' or chromosome == 'M':
                        chromosome = 'chr' + chromosome
                    elif chromosome in genome_versions_chromosome_dict[genome_version]:
                        chromosome = genome_versions_chromosome_dict[chromosome]
                if chromosome in chromosomes:
                    chrom_seq = ''
                    flag = 1
                else:
                    flag = 0
            elif flag:
                chrom_seq += line

    if write_intermediate:
        with open(tmpdir/'INPUT_regions_sequences.fasta', 'w') as outfile:
            logger.info(f"Sequences for input regions were written in {tmpdir/'INPUT_regions_sequences.fasta'} file")
            SeqIO.write(fasta_seqs_list, outfile, "fasta")

    return fasta_seqs_list


if __name__ == '__main__':
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)
    coord_list = ['chr1\t1471665\t1472189\tENSG00000160072\tENSE00003889014\t1-2\t1', 'chr1\t1477174\t1477450\tENSG00000160072\tENSE00003467707\t1-2\t1',\
            'chr1\t1478949\t1479208\tENSG00000160072\tENSE00003608502\t4-5\t1', 'chr1\t1480767\t1481036\tENSG00000160072\tENSE00003474888\t4-5\t1',\
            'chr1\t10\t90\ttest_g\texo1\texon\t1', 'chr1\t100\t200\ttest_g\texo2\texon\t1',\
            'chr1\t400\t500\ttest_g\texo3\texon\t1']
    #genome = '/home/a_vikhorev/Data/DataBases/Hg38_GCA_000001405.15_BWA/GCA_000001405.15.fasta'
    genome = '/home/a_vikhorev/OligoDesigner/tests/FileTest_genome.fasta'
    coord_list = ['chr1\t4\t15\ttest_g\texo3\texon\t1', 'chr3\t10\t20\ttest_g\texo3\texon\t1']
    import time
    start = time.time() ## точка отсчета времени
    f = get_fasta_seqs_from_bed(coord_list, genome, 'GRCh38.p14', True, './', logger)
    print(f)
    end = time.time() - start
    print(end)
