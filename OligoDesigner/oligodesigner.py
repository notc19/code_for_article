###########################################
# imports
###########################################

import os
import logging
from logging import Logger
from datetime import datetime
import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import DNAFMIndex

import get_regions_for_mutations as grfm
import get_regions_for_genes as grfg
import get_sequences_from_bed as gsfb
import dna_fm_index as dfmi
import download_annot_and_genome as dag
import merge_and_filter_coordinates as mafc
import merge_and_filter_sequences as mafs
import run_oligo_design as rod
import run_oligo_postprocessing as rop
import oligonucleotide_analysis as olan
import utils
import models

def create_logger(timestamp: datetime):
    """
    Функция создает logger и хэндлеры, которые передаются в listener в модуле run_oligo_design.
    :param timestamp: время запуска программы
    return: объект класса Logger и Logger handlers
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

def parse_args(
    default_input_genes: str = None,
    default_input_muts: str = None,
    default_input_bed: str = None,
    default_input_fasta: str = None,
    default_output: str = None,
    default_tmpdir: str = None,
    default_annot: str = None,
    default_genome: str = None,
    default_genome_index: str = None,
) -> argparse.Namespace:
    """
    Эта функция парсит аргументы, прнимаемые на вход программой, и возвращает их в виде переменной args_parser.parse_args()
    """
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)   
    args_parser.add_argument('-i','--input_genes',help='path to input file with genes',
                             type=str,required=False, default=default_input_genes)
    args_parser.add_argument('-m','--input_muts', help='path to input file with mutations',
                             type=str,required=False, default=default_input_muts)
    args_parser.add_argument('-b','--input_bed', help='path to input file with target regions in bed format',
                             type=str,required=False, default=default_input_bed)
    args_parser.add_argument('-f','--input_fasta', help='path to input file with target sequences in fasta format',
                             type=str,required=False, default=default_input_fasta)
    
    args_parser.add_argument('--threads', help='Number of threads',
                             type=int, required=False, default=10)    
    args_parser.add_argument('-o','--output', help='path of output folder',
                             type=str,required=False, default=default_output)
    args_parser.add_argument('-t','--tmpdir', help='path to directory with temporary files',
                             type=str,required=False,default=default_tmpdir)
    args_parser.add_argument('-w','--write_intermediate', help='write intermediate results or not',
                             action=argparse.BooleanOptionalAction, default=True)                             

    args_parser.add_argument('-a','--annot', help='path to file with genome annotation in GTF format',
                             type=str, required=False,default=default_annot)
    args_parser.add_argument('--annotation_release', help='release of the reference genome',
                             type=str, required=False, default='current')    
    args_parser.add_argument('-g','--genome', help='path to file with reference genome in fasta format ',
                             type=str, required=False, default=default_genome)
    args_parser.add_argument('--genome_index', help='path to file with genome index in awfmi format (build by DNAFMIndex)',
                             type=str, required=False, default=default_genome_index)
    args_parser.add_argument('-v','--genome_version', help='version of the reference genome',
                             type=str, required=False, default='GRCh38')
    args_parser.add_argument('-s','--species', help='Analyzed species',
                             type=str, required=False, default='homo_sapiens')

    args_parser.add_argument('--download_annotation', help='Download annotation file from Ensembl or not',
                             action=argparse.BooleanOptionalAction, default=False)
    args_parser.add_argument('--download_genome', help='Download genome file from Ensembl or not',
                             action=argparse.BooleanOptionalAction, default=False)    

    args_parser.add_argument('-l','--length', help='Length of the designed olignucleotide bait',
                             type=int, required=False, default=100)
    args_parser.add_argument('--max_i', help='Maximum intersection of oligos in region',
                             type=int, required=False, default=50)
    args_parser.add_argument('--min_i', help='Minimum intersection of oligos in region',
                             type=int, required=False, default=-20)
    args_parser.add_argument('-k','--mutation_flank', help='Flank length from mutation',
                             type=int, required=False, default=100)
    args_parser.add_argument('-r','--gene_regions', help='Type of gene regions to use [full/exons/exons_utr]',
                             type=str, required=False, default='exon')
    args_parser.add_argument('-j','--intron_depth', help='Depth of introns to include when only exons regions were chosen',
                             type=int, required=False, default=50)

    args_parser.add_argument('--process_coordinates', help='Whether to process coordinates to remove duplicate and nested ones and merge intersecting ones',
                             action=argparse.BooleanOptionalAction, default=True)
    args_parser.add_argument('--process_sequences', help='Whether to process sequences to cluster highly identical ones',
                             action=argparse.BooleanOptionalAction, default=False)
    args_parser.add_argument('--target_max_id', help='Maximum allowed identity for input sequences. Sequences with greater identity will be clustered using vsearch',
                             type=float, required=False, default=0.999)
    args_parser.add_argument('--vsearch_binary', help='Absolute path to vsearch binary file',
                             type=str, required=False, default='')
    args_parser.add_argument('--analytics', help='Whether to perform the analysis of the designed oligonucleotides',
                             action=argparse.BooleanOptionalAction, default=True)
    args_parser.add_argument('--postprocess_oligos', help='Whether to perform postprocessing of the designed oligonucleotides',
                             action=argparse.BooleanOptionalAction, default=False) 
    
    args_parser.add_argument('--predict_blast_offtarget', help='Whether to predict off-target Nbases based on BLAST results',
                             action=argparse.BooleanOptionalAction, default=False)
    args_parser.add_argument('--blast_binary', help='Path to blastn binary',
                             type=str, required=False, default='/home/a_makarova/miniconda3/envs/od/bin/blastn')
    args_parser.add_argument('--blast_db', help='Path to BLAST database',
                             type=str, required=False, default='/home/a_makarova/Data/blastdb_GCA/hg38')
    return args_parser.parse_args()

def parse_input(inpt:str, logger:Logger, feature:str) -> list | SeqRecord:
    """
    Эта функция парсит подаваемый на вход файл с мутациями/генами/регионами/последовательностями. Возвращает или список
    с фичами (если на вход подан файл с мутациями/генами/регионами) или объект SeqRecord с последовательностями
    :param inpt: Путь до входного файла
    :param logger: Объект класса Logger
    :param feature: описание входного файла
    """
    try:
        if feature != 'sequences in fasta':
            with open(inpt) as infile:
                logger.info(f'The File with input {feature}: {inpt}')
                entries = [line.strip() for line in infile]
                logger.info(f'{len(entries)} entries were found in {inpt}')
                assert len(entries) > 0, (
                    f'File with input {feature} was provided, but it is empty!'
                )
                return entries
        else:
            logger.info(f'The File with input {feature}: {inpt}')
            try:
                entries = list(SeqIO.parse(inpt, 'fasta'))
                logger.info(f'{len(entries)} entries were found in {inpt}')
                assert len(entries) > 0, (
                    f'File with input {feature} was provided, but it is empty!'
                )
                return entries
            except Exception as e:
                logger.error(f'Something went wrong reading the FASTA file: {e}')
                exit(1)
    except FileNotFoundError:
        logger.error(f'File {inpt} does not exist!')
        return False

def check_dir(dir_path:str, logger:Logger, dir_name:str) -> Path:
    """
    Эта функция проверяет, существует ли директория. Если не существует, то создает ее.
    :param dir_path: путь до директории
    :param logger: объект класса Logger
    :param dir_name: название директории
    """
    if os.path.exists(dir_path):
        logger.info(f'Directory {dir_name}: {dir_path} exists')
    else:
        try:
            os.mkdir(dir_path)
            logger.info(f'The directory {dir_name}: {dir_path} was created!')
        except OSError as e:
            logger.error(f'Cannot create {dir_path} directory: {e}')
            raise

    return Path(dir_path)

def check_file(file:str, logger:Logger, file_name:str) -> Path:
    """
    Эта функция проверяет, существует ли файл. Если не существует, то создает его.
    :param file: путь к файлу
    :param logger: объект класса Logger
    :param file_name: название файла
    """
    if os.path.exists(file):
        logger.info(f'The {file_name} {file} was found!')
    else:
        logger.error(f'ERROR! The {file_name} {file} DOES NOT EXIST!')
        raise FileNotFoundError(file)

    return Path(file)

def main():
    """
    Функция, которая запускает программу!
    """
    # Setup logging
    timestamp = datetime.now()
    logger, log_handlers = create_logger(timestamp)
    logger.info(f'######### Start OligoDesigner at {timestamp.day}.{timestamp.month}.{timestamp.year} {timestamp.hour}:{timestamp.minute}:{timestamp.second} #########')

    # Read params
    params = parse_args()
    if params.analytics and not params.process_coordinates:
        logger.info("Parameter --process_coordinates should be 'True' to perform the panel analysis")
        exit(1)
    logger.info(f'OligoDesigner was launched for {params.gene_regions} regions, intron depth={params.intron_depth}.')

    # Read input files
    input_genes_list = parse_input(params.input_genes, logger, 'genes') if params.input_genes else []
    input_muts_list  = parse_input(params.input_muts, logger, 'mutations') if params.input_muts else []
    input_bed_list   = parse_input(params.input_bed, logger, 'coordinates in bed format') if params.input_bed else []
    input_fasta_list = parse_input(params.input_fasta, logger, 'sequences in fasta') if params.input_fasta else []

    # Terminate program if there is no input files
    if not any([input_genes_list, input_muts_list, input_bed_list, input_fasta_list]):
        logger.error('There are no input files! Exiting...')
        exit(1)

    # Check output directories
    outdir = check_dir(params.output, logger, 'for output files')
    tmpdir = check_dir(params.tmpdir, logger, 'for temporary files')

    # Prepare annotation of the reference genome
    if input_genes_list or input_muts_list:
        if not params.download_annotation:
            annot_file = check_file(params.annot, logger, 'GTF file with annotation')
        else:
            logger.info('Downloading annotation from Ensembl')
            annot_file = dag.download_from_ensembl(
                'annotation', params.species, params.annotation_release,
                params.genome_version, tmpdir, logger,
            )
            logger.info(f'Annotation file was downloaded! {annot_file}')
        with open(annot_file) as infile:
            logger.info(f'Processing annotation file {annot_file}')
            annot = [line.strip() for line in infile if not line.startswith('#')]
            if input_genes_list:
                annot_for_genes, gene_name_dict = utils.prepare_annot_for_genes(annot, logger)
            if input_muts_list:
                annot_for_muts = utils.prepare_annot_for_muts(input_muts_list, annot, logger)
            logger.info('Annotation file was processed')

    # Prepare reference genome
    if input_genes_list or input_muts_list or input_bed_list:
        if not params.download_genome:
            genome_file = check_file(params.genome, logger, 'FASTA file with reference genome')
        else:
            logger.info('Downloading genome from Ensembl')
            genome_file = dag.download_from_ensembl(
                'genome', params.species, params.annotation_release,
                params.genome_version, tmpdir, logger,
            )
            logger.info('Genome was downloaded')

    if params.genome_index:
        check_file(params.genome_index, logger, 'file with genome FM-index')
    else:
        logger.info(f'Creating FM-index from genome {params.genome}.')
        dfmi.create_fm_index_from_genome(params.genome, None)
        logger.info(f'FM-index {params.genome}.awfmi was created from genome.')
        params.genome_index = params.genome+'.awfmi'

    # Get BED regions for input genes
    if input_genes_list:
        input_genes_bed_list = grfg.get_regions_for_genes(
            input_genes_list, annot_for_genes, gene_name_dict,
            params.gene_regions, params.intron_depth,
            params.write_intermediate, tmpdir, logger,
            params.genome_version, params.length,
        )
        input_bed_list.extend(input_genes_bed_list)

    # Get BED regions for input mutations
    if input_muts_list:
        flank_length = int(params.length * 1.7 / 2) if not params.mutation_flank else params.mutation_flank
        input_muts_bed_list = grfm.get_regions_for_muts(
            input_muts_list, annot_for_muts, flank_length,
            params.write_intermediate, tmpdir, logger,
            params.genome_version,
        )
        input_bed_list.extend(input_muts_bed_list)

    # Fix coordinates: merge intersecting, filter duplicated and nested
    if params.process_coordinates:
        input_bed_list = mafc.process_coordinates(input_bed_list, params.write_intermediate, tmpdir, logger)
        logger.info(f'Length of processed input_bed_list: {len(input_bed_list)}')

    # Get FASTA sequences for input regions
    if input_bed_list:
        input_fasta_intermed_list = gsfb.get_fasta_seqs_from_bed(
            input_bed_list, genome_file, params.genome_version,
            params.write_intermediate, tmpdir, logger,
        )
        input_fasta_list.extend(input_fasta_intermed_list)

    # Fix sequences: merge intersecting and cluster highly identical ones
    #input_fasta_list = mafs.process_sequences(input_fasta_list, outdir, logger, params.vsearch_binary, params.target_max_id)

    # Create oligonucleotides!
    oligos = rod.run_oligo_design(input_fasta_list, params.threads, params, logger, log_handlers)

    seq_to_write_list = [SeqRecord(Seq(olig[-4]), id = f'olig_{c+1}', description=' '.join(olig[-3:])) \
                         for c, olig in enumerate([o.split(' ') for o in list(map(str,oligos))])]
    SeqIO.write(seq_to_write_list, f'{tmpdir}/OUTPUT_Oligonucleotides_RAW.fasta', "fasta")
    logger.info(f'Oligs were written in OUTPUT_Oligonucleotides_RAW.fasta file!')

    if params.postprocess_oligos:
        oligos = rop.postprocess_oligos(tmpdir, outdir, params.vsearch_binary, params.target_max_id, logger)

    # Write oligonucleotides to excel and bed files
    oligos_df = pd.DataFrame([o.split(' ') for o in list(map(str,oligos))], 
                             columns=[
                                'Region_Name', 'GC', 'Entropy', 'Tm', 'Masked', 
                                '%1-nt_repeats', '%2-nt_repeats', '%3-nt_repeats', 
                                'Off-target_100%', 'Intersection', 'Length', 'Qual', 
                                'Seq', 'Chromosome', 'Start', 'End'
                            ])
    numeric_cols = [
        'GC', 'Entropy', 'Tm', 'Masked',
        '%1-nt_repeats', '%2-nt_repeats', '%3-nt_repeats',
        'Off-target_100%', 'Intersection', 'Length', 'Qual',
        'Start', 'End'
    ]
    oligos_df[numeric_cols] = oligos_df[numeric_cols].apply(pd.to_numeric)
    
    # Perform BLAST analysis and predict Off-target Nbases
    if params.predict_blast_offtarget:
        logger.info("Predicting BLAST Off-target Nbases...")
        models.init_models()
        oligos_df[['BLAST_offtarget_hits', 'BLAST_offt_Nbases']] = oligos_df['Seq'].apply(
            models.get_blast_offtarget, args=(params.blast_binary, params.blast_db, params.threads)
        ).apply(pd.Series)
            
    oligos_df.insert(0, 'Olig_id', [f'olig_{i+1}' for i in range(oligos_df.shape[0])])
    oligos_df.to_excel(f'{outdir}/OUTPUT_Oligonucleotides.xlsx', index=False)
    logger.info(f'Oligonucleotides were written in {outdir}/OUTPUT_Oligonucleotides.xlsx file')

    oligos_bed = oligos_df[['Chromosome', 'Start', 'End', 'Olig_id', 'Region_Name']]
    oligos_bed.to_csv(f'{outdir}/OUTPUT_Oligonucleotides.bed', sep='\t', index=False, header=False)
    
    # Perform analysis of designed oligonucleotides
    if params.analytics:
        models.init_models()
        fmindex = DNAFMIndex.read_index_from_file(params.genome_index, False)
        regions_qual_dict = rod.calculate_regions_qual(input_fasta_list, fmindex, logger)
        logger.info("Panel analysis results:")
        olan.analyze_oligos_draw_plots(oligos_df, regions_qual_dict, tmpdir, logger, outdir)

    # The end of the program!
    timestamp = datetime.now()
    logger.info(f'######### Finish OligoDesigner at {timestamp.day}.{timestamp.month}.{timestamp.year} {timestamp.hour}:{timestamp.minute}:{timestamp.second} #########')

if __name__ == '__main__':
    import sys
    _DEV = '--dev' in sys.argv
    # python3 oligodesigner.py --dev

    if _DEV:
        sys.argv.remove('--dev')
        _defaults = dict(
            default_input_genes='../tests/FileTest_input_genes_ATLAS',
            default_input_muts='',
            default_input_bed='',
            default_input_fasta='',
            default_output='../tests/6genes_results',
            default_tmpdir='../tests/6genes_results',
            default_annot='/home/a_makarova/Data/gencode.v49.annotation.gtf',
            default_genome='/home/a_makarova/Data/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna',
            default_genome_index='/home/a_makarova/Data/GRCh38.primary_assembly.genome.fa.awfmi',
        )
    else:
        _defaults = {}   # all path args default to None, supply them via CLI

    _original_parse_args = parse_args
    parse_args = lambda: _original_parse_args(**_defaults)

    main()
