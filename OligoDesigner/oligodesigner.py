###########################################
# imports
###########################################

import os
import sys
import logging
import time
from timeit import default_timer as timer
from datetime import datetime
import argparse

import get_regions_for_mutations as grfm
import get_regions_for_genes as grfg
import get_fasta_seq_from_bed as gfsfb
import utils

def create_logger(timestamp):
    """
    This is the docstring!
    """
    logging.basicConfig(level=logging.INFO)

    fh = logging.FileHandler(f'log_oligodesigner_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}.txt', mode='w') #-{timestamp.minute}
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    logger = logging.getLogger()
    logger.addHandler(fh)

    return logger

def parse_args():
    """
    This is the docstring!
    """
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('-i','--input_genes',help='path to input file with genes',type=str,required=False, default='tests/FileTest_input_genes')
    args_parser.add_argument('-m','--input_muts', help='path to input file with mutations',type=str,required=False, default='tests/FileTest_input_muts_hgvs') #tests/FileTest_input_muts_hgvs
    args_parser.add_argument('-b','--input_bed', help='path to input file with tartget regions in bed format',type=str,required=False, default=False)
    args_parser.add_argument('-f','--input_fasta', help='path to input file with tartget sequences in fasta format',type=str,required=False, default=False)
    args_parser.add_argument('-o','--output', help='path of output folder',type=str,required=False, default="./tests/test_reuslts")
    args_parser.add_argument('-t','--tmpdir', help='path to directory with temporary files',type=str,required=False,default='./tests/test_reuslts')
    args_parser.add_argument('-a','--annot', help='path to file with genome annotation in GTF format',required=False,default='/home/a_vikhorev/Data/Genome/Homo_sapiens.GRCh38.109.gtf')
    args_parser.add_argument('-g','--genome', help='path to file with reference genome in fasta format',required=False,default='/home/a_vikhorev/Data/DataBases/Hg38_GCA_000001405.15_BWA/GCA_000001405.15.fasta')
    args_parser.add_argument('-l','--bait_length', help='Length of the designed olignucleotide bait', required=False, default=100)
    args_parser.add_argument('-k','--mutation_flank', help='Flank length from mutation', required=False, default=0)
    args_parser.add_argument('-w','--write_intermediate', help='write intermediate results or not', required=False, default=True)
    args_parser.add_argument('-v','--genome_version', help='version of the reference genome', required=False, default='GRCh38.p14')
    args_parser.add_argument('-r','--gene_regions', help='Type of gene regions to use [full/exons/exons_utr]', required=False, default='exon')
    args_parser.add_argument('-j','--intron_depth', help='Depth of introns to include when only exons regions were chosen', required=False, default=50)
    
    return args_parser.parse_args()

def parse_input(inpt, logger, feature):
    try:
        with open(inpt) as infile:
            logger.info(f'The File with input {feature}: {inpt}')
            input = [line.strip() for line in infile.readlines()]
            logger.info(f'{len(input)} entries were found in {inpt} file')

            assert len(input) > 0

            return input

    except FileNotFoundError as e:
        logger.error(f'File {inpt} does not exist!!!!')
        return False

def check_dir(dir_path, logger, dir_name):
    if os.path.exists(dir_path):
        logger.info(f'The directory {dir_name}: {dir_path} exists')
    else:
        try:
            os.mkdir(dir_path)
            logger.info(f'The directory {dir_name}: {dir_path} was created!')
        except:
            logger.error(f'Cannot create {dir_path} directory!')
            exit(1)

def check_file(file, logger, file_name):
    if os.path.exists(file):
        logger.info(f'The {file_name} {file} was found!')
    else:
        logger.error(f'ERROR! The {file_name} {file} DOES NOT EXIST!')
        raise FileNotFoundError
        exit(1)

def main():
    # Setup logging
    timestamp = datetime.now()
    logger = create_logger(timestamp)
    logger.info(f'######### Start OligoDesigner at {timestamp.day}.{timestamp.month}.{timestamp.year} {timestamp.hour}:{timestamp.minute}:{timestamp.second} #########')

    # Reads params
    logger.info('The OligoDesigner was launched with the following parameters:')
    params = parse_args()

    input_genes_list = parse_input(params.input_genes, logger, 'genes') if params.input_genes else []
    input_muts_list = parse_input(params.input_muts, logger, 'mutations') if params.input_muts else []
    input_bed_list = parse_input(params.input_bed, logger, 'coordinates in bed foramt') if params.input_bed else []
    input_fasta_list = parse_input(params.input_fasta, logger, 'sequences in fasta') if params.input_fasta else []

    check_dir(params.output, logger, 'for output files')
    check_dir(params.tmpdir, logger, 'for temporary files')

    if not input_bed_list and not input_fasta_list:
        check_file(params.annot, logger, 'GTF file with annotation')

    if input_genes_list or input_muts_list:
        with open(params.annot) as infile:
            annot = [line.strip() for line in infile.readlines() if not line.startswith('#')]
            if input_genes_list:
                annot_for_genes, gene_name_dict = utils.prepare_annot_for_genes(annot, logger)
            if input_muts_list:
                annot_for_muts = utils.prepare_annot_for_muts(input_muts_list, annot, logger)

    if not input_fasta_list:
        check_file(params.genome, logger, 'FASTA file with reference genome')

    # Get bed regions for input genes and/or mutations
    if input_genes_list:
        input_genes_bed_list = grfg.get_regions_for_genes(input_genes_list, annot_for_genes, gene_name_dict,
                                                             params.gene_regions, params.intron_depth,
                                                             params.write_intermediate, params.tmpdir, logger,
                                                             params.genome_version)
        input_bed_list.extend(input_genes_bed_list)

    if input_muts_list:
        flank_length = int(params.bait_length * 1.7 / 2) if not params.mutation_flank else params.mutation_flank
        input_muts_bed_list = grfm.get_regions_for_muts(input_muts_list, annot_for_muts, flank_length, \
                                                           params.write_intermediate, params.tmpdir, logger, \
                                                           params.genome_version)
        input_bed_list.extend(input_muts_bed_list)

    if input_bed_list:
        input_fasta_intermed_list = gfsfb.get_fasta_seqs_from_bed(input_bed_list, params.genome, params.genome_version,
                                                                    params.write_intermediate, params.tmpdir, logger)
        input_fasta_list.extend(input_fasta_intermed_list)
    
    timestamp = datetime.now()
    logger.info(f'######### Finish OligoDesigner at {timestamp.day}.{timestamp.month}.{timestamp.year} {timestamp.hour}:{timestamp.minute}:{timestamp.second} #########')

if __name__ == '__main__':
    main()