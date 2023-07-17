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

import regions_and_sequences_tools as rstools
print(rstools.__dir__())
def create_logger(timestamp):
    """
    This is the docstring!
    """
    logging.basicConfig(level=logging.INFO)

    fh = logging.FileHandler(f'log_oligodesigner_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}.txt', mode='w') #-{timestamp.minute}
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))

    logger = logging.getLogger()
    logger.addHandler(fh)

    return logger

def parse_args():
    """
    This is the docstring!
    """
    args_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    args_parser.add_argument('-i','--input_genes',help='path to input file with genes',type=str,required=False, default='test')
    args_parser.add_argument('-m','--input_muts', help='path to input file with mutations',type=str,required=False, default='test_muts')
    args_parser.add_argument('-b','--input_bed', help='path to input file with tartget regions in bed format',type=str,required=False, default=False)
    args_parser.add_argument('-f','--input_fasta', help='path to input file with tartget sequences in fasta format',type=str,required=False, default=False)
    args_parser.add_argument('-o','--output', help='path of output folder',type=str,required=False, default="./testdir")
    args_parser.add_argument('-t','--tmpdir', help='path to directory with temporary files',type=str,required=False,default='./tmpdir')
    args_parser.add_argument('-a','--annot', help='path to file with genome annotation in GTF format',required=False,default='/home/a_vikhorev/Data/Genome/Homo_sapiens.GRCh38.109.gtf')
    args_parser.add_argument('-g','--genome', help='path to file with reference genome in fasta format',required=False,default='/home/a_vikhorev/Data/Genome/Homo_sapiens.GRCh38.dna.toplevel.fa')
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
        os.mkdir(dir_path)
        logger.info(f'The directory {dir_name}: {dir_path} was created!')

def check_file(file, logger, file_name):
    if os.path.exists(file):
        logger.info(f'The {file_name} {file}')
    else:
        logger.error(f'ERROR! The {file_name} {file} DOES NOT EXIST!')
        raise FileNotFoundError
        exit(1)

def main():
    # Setup logging
    timestamp = datetime.now()
    logger = create_logger(timestamp)
    logger.info(f'######### Start Pipeline at {timestamp.day}.{timestamp.month}.{timestamp.year} {timestamp.hour}:{timestamp.minute}:{timestamp.second} #########')

    # Reads params
    logger.info('The OligoDesigner was launched with the following parameters:')
    params = parse_args()

    input_genes_list = parse_input(params.input_genes, logger, 'genes') if params.input_genes else []
    input_muts_list = parse_input(params.input_muts, logger, 'mutations') if params.input_muts else []
    input_bed_list = parse_input(params.input_bed, logger, 'coordinates in bed foramt') if params.input_bed else []
    input_fasta_list = parse_input(params.input_fasta, logger, 'sequences in fasta') if params.input_fasta else []

    check_dir(params.output, logger, 'for output files')
    check_dir(params.tmpdir, logger, 'for temporary files')
    check_file(params.annot, logger, 'GTF file with annotation')
    check_file(params.genome, logger, 'FASTA file with reference genome')

    # Get bed regions for input genes and/or mutations
    if input_genes_list:
        input_genes_bed_list = rstools.get_regions_for_genes(input_genes_list, params.annot)
        input_bed_list.extend(input_genes_bed_list)
    if input_muts_list:
        input_muts_bed_list = rstools.get_regions_for_muts(input_muts_list, params.annot)
        input_bed_list.extend(input_muts_bed_list)

    if input_bed_list:
        input_fasta_intermed_list = rstools.get_fasta_seqs_from_bed(input_bed_list, params.genome)
        input_fasta_list.extend(input_fasta_intermed_list)
    
    logger.info(f'######### Finish Pipeline at {timestamp.day}.{timestamp.month}.{timestamp.year} {timestamp.hour}:{timestamp.minute}:{timestamp.second} #########')

if __name__ == '__main__':
    main()