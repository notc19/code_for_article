from ftplib import FTP
import re
import os
import shutil
import gzip

def ftp_download(ftp_link, directory, file_name, dir_output):
    '''
    Download file from ftp server.
    :param ftp_link: Link to ftp server, e.g. 'ftp.ncbi.nlm.nih.gov'.
    :type ftp_link: string
    :param directory: Directory on ftp server, where the file is located, e.g. 'refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers'.
    :type directory: string
    :param file_name: Name of file that should be downloaded from ftp server, e.g. 'GCF_000001405.39_GRCh38.p13_genomic.fna'.
    :type file_name: string
    :param dir_output: Path to directory for downloaded files.
    :type dir_output: string
    :return: Path to downloaded file.
    :rtype: string
    '''
    
    ftp = FTP(ftp_link)
    ftp.login()  # login to ftp server
    ftp.cwd(directory)  # move to directory

    files = ftp.nlst()

    for file in files:
        if re.match(file_name, file):
            file_output = os.path.join(dir_output, file)
            ftp.retrbinary('RETR ' + file, open(file_output, 'wb').write)

    ftp.quit()

    return file_output

def get_ensemble_ftp_parameters(species, annotation_release, genome_assembly, dir_output, logger):
    '''
    Get ensembl parameters
    :param species: _description_
    :type species: _type_
    :param annotation_release: _description_
    :type annotation_release: _type_
    :param genome_assembly: _description_
    :type genome_assembly: _type_
    :param dir_output: _description_
    :type dir_output: _type_
    :return: _description_
    :rtype: _type_
    '''
    ftp_link = 'ftp.ensembl.org'

    if species == 'human':
        species_id = 'homo_sapiens'
    elif species == 'mouse':
        species_id = 'mus_musculus'
    else:
        species_id = species

    if annotation_release == 'current':
        try:
            file_readme = ftp_download(ftp_link, 'pub/', 'current_README', dir_output)
            with open(file_readme, 'r') as handle:
                for line in handle:
                    if line.startswith('Ensembl Release'):
                        annotation_release = line.strip().split(' ')[2]
            os.remove(file_readme)
        except:
            logger.error(f'Cannot retrieve annotation release! You passed annotation_release={annotation_release}. Check the correctness! Exiting')
            exit(1)

    ftp_directory_gtf = f'pub/release-{annotation_release}/gtf/{species_id}/'
    ftp_directory_fasta = f'pub/release-{annotation_release}/fasta/{species_id}/dna/'
    ftp_file_gtf = f'{species_id.capitalize()}.{genome_assembly}.{annotation_release}.gtf.gz'
    ftp_file_fasta = f'{species_id.capitalize()}.{genome_assembly}.dna.toplevel.fa.gz'

    return ftp_link, ftp_directory_gtf, ftp_directory_fasta, ftp_file_gtf, ftp_file_fasta

def decompress_gzip(file_gzip):
    '''
    Decompress zip files.
    :param file_gzip: Path to zipped file.
    :type file_gzip: string
    :return: Path to unzipped file.
    :rtype: string
    '''
    file_output = file_gzip.split('.gz')[0]
    with gzip.open(file_gzip, 'rb') as f_in:
        with open(file_output, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_gzip)
    return file_output

def download_from_ensembl(typ, species, annotation_release, genome_assembly, dir_output, logger):
    
    ftp_link, ftp_directory_gtf, ftp_directory_fasta, ftp_file_gtf, ftp_file_fasta = \
        get_ensemble_ftp_parameters(species, annotation_release, genome_assembly, dir_output, logger)

    if typ == 'annotation':
        try:
            file = ftp_download(ftp_link, ftp_directory_gtf, ftp_file_gtf, dir_output)
        except:
            logger.error(f'Cannot downlowd {ftp_link}/{ftp_directory_gtf}/{ftp_file_gtf}! Check the correctness of the link')
            exit(1)
    elif typ == 'genome':
        try:
            file = ftp_download(ftp_link, ftp_directory_fasta, ftp_file_fasta, dir_output)
        except:
            logger.error(f'Cannot downlowd {ftp_link}/{ftp_directory_gtf}/{ftp_file_gtf}! Check the correctness of the link')
            exit(1)

    file = decompress_gzip(file)

    return file

if __name__ == '__main__':
    typ = 'genome'
    species = 'homo_sapiens'
    annotation_release = 'current'
    genome_assembly = 'GRCh38'
    dir_output = './'

    annot = download_from_ensembl(typ, species, annotation_release, genome_assembly, dir_output)