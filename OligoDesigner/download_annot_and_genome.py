from ftplib import FTP
from pathlib import Path
import logging
import re
import os
import shutil
import gzip

def ftp_download(ftp_link:str, directory:str, file_name:str, dir_output:str):
    '''
    Скачивает файл с ftp сервера.
    :param ftp_link: Ссылка на ftp сервер, например 'ftp.ncbi.nlm.nih.gov'.
    :param directory: Directory on ftp server, where the file is located, e.g. 'refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers'.
    :param file_name: Name of file that should be downloaded from ftp server, e.g. 'GCF_000001405.39_GRCh38.p13_genomic.fna'.
    :param dir_output: Path to directory for downloaded files.
    :return: Path to downloaded file.
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

def get_ensemble_ftp_parameters(species:str, annotation_release:str, 
                                genome_assembly:str, dir_output:str, logger:logging):
    '''
    Эта функция создает параметры, которые необходимы для скачивания файла из Ensembl
    :param species: название вида (например, homo_sapiens)
    :param annotation_release: релиз аннотации (например, release-110)
    :param genome_assembly: версия сборки генома (например, GRCh38)
    :param dir_output: директория для скачанного файла (скачивается и удаляется readme с Ensembl 
                       для поиска current версии аннотации)
    :return: строки с путем для скачивания файлов с ftp Ensembl
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

def decompress_gzip(file_gzip:str) -> str:
    '''
    Распаковывает gz архив
    :param file_gzip: Путь к архивированному файлу
    :return: Путь к разархивированному файлу
    '''
    file_output = file_gzip.split('.gz')[0]
    with gzip.open(file_gzip, 'rb') as f_in:
        with open(file_output, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_gzip)
    return file_output

def download_from_ensembl(typ:str, species:str, annotation_release:str, 
                          genome_assembly:str, dir_output:Path, logger:logging) -> Path:
    """
    Эта функция выполняет скачивание из Ensembl генома или аннотации, используя на входе 
    название вида, релиз Ensembl и версию сборки генома.
    :param typ: Тип файла, который нужно скачать (геном/аннотация)
    :param species: Вид, для которого нужно скачать геном/аннотацию
    :param annotation_release: Релиз аннотации
    :param genome_assembly: версия сборки генома
    :param dir_output: Объект Path с путем до директории, в которую будет сохранен файл
    :param logger: объект класса logging
    :return: Объект Path с путем к скачанному файлу
    """
    ftp_link, ftp_directory_gtf, ftp_directory_fasta, ftp_file_gtf, ftp_file_fasta = \
        get_ensemble_ftp_parameters(species, annotation_release, genome_assembly, dir_output, logger)

    if typ == 'annotation':
        try:
            file = ftp_download(ftp_link, ftp_directory_gtf, ftp_file_gtf, dir_output)
        except Exception as e:
            logger.error(f'Cannot downlowd {ftp_link}/{ftp_directory_gtf}/{ftp_file_gtf}! Check the correctness of the link')
            logger.error(e)
            exit(1)
    elif typ == 'genome':
        try:
            file = ftp_download(ftp_link, ftp_directory_fasta, ftp_file_fasta, dir_output)
        except Exception as e:
            logger.error(f'Cannot downlowd {ftp_link}/{ftp_directory_gtf}/{ftp_file_gtf}! Check the correctness of the link')
            logger.error(e)
            exit(1)

    file = decompress_gzip(file)

    return Path(file)

if __name__ == '__main__':
    typ = 'genome'
    species = 'homo_sapiens'
    annotation_release = 'current'
    genome_assembly = 'GRCh38'
    dir_output = './'

    annot = download_from_ensembl(typ, species, annotation_release, genome_assembly, dir_output)