from ftplib import FTP
from pathlib import Path
import logging
import re
import os
import shutil
import gzip

# ---------------------------------------------------------------------------
# Ensembl
# ---------------------------------------------------------------------------

def ftp_download(ftp_link: str, directory: str, file_name: str, dir_output: str | Path) -> str:
    """
    Скачивает файл с ftp сервера.
    :param ftp_link:    Ссылка на ftp сервер, например 'ftp.ncbi.nlm.nih.gov'.
    :param directory:   Directory on ftp server, where the file is located, e.g. 'refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers'.
    :param file_name:   Name of file that should be downloaded from ftp server, e.g. 'GCF_000001405.40_GRCh38.p14_genomic.fna.gz'
    :param dir_output:  Path to directory for downloaded files.
    :return: Path to downloaded file.
    """
    dir_output = str(dir_output)

    ftp = FTP(ftp_link)
    ftp.login()
    ftp.cwd(directory)

    files = ftp.nlst()

    file_output = None
    for file in files:
        if re.match(file_name, file):
            file_output = os.path.join(dir_output, file)
            with open(file_output, 'wb') as fh:
                ftp.retrbinary('RETR ' + file, fh.write)

    ftp.quit()

    if file_output is None:
        raise FileNotFoundError(f"No file matching '{file_name}' found in ftp://{ftp_link}/{directory}")

    return file_output


def get_ensemble_ftp_parameters(species: str, annotation_release: str,
                                genome_assembly: str, dir_output: str | Path,
                                logger: logging.Logger) -> tuple:
    """
    Эта функция создает параметры, которые необходимы для скачивания файла из Ensembl
    :param species:             название вида (например, homo_sapiens)
    :param annotation_release:  релиз аннотации (например, release-110)
    :param genome_assembly:     версия сборки генома (например, GRCh38)
    :param dir_output:          директория для скачанного файла (скачивается и удаляется readme с Ensembl для поиска current версии аннотации)
    :param logger:              объект класса Logger
    :return: строки с путем для скачивания файлов с ftp Ensembl, i.e. 
            (ftp_link, ftp_directory_gtf, ftp_directory_fasta, ftp_file_gtf, ftp_file_fasta)
    """
    ftp_link = 'ftp.ensembl.org'

    if species == 'human':
        species_id = 'homo_sapiens'
    elif species == 'mouse':
        species_id = 'mus_musculus'
    else:
        species_id = species

    if annotation_release == 'current':
        try:
            file_readme = ftp_download(ftp_link, 'pub/', 'current_README',
                                       dir_output)
            with open(file_readme, 'r') as handle:
                for line in handle:
                    if line.startswith('Ensembl Release'):
                        try:
                            annotation_release = line.strip().split(' ')[2]
                        except IndexError as e:
                            logger.error(
                                f'Error {e}. Unexpected README format; could not parse release number from line: {line.strip()!r}'
                            )
                            raise
            os.remove(file_readme)
        except Exception as e:
            logger.error(
                f'Cannot retrieve annotation_release={annotation_release!r}. Check the correctness! Error: {e}'
            )
            raise

    species_id_title = '_'.join(w.capitalize() for w in species_id.split('_'))

    ftp_directory_gtf   = f'pub/release-{annotation_release}/gtf/{species_id}/'
    ftp_directory_fasta = f'pub/release-{annotation_release}/fasta/{species_id}/dna/'
    ftp_file_gtf   = f'{species_id_title}.{genome_assembly}.{annotation_release}.gtf.gz'
    ftp_file_fasta = f'{species_id_title}.{genome_assembly}.dna.toplevel.fa.gz'

    return ftp_link, ftp_directory_gtf, ftp_directory_fasta, ftp_file_gtf, ftp_file_fasta


def decompress_gzip(file_gzip: str) -> str:
    """
    Распаковывает .gz архив
    :param file_gzip: Путь к архивированному файлу
    :return:          Путь к разархивированному файлу
    """
    file_output = file_gzip.split('.gz')[0]
    with gzip.open(file_gzip, 'rb') as f_in:
        with open(file_output, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_gzip)
    return file_output


def download_from_ensembl(typ: str, species: str, annotation_release: str,
                          genome_assembly: str, dir_output: str | Path,
                          logger: logging.Logger) -> Path:
    """
    Эта функция выполняет скачивание из Ensembl генома или аннотации, используя на входе 
    название вида, релиз Ensembl и версию сборки генома.
    :param typ:                'annotation' (GTF) or 'genome' (FASTA)
    :param species:            'human', 'mouse', or a raw Ensembl species ID
    :param annotation_release: Ensembl release number or 'current'
    :param genome_assembly:    Genome assembly version (e.g. 'GRCh38')
    :param dir_output:         Объект Path с путем до директории, в которую будет сохранен файл
    :param logger:             объект класса Logger
    :return:                   Объект Path с путем к скачанному файлу
    """
    ftp_link, ftp_directory_gtf, ftp_directory_fasta, ftp_file_gtf, ftp_file_fasta = \
        get_ensemble_ftp_parameters(species, annotation_release, genome_assembly, dir_output, logger)

    if typ == 'annotation':
        try:
            file = ftp_download(ftp_link, ftp_directory_gtf,
                                ftp_file_gtf, dir_output)
        except Exception as e:
            logger.error(
                f'Cannot download ftp://{ftp_link}/{ftp_directory_gtf}/{ftp_file_gtf}! Check the correctness of the link.'
            )
            logger.error(e)
            exit(1)

    elif typ == 'genome':
        try:
            file = ftp_download(ftp_link, ftp_directory_fasta,
                                ftp_file_fasta, dir_output)
        except Exception as e:
            # FIX: was incorrectly logging the GTF path for genome failures.
            logger.error(
                f'Cannot download ftp://{ftp_link}/{ftp_directory_fasta}/{ftp_file_fasta}! Check the correctness of the link.'
            )
            logger.error(e)
            exit(1)

    else:
        logger.error(f"Unknown typ={typ!r}. Expected 'annotation' or 'genome'.")
        exit(1)

    file = decompress_gzip(file)
    return Path(file)


# ---------------------------------------------------------------------------
# GENCODE
# ---------------------------------------------------------------------------

_GENCODE_SPECIES_META = {
    'human': {
        'ftp_subdir':       'Gencode_human',
        'genome_assembly':  'GRCh38',
        'release_prefix':   'v',        # e.g. v49
        'gtf_re':           re.compile(r'gencode\.(v\d+)\.primary_assembly\.annotation\.gtf\.gz'),
    },
    'mouse': {
        'ftp_subdir':       'Gencode_mouse',
        'genome_assembly':  'GRCm39',
        'release_prefix':   'vM',       # e.g. vM38
        'gtf_re':           re.compile(r'gencode\.(vM\d+)\.primary_assembly\.annotation\.gtf\.gz'),
    },
}


def get_gencode_ftp_parameters(species: str, annotation_release: str,
                               logger: logging.Logger) -> tuple:
    """
    Build FTP parameters required to download files from the GENCODE FTP.

    Unlike Ensembl, GENCODE stores both the GTF and the genome FASTA in the
    same release directory, so a single ftp_directory is returned.

    Supported species: 'human' (GRCh38, releases v1–v49+) and
                       'mouse' (GRCm39, releases vM1–vM38+).

    Pass annotation_release='current' to automatically resolve to the latest
    stable release via the GENCODE latest_release/ FTP directory.  The actual
    release number is parsed from the GTF filename present in that directory
    (e.g. 'gencode.v49.primary_assembly.annotation.gtf.gz' → '49') — no
    README download is required.

    Pass an explicit release number string (e.g. '49' for human, 'M38' for mouse) to pin to a specific release.

    :param species:            'human' or 'mouse'.
    :param annotation_release: Release number string or 'current'.
    :param logger:             Logger instance.
    :return: (ftp_link, ftp_directory, ftp_file_gtf, ftp_file_fasta, resolved_release)
             resolved_release is the numeric string actually used (useful when 'current' was passed).
    """
    if species not in _GENCODE_SPECIES_META:
        logger.error(
            f"Unsupported species={species!r} for GENCODE. Supported: {list(_GENCODE_SPECIES_META)}"
        )
        exit(1)

    meta             = _GENCODE_SPECIES_META[species]
    ftp_link         = 'ftp.ebi.ac.uk'
    ftp_subdir       = meta['ftp_subdir']
    genome_assembly  = meta['genome_assembly']
    release_prefix   = meta['release_prefix']
    gtf_re           = meta['gtf_re']

    if annotation_release == 'current':
        ftp_directory = f'pub/databases/gencode/{ftp_subdir}/latest_release/'
        try:
            ftp = FTP(ftp_link)
            ftp.login()
            ftp.cwd(ftp_directory)
            files = ftp.nlst()
            ftp.quit()
        except Exception as e:
            logger.error(
                f'Cannot list GENCODE FTP directory ftp://{ftp_link}/{ftp_directory}. Error: {e}'
            )
            exit(1)

        resolved_release = None
        for fname in files:
            m = gtf_re.match(fname)
            if m:
                resolved_release = m.group(1).lstrip('vM').lstrip('v')
                break
        if resolved_release is None:
            logger.error(
                f'Could not determine current GENCODE release for {species!r} from directory listing of '
                f'ftp://{ftp_link}/{ftp_directory}'
            )
            exit(1)
        logger.info(
            f"Resolved GENCODE 'current' release for {species!r} to: {release_prefix}{resolved_release}"
        )
    else:
        resolved_release = annotation_release.lstrip('vM').lstrip('v')
        ftp_directory = (
            f'pub/databases/gencode/{ftp_subdir}/release_{resolved_release}/'
        )

    version_tag = f'{release_prefix}{resolved_release}'
    ftp_file_gtf = f'gencode.{version_tag}.primary_assembly.annotation.gtf.gz'
    ftp_file_fasta = f'{genome_assembly}.primary_assembly.genome.fa.gz'

    return ftp_link, ftp_directory, ftp_file_gtf, ftp_file_fasta, resolved_release


def download_from_gencode(typ: str, species: str, annotation_release: str, dir_output: str | Path,
                          logger: logging.Logger) -> Path:
    """
    Download a genome or annotation file from the GENCODE FTP.
    Mirrors the API of download_from_ensembl().  The genome_assembly argument
    is not needed here because GENCODE encodes the assembly in the filename
    and it is determined automatically from the species.

    :param typ:                'annotation' (GTF) or 'genome' (FASTA).
    :param species:            'human' or 'mouse'.
    :param annotation_release: GENCODE release number (e.g. '49' for human, 'M38' for mouse) or 'current'.
    :param dir_output:         Directory where the decompressed file is saved.
    :param logger:             Logger object.
    :return:                   Path to the decompressed file.
    """
    ftp_link, ftp_directory, ftp_file_gtf, ftp_file_fasta, resolved_release = \
        get_gencode_ftp_parameters(species, annotation_release, logger)

    if typ == 'annotation':
        target_file = ftp_file_gtf
    elif typ == 'genome':
        target_file = ftp_file_fasta
    else:
        logger.error(f"Unknown typ={typ!r}. Expected 'annotation' or 'genome'.")
        exit(1)

    try:
        file = ftp_download(ftp_link, ftp_directory, target_file, dir_output)
    except Exception as e:
        logger.error(
            f'Cannot download ftp://{ftp_link}/{ftp_directory}/{target_file}! '
            f'Check the correctness of the link.'
        )
        logger.error(e)
        exit(1)

    file = decompress_gzip(file)
    return Path(file)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    logger = logging.getLogger(__name__)

    output_path = download_from_ensembl(
        typ='genome',
        species='homo_sapiens',
        annotation_release='current',
        genome_assembly='GRCh38',
        dir_output='./',
        logger=logger,
    )
    logger.info(f'Ensembl genome saved to: {output_path}')

    gencode_path = download_from_gencode(
        typ='annotation',
        species='human',
        annotation_release='current',
        dir_output='./',
        logger=logger,
    )
    logger.info(f'GENCODE annotation saved to: {gencode_path}')
