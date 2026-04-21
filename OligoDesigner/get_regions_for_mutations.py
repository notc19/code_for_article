from logging import Logger
from pathlib import Path
import re
from constants import genome_versions_chromosome_dict, rev_dict

def _strip_version(entry: str) -> str:
    """
    Strip the version suffix from a versioned Ensembl ID/gene name.
    "ENST00000302118.5"  →  "ENST00000302118"
    "BRCA1.3"            →  "BRCA1"
    "ENST00000302118"    →  "ENST00000302118"
    :param entry: versioned or unversioned entry
    :return: unversioned entry
    """
    return entry.split(".")[0]

def build_unversioned_muts_index(annot_for_muts: dict) -> dict[str, str]:
    """
    Build a lookup index {unversioned_transcript_id: versioned_transcript_id}
    from the dict returned by utils.prepare_annot_for_muts.
    :param annot_for_muts: dict {versioned_tid: [chrom, exons, start_codon, strand]}
    :return: dict {unversioned_tid: versioned_tid}
    """
    return {_strip_version(tid): tid for tid in annot_for_muts}


def get_gene_id_for_transcript(annot_for_genes: dict) -> dict[str, str]:
    """
    Build index {unversioned_transcript_id: gene_id} from annot (utils.prepare_annot_for_genes).
    Resolves a transcript ID (from an unknown annotation version) to its gene,
    so the MANE_Select/basic transcript for that gene can be used.

    annot_for_genes values:
      [0] chrom  [1] gene_coords  [2] exon_list  [3] utr_list  [4] exon_id_list
      [5] mane_select  [6] mane_plus_clinical  [7] transcript_list  [8] selected_tid
    :param annot_for_genes: dict {gene_id: [<annot_for_genes values>]}
    :return: dict {unversioned_tid: gene_id}
    """
    tid_to_gene: dict[str, str] = {}
    for gene_id, info in annot_for_genes.items():
        for tid in info[7]:
            tid_to_gene[_strip_version(tid)] = gene_id
    return tid_to_gene


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
    i = re.match(r'\d+', pos_and_nucl).end()
    pos = int(pos_and_nucl[:i])
    ref = pos_and_nucl[i:].split('>')[0]

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref + '\t' + ':'.join(hgvsg) + '\t' + 'snp'

def get_pos_c_pos_m_ref(pos_and_nucl:str, sign:str) -> list[int, int, str]:
    """
    Эта функция извлекает позицию, смещение позиции и референсный нуклеотид из мутации в формате hgvsc (hgvsc.split()[1])
    :param pos_and_nucl: часть мутации в формате hgvsc, содержащая позиции и референсный нуклеотид
    :param sign: знак, определяющий позицию в интроне относительно экзона (+ или -)
    :return: Список с координатой в экзоне (pos_c), смещением координаты в интрон (pos_m), референсным нуклеотидом (ref)
    """
    pos_c = int(pos_and_nucl.split(sign)[0])
    j = re.match(r'\d+', pos_and_nucl.split(sign)[1]).end()
    ref = pos_and_nucl.split(sign)[1][j:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
        else 'del'+pos_and_nucl.split('del')[1] if 'ins' not in pos_and_nucl \
        else 'ins'+pos_and_nucl.split('ins')[1]
    pos_m = int(pos_and_nucl.split(sign)[1][:j])
    if '_' in str(pos_and_nucl):
        pos_c = int(str(pos_and_nucl).split('_')[0].split(sign)[0])
        pos_m = int(str(pos_and_nucl).split('_')[0].split(sign)[1][:j])
    
    return pos_c, pos_m, ref

def get_pos(strand:str, start_codon_coord:int, exons:list, pos_c:int, ref:str) -> tuple[int, str]:
    """
    Эта функция находит координату в геноме для мутации в формате hgvsc
    :param strand: цепь, на которой находится мутация
    :param start_codon_coord: координата старт-кодона
    :param exons: список с экзонами гена
    :param pos_c: позиция в экзоне
    :param ref: референсный нуклеотид
    :return: Возвращает позицию в геноме и референсный нуклеотид
    """
    exons = [list(e) for e in exons]
    i = 0
    if strand == '+':
        while True:
            #added bounds check to prevent infinite loop / IndexError when pos_c
            if i >= len(exons):
                raise ValueError(
                    f'pos_c exceeds total coding length of transcript (remaining pos_c={pos_c}, exon index={i})'
                )
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
            if i >= len(exons):
                raise ValueError(
                    f'pos_c exceeds total coding length of transcript (remaining pos_c={pos_c}, exon index={i})'
                )
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

def _resolve_transcript(
    input_tid: str,
    annot_for_muts: dict,
    unversioned_index: dict[str, str],
    tid_to_gene: dict[str, str],
    annot_for_genes: dict,
    logger: Logger,
) -> tuple[str | None, str | None, str | None]:
    """
    Resolve an input transcript ID (versioned or unversioned, from any annotation
    version) to a transcript entry present in annot_for_muts.

    1. Exact match in annot_for_muts (versioned key).
    2. Unversioned match via unversioned_index (strips ".N" from input_tid).
    3. MANE_Select fallback: look up the gene via tid_to_gene, then use the
       MANE_Select transcript of that gene (from annot_for_genes slot [5]).
       If MANE_Select is absent, use the gene's selected_transcript_id (slot [8]).

    :param input_tid: transcript ID as it appears in the mutation string
    :param annot_for_muts: dict {versioned_tid: [chrom, exons, start_codon, strand]}
    :param unversioned_index: dict {unversioned_tid: versioned_tid} (from build_unversioned_muts_index)
    :param tid_to_gene: dict {unversioned_tid: gene_id} (from get_gene_id_for_transcript)
    :param annot_for_genes: dict {gene_id: [...]} from utils.prepare_annot_for_genes
    :param logger: Logger instance
    :return: tuple (resolved_versioned_tid, gene_id, resolution_note)
    """
    # Step 1: exact match (input already versioned and present)
    if input_tid in annot_for_muts:
        gene_id = tid_to_gene.get(_strip_version(input_tid))
        return input_tid, gene_id, "exact"

    # Step 2: unversioned match (strip version from input, look up in index)
    unversioned_input = _strip_version(input_tid)
    if unversioned_input in unversioned_index:
        resolved = unversioned_index[unversioned_input]
        gene_id = tid_to_gene.get(unversioned_input)
        return resolved, gene_id, "unversioned"

    # Step 3: MANE_Select fallback via gene index
    gene_id = tid_to_gene.get(unversioned_input)
    if gene_id and gene_id in annot_for_genes:
        gene_info = annot_for_genes[gene_id]
        mane_tids = gene_info[5]                # info[5]: mane_select_transcripts
        selected_tid = gene_info[8]             # info[8]: selected_transcript_id

        # Prefer MANE_Select, fall back to the gene's selected transcript.
        candidates = mane_tids if mane_tids else ([selected_tid] if selected_tid else [])
        for candidate in candidates:
            unversioned_candidate = _strip_version(candidate)
            if candidate in annot_for_muts:
                logger.info(
                    f'{input_tid}: not found in annotation; using MANE/canonical '
                    f'transcript {candidate} for gene {gene_id}.'
                )
                return candidate, gene_id, "mane_fallback"
            if unversioned_candidate in unversioned_index:
                resolved = unversioned_index[unversioned_candidate]
                logger.info(
                    f'{input_tid}: not found in annotation; using MANE/canonical '
                    f'transcript {resolved} for gene {gene_id}.'
                )
                return resolved, gene_id, "mane_fallback"

    return None, None, None


def get_coords_for_hgvs_c(
    hgvsc: list[str],
    flank_length: int,
    annot_for_muts: dict,
    logger: Logger,
    unversioned_index: dict[str, str] | None = None,
    tid_to_gene: dict[str, str] | None = None,
    annot_for_genes: dict | None = None,
) -> str | None:
    """
    Accept a mutation in HGVS c./n. notation and return a BED-format region
    flanking the mutation.

    Transcript ID resolution
    ------------------------
    The input transcript ID (hgvsc[0]) may be versioned or unversioned and may
    come from a different annotation version than annot_for_muts.  Resolution
    proceeds in three steps (see _resolve_transcript):
      1. Exact match in annot_for_muts.
      2. Unversioned match (strip ".N" suffix).
      3. MANE_Select / canonical fallback for the same gene.

    The BED name field (column 4) always carries the *unversioned* transcript ID
    that was actually used for coordinate calculation, plus the gene ID if known,
    in the format GENE_ID | ENST_unversioned:c.notation.

    :param hgvsc: mutation split on ':', e.g. ['ENST00000357654.9', 'c.1234A>T']
    :param flank_length: flanking region length in bp
    :param annot_for_muts: dict {versioned_tid: [chrom, exons, start_codon, strand]}
    :param logger: object of class Logger
    :param unversioned_index: optional pre-built {unversioned_tid: versioned_tid};
                              built on-the-fly if None (inefficient for large inputs)
    :param tid_to_gene: optional pre-built {unversioned_tid: gene_id}
    :param annot_for_genes: optional gene annotation dict from prepare_annot_for_genes
    :return: BED-format string or None if the transcript cannot be resolved
    """
    # Build indices on-the-fly if not pre-supplied (convenience for single calls).
    if unversioned_index is None:
        unversioned_index = build_unversioned_muts_index(annot_for_muts)
    if tid_to_gene is None and annot_for_genes is not None:
        tid_to_gene = get_gene_id_for_transcript(annot_for_genes)
    if tid_to_gene is None:
        tid_to_gene = {}

    input_tid = hgvsc[0]
    trans, gene_id, resolution = _resolve_transcript(
        input_tid, annot_for_muts, unversioned_index, tid_to_gene, annot_for_genes or {}, logger
    )

    if trans is None:
        logger.warning(
            f'Transcript {input_tid} not present in provided annotation and could not be resolved via MANE fallback. Skipping {":".join(hgvsc)}.'
        )
        return None

    chrom = annot_for_muts[trans][0]
    exons = [list(e) for e in annot_for_muts[trans][1]]
    start_codon_coord = int(annot_for_muts[trans][2])
    strand = annot_for_muts[trans][3]

    # Build the canonical name for the BED name field:
    # always use unversioned transcript ID; prepend gene_id if known.
    unversioned_trans = _strip_version(trans)
    bed_name_prefix = f'{gene_id}|{unversioned_trans}' if gene_id else unversioned_trans
    # Reconstruct the mutation string with the resolved (unversioned) transcript ID.
    resolved_mut_str = f'{bed_name_prefix}:{hgvsc[1]}'

    pos_and_nucl = hgvsc[1].split('.')[1]
    if pos_and_nucl[0].isdigit():
        if not '-' in pos_and_nucl and not '+' in pos_and_nucl:
            j = re.match(r'\d+', pos_and_nucl).end()
            ref = pos_and_nucl[j:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
                else 'del'+pos_and_nucl[1:].split('del')[1] if 'ins' not in pos_and_nucl \
                else 'ins'+pos_and_nucl.split('ins')[1]
            pos_c = int(pos_and_nucl[:j]) - 3
            if '_' in str(pos_and_nucl[:j]):
                pos_c = int(str(pos_and_nucl[:j]).split('_')[0]) - 3
            try:
                pos, ref = get_pos(strand, start_codon_coord, exons, pos_c, ref)
            except ValueError as e:
                logger.warning(f'Skipping {":".join(hgvsc)}: {e}')
                return None

        elif '-' in pos_and_nucl:
            pos_c, pos_m, ref = get_pos_c_pos_m_ref(pos_and_nucl, '-')
            try:
                pos, ref = get_pos(strand, start_codon_coord, exons, pos_c, ref)
            except ValueError as e:
                logger.warning(f'Skipping {":".join(hgvsc)}: {e}')
                return None
            pos = pos - 2 - pos_m

        elif '+' in pos_and_nucl:
            pos_c, pos_m, ref = get_pos_c_pos_m_ref(pos_and_nucl, '+')
            try:
                pos, ref = get_pos(strand, start_codon_coord, exons, pos_c, ref)
            except ValueError as e:
                logger.warning(f'Skipping {":".join(hgvsc)}: {e}')
                return None
            pos = pos - 2 + pos_m

    elif pos_and_nucl.startswith('-'):
        if not '-' in pos_and_nucl[1:] and not '+' in pos_and_nucl[1:]:
            j = re.match(r'\d+', pos_and_nucl[1:]).end()
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
            j = re.match(r'\d+', pos_and_nucl[1:]).end()
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
        logger.warning(
            f'This format of HGVS mutation ({input_tid}:{hgvsc[1]}) is not supported. Skipping this variant.'
        )
        return None

    start = pos - flank_length
    end = pos + flank_length

    return f'{chrom}\t{start}\t{end}\t{ref}\t{resolved_mut_str}\tsnp'

def fix_chromosome_name_and_check_coords(coord_list:list, genome_version:str, logger:Logger) -> list:
    """
    Эта функция проверяет список регионов и фиксит при необходимости названия хромосом дял человека
    (приводит все хромосомы к названию типа chr1)
    :param coord_list: Список координат регионов, фланкирующих мутации
    :param genome_version: Версия генома (GRCh38/37)
    :param logger: Объект класса Logger
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
        except (IndexError, KeyError):
            logger.warning(f'Something wrong with {coord} variant! It has no ref allele! Skipping this variant!')
    
    return fixed_coord_list

def get_regions_for_muts(
    muts_list: list,
    annot_for_muts: dict,
    flank_length: int,
    write_intermediate: bool,
    tmpdir: Path,
    logger: Logger,
    genome_version: str,
    annot_for_genes: dict | None = None,
) -> list:
    """
    Find genomic regions flanking each mutation in muts_list.

    Transcript ID resolution for HGVS c./n. mutations
    ---------------------------------------------------
    Input transcript IDs may be versioned or unversioned and may originate from
    any annotation version.  Resolution proceeds as follows:
      1. Exact match in annot_for_muts.
      2. Unversioned match (strip ".N" suffix from input ID).
      3. MANE_Select / canonical fallback: if the transcript is known to belong
         to a gene (via annot_for_genes), the MANE_Select transcript of that gene
         is used instead.  Requires annot_for_genes to be supplied.

    The BED name field for resolved c./n. mutations is formatted as: GENE_ID|ENST_unversioned:c.notation

    :param muts_list: list of mutation strings (VCF, HGVS g., HGVS c./n. formats)
    :param annot_for_muts: dict {versioned_tid: [chrom, exons, start_codon, strand]} from utils.prepare_annot_for_muts
    :param flank_length: flanking region length in bp
    :param write_intermediate: whether to write an intermediate BED file
    :param tmpdir: Path to intermediate results directory
    :param logger: Logger instance
    :param genome_version: genome assembly version for chromosome name normalisation
    :param annot_for_genes: optional dict from utils.prepare_annot_for_genes;
                            enables MANE_Select fallback for unknown transcript versions
    :return: list of BED-format region strings
    """
    # Build lookup indices once for the whole list (avoids O(n*m) repeated scans).
    unversioned_index = build_unversioned_muts_index(annot_for_muts)
    tid_to_gene = get_gene_id_for_transcript(annot_for_genes) if annot_for_genes else {}

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
            bed_coords = get_coords_for_hgvs_c(
                mut.split(':'), flank_length, annot_for_muts, logger,
                unversioned_index=unversioned_index,
                tid_to_gene=tid_to_gene,
                annot_for_genes=annot_for_genes,
            )
            if bed_coords:
                bed_coords_list.append(bed_coords)
        elif len(mut.split(':')) > 1 and 'n.' in mut.split(':')[1]:
            hgvs_counter += 1
            bed_coords = get_coords_for_hgvs_c(
                mut.split(':'), flank_length, annot_for_muts, logger,
                unversioned_index=unversioned_index,
                tid_to_gene=tid_to_gene,
                annot_for_genes=annot_for_genes,
            )
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
        logger.error('ERROR! No mutations were found! Check input files! Exiting...')
        exit(1)

    mut_found_percent = len(bed_coords_list)/(vcf_counter + hgvs_counter)
    logger.info(f'Genome coordinates for {len(bed_coords_list)} were found: processed ({round(mut_found_percent*100)}%) mutations')
    if mut_found_percent < 0.8:
        logger.warning(f'Less than 80% mutations were found in provided annotation! Please check input_muts and genome annotation.')

    if write_intermediate:
        with open(f'{tmpdir}/INPUT_coords_for_mutations.bed', 'w') as outfile:
            logger.info(f"Coordinates for input mutations were written in {tmpdir}/INPUT_coords_for_mutations.bed file")
            outfile.write('\n'.join(bed_coords_list))

    return bed_coords_list

if __name__ == '__main__':
    import utils
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    input_muts_file = Path('../tests/FileTest_input_muts')
    annot_file      = Path('../tests/FileTest_cardio_gencode_annot.gtf')
    tmpdir          = Path('../tests/test_results')
    write_intermediate = True
    flank_length    = 10
    genome_version  = 'GRCh38'

    with open(input_muts_file) as infile:
        input_muts_list = [line.strip() for line in infile.readlines()]

    with open(annot_file) as infile:
        annot_lines = [line.strip() for line in infile.readlines() if not line.startswith('#')]

    # Build both annotation dicts from the same GTF.
    annot_for_muts  = utils.prepare_annot_for_muts(input_muts_list, annot_lines, logger)
    annot_for_genes, gene_name_dict = utils.prepare_annot_for_genes(annot_lines, logger)

    result = get_regions_for_muts(
        input_muts_list, annot_for_muts, flank_length,
        write_intermediate, tmpdir, logger, genome_version,
        annot_for_genes=annot_for_genes,
    )
    print('\n'.join(result))
