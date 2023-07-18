genome_versions_chromosome_dict = {
    'GRCh38.p14': {
        'NC_000001.11': 'chr1',
        'NC_000002.12': 'chr2',
        'NC_000003.12': 'chr3',
        'NC_000004.12': 'chr4',
        'NC_000005.10': 'chr5',
        'NC_000006.12': 'chr6',
        'NC_000007.14': 'chr7',
        'NC_000008.11': 'chr8',
        'NC_000009.12': 'chr9',
        'NC_000010.11': 'chr10',
        'NC_000011.10': 'chr11',
        'NC_000012.12': 'chr12',
        'NC_000013.11': 'chr13',
        'NC_000014.9': 'chr14',
        'NC_000015.10': 'chr15',
        'NC_000016.10': 'chr16',
        'NC_000017.11': 'chr17',
        'NC_000018.10': 'chr18',
        'NC_000019.10': 'chr19',
        'NC_000020.11': 'chr20',
        'NC_000021.9': 'chr21',
        'NC_000022.11': 'chr22',
        'NC_000023.11': 'chrX',
        'NC_000024.10': 'chrY'
    },
    'GRCh37.p13': {
        'NC_000001.10': 'chr1',
        'NC_000002.11': 'chr2',
        'NC_000003.11': 'chr3',
        'NC_000004.11': 'chr4',
        'NC_000005.9': 'chr5',
        'NC_000006.11': 'chr6',
        'NC_000007.13': 'chr7',
        'NC_000008.10': 'chr8',
        'NC_000009.11': 'chr9',
        'NC_000010.10': 'chr10',
        'NC_000011.9': 'chr11',
        'NC_000012.11': 'chr12',
        'NC_000013.10': 'chr13',
        'NC_000014.8': 'chr14',
        'NC_000015.9': 'chr15',
        'NC_000016.9': 'chr16',
        'NC_000017.10': 'chr17',
        'NC_000018.9': 'chr18',
        'NC_000019.9': 'chr19',
        'NC_000020.10': 'chr20',
        'NC_000021.8': 'chr21',
        'NC_000022.10': 'chr22',
        'NC_000023.10': 'chrX',
        'NC_000024.9': 'chrY'
    }
}

def get_coords_for_vcf(vcf:str, flank_length:int):
    """
    """
    chrom = vcf[0]
    rsid = vcf[1]
    pos = int(vcf[2])
    ref = vcf[3]

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref + '\t' + rsid

def get_coords_for_hgvs_g(hgvsg:str, flank_length:int):
    """
    """
    chrom = hgvsg[0]
    pos_and_nucl = hgvsg[1].split('.')[1] if '.' in hgvsg[1] else hgvsg[1]
    i = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
    pos = int(pos_and_nucl[:i])
    ref = pos_and_nucl[i:].split('>')[0]

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref

def fix_chromosome_name(coord_list, genome_version, logger):
    """
    """
    fixed_coord_list = []
    for coord in coord_list:
        coord = coord.split()
        chr_in_input = coord[0]
        if coord[0][:3] != 'chr':
            if coord[0].isdigit():
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
        fixed_coord_list.append('\t'.join([fixed_chr, coord[1], coord[2], coord[3], chr_in_input]))
    
    return fixed_coord_list

def get_coords_for_hgvs_c(hgvsc:str, flank_length:int, annot, logger):
    
    if hgvsc[0] in annot:
        trans = hgvsc[0]
        chrom = annot[trans][0]
        exons = annot[trans][1]
        pos_and_nucl = hgvsc[1].split('.')[1]
        if not '-' in hgvsc[1] and not '+' in hgvsc[1] and not '*' in hgvsc[1]:
            j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
            ref = pos_and_nucl[j:].split('>')[0]
            pos_c = int(pos_and_nucl[:j])
            if '_' in str(pos_and_nucl[:j]):
                pos_c = int(str(pos_and_nucl[:j]).split('_')[0])
            i = 0
            while True:
                try:
                    if pos_c > (int(exons[i][1]) - int(exons[i][0])):
                        pos_c -= (int(exons[i][1]) - int(exons[i][0]))

                    elif pos_c <= (int(exons[i][1]) - int(exons[i][0])):
                        pos = int(exons[i][0]) + pos_c
                        break
                except IndexError:
                    print(j, pos_and_nucl[:j], pos_c, pos_and_nucl)
                i += 1
        elif pos_and_nucl.startswith('-'):
            j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
            ref = pos_and_nucl[j:].split('>')[0]
            pos = int(exons[0][0]) - int(pos_and_nucl[1:j])
        elif pos_and_nucl.startswith('*'):
            j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
            ref = pos_and_nucl[j:].split('>')[0]
            pos = int(exons[1][1]) + int(pos_and_nucl[1:j])
        else:
            logger.warning(f'This format of HGVS mutation ({trans +":"+hgvsc[1]}) is not supproted yet. Skipping this variant')
            return None
    else:
        logger.warning(f'Transcript {hgvsc[0]} not present in provided annotation!!! Skipping {hgvsc[0] +":"+hgvsc[1]} variant!')
        return None

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref

def get_regions_for_muts(muts_list:list, annot:dict, flank_length:int, write_intermediate:bool, tmpdir:str, \
                         logger, genome_version:str):
    """
    """
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
            bed_coords = get_coords_for_hgvs_c(mut.split(':'), flank_length, annot, logger)
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

    # HERE TO FIX CHROMOSOME NAMES
    if genome_version in ['GRCh38.p14', 'GRCh37.p13']:
        bed_coords_list = fix_chromosome_name(bed_coords_list, genome_version, logger)

    if write_intermediate:
        with open(tmpdir+'/INPUT_coords_for_mutations.bed', 'w') as outfile:
            logger.info(f"Coordinates for input mutations were written in {tmpdir+'/INPUT_coords_for_mutations.bed'} file")
            outfile.write('\n'.join(bed_coords_list))

    return bed_coords_list

def get_regions_for_genes(genes_list, annot):
    return genes_list

def get_fasta_seqs_from_bed(coord_list, genome, genome_version):
    """
    """



    return coord_list

if __name__ == '__main__':
    #print(get_regions_for_muts(['#', 'chr1\t.\t100000\tA'],'asd',85, True, './tmpdir'))

    hgvsc = ['test','c.209A>T']
    flank_length = 100
    annot = {'test':['1', [[10,90],[100,200],[400,500]]]}
    print(get_coords_for_hgvs_c(hgvsc, flank_length, annot, ""))