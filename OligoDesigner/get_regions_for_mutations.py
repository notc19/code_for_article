from constants import genome_versions_chromosome_dict

def get_coords_for_vcf(vcf:list, flank_length:int) -> str:
    """
    """
    chrom = vcf[0]
    rsid = vcf[1]
    pos = int(vcf[2])
    ref = vcf[3]

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref + '\t' + rsid

def get_coords_for_hgvs_g(hgvsg:list, flank_length:int) -> str:
    """
    """
    chrom = hgvsg[0]
    pos_and_nucl = hgvsg[1].split('.')[1] if '.' in hgvsg[1] else hgvsg[1]
    i = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
    pos = int(pos_and_nucl[:i])
    ref = pos_and_nucl[i:].split('>')[0]

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref + '\t' + ':'.join(hgvsg)

def get_coords_for_hgvs_c_intr(pos_and_nucl, exons, sign):
    pos_c = int(pos_and_nucl.split(sign)[0])
    j = len(['_' for c in range(len(pos_and_nucl.split(sign)[1])) if pos_and_nucl.split(sign)[1][c].isdigit()]) # index of the last digit +1
    ref = pos_and_nucl.split(sign)[1][j:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
        else 'del'+pos_and_nucl.split('del')[1] if 'ins' not in pos_and_nucl \
        else 'ins'+pos_and_nucl.split('ins')[1]
    pos_m = int(pos_and_nucl.split(sign)[1][:j])
    if '_' in str(pos_and_nucl):
        pos_c = int(str(pos_and_nucl).split('_')[0].split(sign)[0])
        pos_m = int(str(pos_and_nucl).split('_')[0].split(sign)[1])
    i = 0
    while True:
        if pos_c > (int(exons[i][1]) - int(exons[i][0])):
            pos_c -= (int(exons[i][1]) - int(exons[i][0]))
            i += 1
        elif pos_c <= (int(exons[i][1]) - int(exons[i][0])):
            pos = int(exons[i][0]) + pos_c - pos_m if sign == '-' else int(exons[i][0]) + pos_c + pos_m
            break
    return pos, ref

def get_coords_for_hgvs_c_utr_intr(pos_and_nucl, exons, sign, orient):
    pos_c = int(pos_and_nucl[1:].split(sign)[0])
    j = len(['_' for c in range(len(pos_and_nucl[1:].split(sign)[1])) if pos_and_nucl[1:].split(sign)[1][c].isdigit()]) # index of the last digit +1
    ref = pos_and_nucl[1:].split(sign)[1][j:].split('>')[0] if 'del' not in pos_and_nucl[1:] and 'ins' not in pos_and_nucl[1:] \
        else 'del'+pos_and_nucl[1:].split('del')[1] if 'ins' not in pos_and_nucl[1:] \
        else 'ins'+pos_and_nucl[1:].split('ins')[1]
    pos_m = int(pos_and_nucl[1:].split(sign)[1][:j])

    if '_' in str(pos_and_nucl[1:]):
        pos_c = int(str(pos_and_nucl[1:]).split('_')[0].split(sign)[0])
        pos_m = int(str(pos_and_nucl[1:]).split('_')[0].split(sign)[1])
        print(pos_c, pos_m)
    pos_m *= -1 if sign == '-' else 1
    pos = int(exons[-1][1]) + pos_c + pos_m if orient == '*' else int(exons[0][0]) - pos_c + pos_m
    print(pos_c, pos_m, pos)
    return pos, ref

def get_coords_for_hgvs_c(hgvsc:str, flank_length:int, annot, logger):
    
    if hgvsc[0] in annot:
        trans = hgvsc[0]
        chrom = annot[trans][0]
        exons = annot[trans][1]
        pos_and_nucl = hgvsc[1].split('.')[1]

        if pos_and_nucl[0].isdigit():
            if not '-' in pos_and_nucl and not '+' in pos_and_nucl:
                j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
                ref = pos_and_nucl[j:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
                    else 'del'+pos_and_nucl[1:].split('del')[1] if 'ins' not in pos_and_nucl \
                    else 'ins'+pos_and_nucl.split('ins')[1]
                pos_c = int(pos_and_nucl[:j])
                if '_' in str(pos_and_nucl[:j]):
                    pos_c = int(str(pos_and_nucl[:j]).split('_')[0])
                i = 0
                while True:
                    if pos_c > (int(exons[i][1]) - int(exons[i][0])):
                        pos_c -= (int(exons[i][1]) - int(exons[i][0]))
                        i += 1
                    elif pos_c <= (int(exons[i][1]) - int(exons[i][0])):
                        pos = int(exons[i][0]) + pos_c
                        break
                    
            elif '-' in pos_and_nucl:
                pos, ref = get_coords_for_hgvs_c_intr(pos_and_nucl, exons, '-')

            elif '+' in pos_and_nucl:
                pos, ref = get_coords_for_hgvs_c_intr(pos_and_nucl, exons, '+')

        elif pos_and_nucl.startswith('-'):
            if not '-' in pos_and_nucl[1:] and not '+' in pos_and_nucl[1:]:
                j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
                ref = pos_and_nucl[j+1:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
                    else 'del'+pos_and_nucl[1:].split('del')[1] if 'ins' not in pos_and_nucl \
                    else 'ins'+pos_and_nucl.split('ins')[1]
                if '_' in str(pos_and_nucl):
                    pos_c = int(str(pos_and_nucl[1:]).split('_')[0])
                else:
                    pos_c = int(pos_and_nucl[1:j+1])
                pos = int(exons[0][0]) - pos_c
            elif '-' in pos_and_nucl[1:]:
                pos, ref = get_coords_for_hgvs_c_utr_intr(pos_and_nucl, exons, '-', '-')
            elif '+' in pos_and_nucl[1:]:
                pos, ref = get_coords_for_hgvs_c_utr_intr(pos_and_nucl, exons, '+', '-')

        elif pos_and_nucl.startswith('*'):
            if not '-' in pos_and_nucl[1:] and not '+' in pos_and_nucl[1:]:
                j = len(['_' for c in range(len(pos_and_nucl)) if pos_and_nucl[c].isdigit()]) # index of the last digit +1
                ref = pos_and_nucl[j+1:].split('>')[0] if 'del' not in pos_and_nucl and 'ins' not in pos_and_nucl \
                    else 'del'+pos_and_nucl[1:].split('del')[1] if 'ins' not in pos_and_nucl \
                    else 'ins'+pos_and_nucl.split('ins')[1]
                if '_' in str(pos_and_nucl):
                    pos_c = int(str(pos_and_nucl[1:]).split('_')[0])
                else:
                    pos_c = int(pos_and_nucl[1:j+1])
                pos = int(exons[-1][1]) + pos_c
            elif '-' in pos_and_nucl[1:]:
                pos, ref = get_coords_for_hgvs_c_utr_intr(pos_and_nucl, exons, '-', '*')
            elif '+' in pos_and_nucl[1:]:
                pos, ref = get_coords_for_hgvs_c_utr_intr(pos_and_nucl, exons, '+', '*')
        else:
            logger.warning(f'This format of HGVS mutation ({trans +":"+hgvsc[1]}) is not supproted yet. Skipping this variant')
            return None
    else:
        logger.warning(f'Transcript {hgvsc[0]} not present in provided annotation!!! Skipping {hgvsc[0] +":"+hgvsc[1]} variant!')
        return None

    start = pos - flank_length
    end = pos + flank_length

    return chrom + '\t' + str(start) + '\t' + str(end) + '\t' + ref + '\t' + ':'.join(hgvsc)

def fix_chromosome_name_and_check_coords(coord_list:list, genome_version:str, logger) -> list:
    """
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
            fixed_coord_list.append('\t'.join([fixed_chr, coord[1], coord[2], coord[3], coord[4], chr_in_input]))
        except:
            logger.warning(f'Something wrong with {coord} variant! It has no ref allele! Skipping this variant!')
    
    return fixed_coord_list

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

    # fix chromosome names for human genome assemblies
    if genome_version in ['GRCh38.p14', 'GRCh37.p13']:
        bed_coords_list = fix_chromosome_name_and_check_coords(bed_coords_list, genome_version, logger)

    if not bed_coords_list:
        logger.error('ERROR! No mutation were found! Check input files! Exiting...')
        exit(1)

    mut_found_percent = len(bed_coords_list)/(vcf_counter + hgvs_counter)
    logger.info(f'Genome coordinates were found for {len(bed_coords_list)} ({round(mut_found_percent*100)}%) mutations')
    if mut_found_percent < 0.8:
        logger.warning(f'Too little input mutations have coordinates in provided genome assembly! Please check if you use right genome assmbly')

    if write_intermediate:
        with open(tmpdir+'/INPUT_coords_for_mutations.bed', 'w') as outfile:
            logger.info(f"Coordinates for input mutations were written in {tmpdir+'/INPUT_coords_for_mutations.bed'} file")
            outfile.write('\n'.join(bed_coords_list))

    return bed_coords_list


if __name__ == '__main__':
    # Simple test
    hgvsc = ['test', 'c.*10+3A>G'] 
    flank_length = 100
    annot = {'test':['1', [[10,90],[100,200],[400,500]]]}
    print(get_coords_for_hgvs_c(hgvsc, flank_length, annot, ""))