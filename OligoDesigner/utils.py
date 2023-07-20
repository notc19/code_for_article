def prepare_annot_for_muts(input_muts_list, annot, logger):
    prepare_annot_for_mut_flag = 0
    for mut in input_muts_list:
        if len(mut.split(':')) > 1:
            if 'c.' in mut.split(':')[1] or 'n.' in mut.split(':')[1]:
                prepare_annot_for_mut_flag = 1

    if not prepare_annot_for_mut_flag:
         return None

    #logger.info('Preparing annotation file for mutations')
    transcript_dict = {}
    for i in range(len(annot)):
        if annot[i].split('\t')[2] == 'transcript':
            transcript_id = annot[i].split('";')[2].split('"')[1] + '.' + annot[i].split('";')[3].split('"')[1]
            transcript_chr = annot[i].split('\t')[0]
            exon_list = []
            k = 1
            while True:
                if (i + k == (len(annot))) or (annot[i+k].split('\t')[2] == 'transcript'):
                    break
                elif annot[i+k].split('\t')[2] == 'exon':
                    exon_list.append([annot[i+k].split('\t')[3], annot[i+k].split('\t')[4]])
                k += 1
            transcript_dict[transcript_id] = [transcript_chr, exon_list]
                
    return transcript_dict

def prepare_annot_for_genes(annot, logger):

    gene_dict, gene_name_dict = {}, {}
    for i in range(len(annot)):
        if annot[i].split('\t')[2] == 'gene':
            gene_chr = annot[i].split('\t')[0]
            gene_coords = [annot[i].split()[3], annot[i].split()[4]]
            gene_id = annot[i].split('";')[0].split('"')[1]
            gene_name = annot[i].split('";')[2].split('"')[1]

            gene_name_dict[gene_name] = gene_id
            exon_list, utr_list = [], []
            exon_id_list, utr_id_list = [], []
            k = 1
            while True:
                if (i + k == (len(annot))) or (annot[i+k].split('\t')[2] == 'gene'):
                    break
                elif annot[i+k].split('\t')[2] == 'exon':
                    exon_id = annot[i+k].split('exon_id "')[1].split('";')[0]
                    if exon_id not in exon_id_list:
                        exon_list.append([annot[i+k].split('\t')[3], annot[i+k].split('\t')[4]])
                        exon_id_list.append(exon_id)
                elif annot[i+k].split('\t')[2] == 'five_prime_utr' or annot[i+k].split('\t')[2] == 'three_prime_utr':
                    utr_coords = [annot[i+k].split('\t')[3], annot[i+k].split('\t')[4]]
                    if utr_coords not in utr_id_list:
                        utr_list.append(utr_coords)

                k += 1
            gene_dict[gene_id] = [gene_chr, gene_coords, exon_list, utr_list, exon_id_list]
            
    return gene_dict, gene_name_dict

if __name__ == '__main__':
    test_gtf = '/home/a_vikhorev/OligoDesigner/tests/FileTest_annot.gtf' #'/home/a_vikhorev/Data/Genome/Homo_sapiens.GRCh38.109.gtf'
    with open(test_gtf) as infile:
        annot = [line.strip() for line in infile.readlines() if not line.startswith('#')]

    d,_ = prepare_annot_for_genes(annot, '') #/home/a_vikhorev/OligoDesigner/tests/test_annot.gtf
    print(d)
    print(_)
