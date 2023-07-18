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
    with open(annot) as infile:
        input = [line.strip() for line in infile.readlines() if not line.startswith('#')]
        for i in range(len(input)):
            if input[i].split('\t')[2] == 'transcript':
                transcript_id = input[i].split('";')[2].split('"')[1] + '.' + input[i].split('";')[3].split('"')[1]
                transcript_chr = input[i].split('\t')[0]
                exon_list = []
                k = 1
                while True:
                    if (i + k == len(input)) or (input[i+k].split('\t')[2] == 'transcript'):
                        break
                    elif input[i+k].split('\t')[2] == 'exon':
                        exon_list.append([input[i+k].split('\t')[3], input[i+k].split('\t')[4]])
                    k += 1
                transcript_dict[transcript_id] = [transcript_chr, exon_list]

    return transcript_dict

def prepare_annot_for_genes(input_genes_list, annot, logger):
    pass

if __name__ == '__main__':
    d = prepare_annot_for_muts(['Chr2:c.123G>C'], '/home/a_vikhorev/Data/Genome/Homo_sapiens.GRCh38.109.gtf', '')
    print(d['ENST00000641515.2'])
