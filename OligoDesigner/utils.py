from pathlib import Path

def prepare_annot_for_muts(input_muts_list:list, annot:Path) -> dict:
    """
    Эта функция парсит аннотацию в формате gtf и возвращает словарь, необходимый для поиска координат для мутаций
    :param input_muts_list: Список мутаций. Нужен чтобы понять, нужно ли вообще парсить аннотацию или нет
    :param annot: Объект Path с путем к аннотации в формате gtf
    :return: Словарь с информацией по каждому транскрипту, необходимой для поиска фланков мутаций
    """
    prepare_annot_for_mut_flag = 0
    for mut in input_muts_list:
        if len(mut.split(':')) > 1:
            if 'c.' in mut.split(':')[1] or 'n.' in mut.split(':')[1]:
                prepare_annot_for_mut_flag = 1

    if not prepare_annot_for_mut_flag:
         return None

    transcript_dict = {}
    for i in range(len(annot)):
        strand = annot[i].split()[6]
        if annot[i].split('\t')[2] == 'transcript':
            transcript_id = annot[i].split('";')[2].split('"')[1] + '.' + annot[i].split('";')[3].split('"')[1]
            transcript_chr = annot[i].split('\t')[0]
            exon_list = []
            k = 1
            start_codon_coord = 0
            while True:
                if (i + k == (len(annot))) or (annot[i+k].split('\t')[2] == 'transcript'):
                    break
                elif annot[i+k].split('\t')[2] == 'exon':
                    exon_list.append([annot[i+k].split('\t')[3], annot[i+k].split('\t')[4]])
                elif annot[i+k].split('\t')[2] == 'start_codon':
                    start_codon_coord = annot[i+k].split('\t')[3]
                k += 1
            if start_codon_coord:
                transcript_dict[transcript_id] = [transcript_chr, exon_list, start_codon_coord, strand]
                
    return transcript_dict

def prepare_annot_for_genes(annot:Path, logger) -> [dict, dict]:
    """
    Эта функция парсит аннотацию в формате gtf и возвращает словарь, необходимый для поиска координат для генов
    :param annot: Объект Path с путем к аннотации в формате gtf
    :return: Словарь с информацией по каждому гену, необходимой для поиска регионов 
             и словарь с соответствием символа гена и его идентификатора
    """
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
            stop_codon_coords = []
            while True:
                if (i + k == (len(annot))) or (annot[i+k].split('\t')[2] == 'gene'):
                    if exon_list:
                        exon_list = sorted(exon_list, key=lambda x: x[0])
                        if stop_codon_coords:
                            if strand == '+':
                                #logger.error(stop_codon_coords)
                                for stop in stop_codon_coords:
                                    for exon in range(len(exon_list)):
                                        if exon_list[exon][0] <= stop <= exon_list[exon][1]:
                                            exon_list[exon][1] = stop
                            elif strand == '-':
                                #logger.error(stop_codon_coords)
                                for stop in stop_codon_coords:
                                    for exon in range(len(exon_list)):
                                        if exon_list[exon][0] <= stop <= exon_list[exon][1]:
                                            exon_list[exon][0] = stop
                    break

                elif annot[i+k].split('\t')[2] == 'exon' and annot[i+k].split('transcript_biotype "')[1].split('";')[0] == 'protein_coding':
                    strand = annot[i+k].split('\t')[6]
                    exon_id = annot[i+k].split('exon_id "')[1].split('";')[0]
                    if exon_id not in exon_id_list:
                        exon_list.append([int(annot[i+k].split('\t')[3]), int(annot[i+k].split('\t')[4])])
                        exon_id_list.append(exon_id)
                elif annot[i+k].split('\t')[2] == 'five_prime_utr' or annot[i+k].split('\t')[2] == 'three_prime_utr':
                    utr_coords = [annot[i+k].split('\t')[3], annot[i+k].split('\t')[4]]
                    if utr_coords not in utr_id_list:
                        utr_list.append(utr_coords)
                elif annot[i+k].split('\t')[2] == 'stop_codon':
                    stop_codon_coords.append(int(annot[i+k].split('\t')[4]))
                    
                k += 1

            logger.warning(f'{gene_id} {exon_list}')

            gene_dict[gene_id] = [gene_chr, gene_coords, exon_list, utr_list, exon_id_list]
            
    return gene_dict, gene_name_dict

if __name__ == '__main__':
    test_gtf = './tests/FileTest_annot.gtf'
    with open(test_gtf) as infile:
        annot = [line.strip() for line in infile.readlines() if not line.startswith('#')]

    annot_test = prepare_annot_for_muts(['Chr2:c.123G>C'], annot, '')

    print(annot_test)
