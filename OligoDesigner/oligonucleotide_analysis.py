import logging
import pandas as pd
from oligo_plot import Designer

def get_genome_coords_for_oligos(oligos_df, input_fasta_list, tmpdir, logger):
    """
    Эта функция находит геномные координаты для каждого олига

    """
    exon_seq_dict, exon_coord_dict = {}, {}
    for input_fasta in input_fasta_list:
        if 'start' in input_fasta.description:
            exon_id = input_fasta.id
            start_coord = input_fasta.description.split('start:')[1].split(' end')[0]
            chrom = input_fasta.description.split('chrom:')[1].split(' start')[0]
            sequence = input_fasta.seq

            exon_coord_dict[exon_id] = [chrom, start_coord]
            exon_seq_dict[exon_id] = sequence

    olig_bed_list = []
    for _, olig in oligos_df.iterrows():
        if olig.Region_Name in exon_seq_dict:
            olig_start_coord = exon_seq_dict[olig.Region_Name].index(olig.Seq) + int(exon_coord_dict[olig.Region_Name][1])
            olig_end_coord = olig_start_coord + int(olig.Length)
            olig_chr = exon_coord_dict[olig.Region_Name][0]
            olig_bed_list.append([olig_chr, str(olig_start_coord), str(olig_end_coord), olig.Olig_id, olig.Region_Name])

    if olig_bed_list:
        with open(tmpdir/'OUTPUT_Oligonucleotides.bed', 'w') as outfile:
            logger.info(f"Oligunucleotide genome coordinates were written in {tmpdir/'OUTPUT_Oligonucleotides.bed'} file")
            outfile.write('\n'.join(['\t'.join(l) for l in olig_bed_list]))

    return olig_bed_list

def analyze_oligos(oligo_df:pd.DataFrame, input_fasta_list:list, qual_dict:dict, tmpdir, logger:logging, outdir) -> None:
    """
    Эта функция выполняет анализ разработанной панели олигонуклеотидов
    :param oligo_df: Датафрейм с олигонуклеотидами
    :param logger: Объект класса logging
    :return: Ничего не возвращает
    """
    oligo_df[["GC", "Entropy", "Tm", "Repeats", "Masked", "Blast_Hits"]] = \
        oligo_df[["GC", "Entropy", "Tm", "Repeats", "Masked", "Blast_Hits"]].apply(pd.to_numeric)

    logger.info(f'The mean GC content in designed oligos is {round(oligo_df["GC"].mean(), 3)}')
    logger.info(f'The mean Shannon Entropy in designed oligos is {round(oligo_df["Entropy"].mean(), 3)}')
    logger.info(f'The mean Melting Temperature in designed oligos is {round(oligo_df["Tm"].mean(), 3)}')
    logger.info(f'The mean Repeat Length in designed oligos is {round(oligo_df["Repeats"].mean(), 3)}')
    logger.info(f'The mean Masked Length in designed oligos is {round(oligo_df["Masked"].mean(), 3)}')
    logger.info(f'The mean Blast Hits in designed oligos is {round(oligo_df["Blast_Hits"].mean(), 3)}')

    olig_bed_list = get_genome_coords_for_oligos(oligo_df, input_fasta_list, tmpdir, logger)

    if olig_bed_list:
        des = Designer(tmpdir/'INPUT_coords_FIXED.bed', tmpdir/'OUTPUT_Oligonucleotides.bed', qual = qual_dict, result_path=outdir)
        des.save()
