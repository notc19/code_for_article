import pandas as pd 
import numpy as np
import math 
import matplotlib.pyplot as plt
from PIL import Image

class Designer():

    def __init__(self,exon_data_path, olig_data_path, result_path, qual = [], \
                  with_distr = True, subplots_in_row = 5, figsize = (5,15), dpi = 1000):
        self.exon_cols = ['crh_name','ex_start', 'ex_end', 'gene_id', 'ex_id', 'type', 'chr_num']
        self.olig_cols = ['crh_name', 'ol_start', 'ol_end', 'ol_id', 'ex_id']
        plt.rcParams.update({'font.size': 2})
        self.exon_data = exon_data_path
        self.olig_data = olig_data_path
        self.with_distr = with_distr
        self.subplots_in_row = subplots_in_row
        self.figsize = figsize # tuple(x,y)
        self.dpi = dpi
        self.result_path = result_path
        self.qual_data = qual

    def data_parser(self):

        def csv_reader(path, names):
            return pd.read_csv(path, sep='\t', header=None, names=names)
        
        
        exon_data = csv_reader(self.exon_data, self.exon_cols)
        exon_data['gene_id'] = [x if (x!='A') and (x!='T') and (x!='G') and (x!='C') else 'SNP' for x in exon_data.gene_id]

        olig_data = csv_reader(self.olig_data, self.olig_cols)
        
        gene_ex_cords = {} # gene_id : [ex_id, ex_start, ex_end]

        for i in exon_data['gene_id'].unique():
            gene_ex_cords[i] = []
            a = exon_data[exon_data['gene_id'] == i]
            for j in a.iloc():
                gene_ex_cords[i].append([j['ex_id'],j['ex_start'],j['ex_end']])
                # print(j)


        ex_olig_cords = {} # ex_id : [ol_id, ol_start, ol_end]

        for i in olig_data['ex_id'].unique():
            ex_olig_cords[i] = []
            a = olig_data[olig_data['ex_id'] == i]
            for j in a.iloc():
                ex_olig_cords[i].append([j[3],j[1],j[2]])
        
        return gene_ex_cords, ex_olig_cords
    
    def plot_splitter(self, num_of_exomes, coef_sp = 2):
        nrow = num_of_exomes
        nfloor = coef_sp
        if nrow > self.subplots_in_row:
            nfloor = coef_sp*math.ceil(nrow/self.subplots_in_row)
            nrow = self.subplots_in_row
        return nfloor, nrow

    def exon_subplot_designer(self, exon_name, e_start, e_end, olig_list, cordx, cordy, axis):
        bottom = []
        len_olig = []
        olig_id = []

        curax = axis[cordy, cordx]

        for i in olig_list:
            bottom.append(i[1]-e_start)
            len_olig.append(i[2]-i[1])
            olig_id.append(i[0])

        for i in range(len(bottom)):
            curax.bar(olig_id[i],len_olig[i], width = 0.2, bottom = bottom[i], snap=False)
            # bottom[i]+=len_olig[i]
        # curax_dis = axis[cordy-1,cordx]

        curax.bar("Exon", e_end - e_start, width = 0.5, bottom = 0, color = 'black', align = 'center')
        curax.set_ylim([-40, e_end-e_start+40])
        curax.tick_params(labelrotation=90, length = 1, width = 1, pad = 1)
        curax.tick_params(axis = 'y', labelrotation=90)
        # curax.set_yticklabels([e_start, e_start + (e_end-e_start)/2, e_end])
        curax.set_title(f"{olig_id[0]} : {olig_id[-1]}\n {exon_name.split(';')[0]}\n [{e_start}:{e_end}]", rotation = 90)

    def quality_subplot_designer(self, qx, qy, axis, exon_name):
        curax = axis[qx,qy]
        curax.invert_xaxis()
        xs = self.qual_data[exon_name]
        ys = list(range(len(xs)))
        #for i in range(len(xs)):
        #    ys[i] = ys[i]*10
        curax.set_ylim([-40, ys[-1]+40])
        curax.tick_params(labelrotation=90, length = 1, width = 1, pad = 1)
        curax.tick_params(axis = 'y', labelrotation=90)
        
        curax.plot(xs,ys)
    
    def exome_container(self, exome_dict, olig_dict):
        for gene, exons in exome_dict.items():
            gene_name = gene
            coef_sp = 2

            if not self.with_distr:
                coef_sp = 1
            
            plot_floors, plot_rows = self.plot_splitter(len(exons), coef_sp=coef_sp)
            fig, ax = plt.subplots(plot_rows,plot_floors, figsize = (5,15))
            
            for id in range(len(exons)):
                cordx = id%self.subplots_in_row
                cordy = id//self.subplots_in_row*2+1 if self.with_distr else id//self.subplots_in_row

                if self.with_distr:
                    qx = cordx 
                    qy = cordy - 1
                    self.quality_subplot_designer(qx = qx, qy = qy, axis = ax, exon_name = exons[id][0])

                # zaplatka(exon[0],exon[1:3], olig_dict[exon[0]])
                self.exon_subplot_designer(exon_name=exons[id][0], e_start=exons[id][1], e_end=exons[id][2],\
                                            olig_list=olig_dict[exons[id][0]], cordx=cordy, cordy=cordx, axis = ax)
            fig.tight_layout()
            # fig.suptitle(gene_name, fontsize = '5')
            curfile = rf"{self.result_path}/{gene_name}_exon_oligo_distr.png"
            plt.savefig(curfile, dpi = self.dpi, format = 'png')  
            self.png_rotater(curfile)
        
    def png_rotater(self, currimg):
        ori = Image.open(currimg)
        ori = ori.rotate(270, expand=True)
        ori.save(currimg, quality = 100)


    def save(self):
        exs, oligs = self.data_parser() #insert paths   
        self.exome_container(exs,oligs)      

