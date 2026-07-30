import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.collections as mcollections
from PIL import Image

HEIGHT_PER_OLIG = 0.15
MIN_ROW_PAIR_HEIGHT = 3.0

class Designer():

    def __init__(self, exon_data_path, olig_data_path, result_path, quality_dict={},
                 with_distr=True, subplots_in_row=6, fig_width=5, dpi=600):
        """
        Инициализация объекта Designer для визуализации олигонуклеотидов и профилей качества.
        Высота фигуры вычисляется автоматически из числа олигов на экзон через compute_figsize();
        ширина задаётся параметром fig_width.

        :param exon_data_path: путь до tsv-файла с координатами экзонов
        :param olig_data_path: путь до tsv-файла с координатами олигонуклеотидов
        :param result_path: путь до директории для сохранения PNG
        :param qual_dict: словарь {sequence_id: [float, ...]} Q-значений для каждой позиции окна;
                     если не передан, профили качества не отображаются
        :param with_distr: если True, под каждым экзоном отображается subplot с профилем качества
        :param subplots_in_row: максимальное количество экзонов в одной строке фигуры
        :param fig_width: ширина фигуры в дюймах; высота вычисляется автоматически
        :param dpi: разрешение сохраняемого PNG
        """
        self.exon_cols = ['crh_name', 'ex_start', 'ex_end', 'gene_id', 'ex_id', 'type', 'chr_num']
        self.olig_cols = ['crh_name', 'ol_start', 'ol_end', 'ol_id', 'ex_id']
        plt.rcParams.update({'font.size': 2})
        self.exon_data = exon_data_path
        self.olig_data = olig_data_path
        self.with_distr = with_distr
        self.subplots_in_row = subplots_in_row
        self.fig_width = fig_width
        self.dpi = dpi
        self.result_path = result_path
        self.qual_data = quality_dict
        self.len_lim = 0

    def compute_figsize(self, exons: list, olig_dict: dict) -> tuple:
        """
        Вычисляет размер фигуры (ширина, высота) в дюймах на основе числа олигов на экзон.
        Высота каждой строки сетки пропорциональна максимальному числу олигов среди экзонов
        в этой строке, умноженному на HEIGHT_PER_OLIG; минимум MIN_ROW_PAIR_HEIGHT дюймов.
        :param exons: список экзонов гена в формате [[ex_id, ex_start, ex_end], ...]
        :param olig_dict: dict {ex_id: [[ol_id, ol_start, ol_end], ...]}
        :return: кортеж (fig_width, fig_height) в дюймах
        """
        n_row_pairs = math.ceil(len(exons) / self.subplots_in_row)
        total_height = 0.0

        for row_pair in range(n_row_pairs):
            row_exons = exons[row_pair * self.subplots_in_row :
                              (row_pair + 1) * self.subplots_in_row]
            max_oligs = max(
                len(olig_dict.get(ex[0], [])) for ex in row_exons
            )
            row_height = max(max_oligs * HEIGHT_PER_OLIG, MIN_ROW_PAIR_HEIGHT)
            total_height += row_height

        return (self.fig_width, total_height)
    
    def data_parser(self):
        """
        Читает TSV-файлы экзонов и олигонуклеотидов и строит словари координат.

        Однонуклеотидные gene_id ('A', 'T', 'G', 'C') заменяются на 'SNP'.

        :return: кортеж (gene_ex_cords, ex_olig_cords), где
                 gene_ex_cords — dict {gene_id: [[ex_id, ex_start, ex_end], ...]},
                 ex_olig_cords — dict {ex_id:   [[ol_id, ol_start, ol_end], ...]}
        """
        def csv_reader(path, names):
            return pd.read_csv(path, sep='\t', header=None, names=names)

        exon_data = csv_reader(self.exon_data, self.exon_cols)
        exon_data['gene_id'] = [
            x if (x != 'A') and (x != 'T') and (x != 'G') and (x != 'C') else 'SNP'
            for x in exon_data.gene_id
        ]

        olig_data = csv_reader(self.olig_data, self.olig_cols)

        gene_ex_cords = {}  # gene_id : [[ex_id, ex_start, ex_end], ...]
        for i in exon_data['gene_id'].unique():
            gene_ex_cords[i] = []
            a = exon_data[exon_data['gene_id'] == i]
            for j in a.itertuples(index=False):
                gene_ex_cords[i].append([j.ex_id, j.ex_start, j.ex_end])

        ex_olig_cords = {}  # ex_id : [[ol_id, ol_start, ol_end], ...]
        for i in olig_data['ex_id'].unique():
            ex_olig_cords[i] = []
            a = olig_data[olig_data['ex_id'] == i]
            for j in a.itertuples(index=False):
                ex_olig_cords[i].append([j.ol_id, j.ol_start, j.ol_end])

        return gene_ex_cords, ex_olig_cords
    
    def quality_weighted_coverage(self, olig_list, normalized_qual_list, seq_start) -> float | None:
        """
        Функция вычисляет quality-weighted coverage score для одного таргетного региона (экзона).
        :param olig_list: [[ol_id, ol_start, ol_end], ...] для региона
        :param normalized_qual_list: [float, ...] список с нормализованными Q for sequence_id
        :param seq_start: начальная координата региона (в геноме) (чтобы найти ol_start -> window index)
        :return: quality-weighted coverage score в диапазоне [0, 1], либо None
        """
        if not normalized_qual_list or not olig_list:
            return None

        weighted_sum = 0.0
        total_length = 0

        for olig in olig_list:
            win_idx = olig[1] - seq_start # ol_start -> window index
            length  = olig[2] - olig[1]   # ol_end - ol_start

            if 0 <= win_idx < len(normalized_qual_list):
                weighted_sum += normalized_qual_list[win_idx] * length
                total_length += length

        return weighted_sum / total_length if total_length > 0 else None
    
    def plot_splitter(self, num_of_exomes, coef_sp=2):
        """
        Вычисляет размеры сетки subplots (строки × столбцы) для заданного числа экзонов.

        Экзоны располагаются по столбцам (не более subplots_in_row в ряду).
        Каждый экзон занимает coef_sp строк (при with_distr=True: строка качества +
        строка баров; при with_distr=False: одна строка).
        Число столбцов фиксировано и равно subplots_in_row.

        :param num_of_exomes: общее число экзонов для отображения
        :param coef_sp: число строк на один экзон (2 с профилем качества, 1 без)
        :return: кортеж (nrow, ncol) — число строк и столбцов сетки subplots
        """
        n_row_groups = math.ceil(num_of_exomes / self.subplots_in_row)
        nrow = coef_sp * n_row_groups   # total rows: each group occupies coef_sp rows
        ncol = self.subplots_in_row     # columns fixed to subplots_in_row
        return nrow, ncol

    def exon_subplot_designer(self, exon_name, score, e_start, e_end, olig_list, cordx, cordy, axis):
        """
        Рисует bar-subplot для одного экзона: столбцы олигонуклеотидов и чёрный столбец экзона.

        Если olig_list пуст (нет олигов для данного экзона), рисуется только столбец экзона,
        а self.len_lim остаётся без изменений.

        :param exon_name: идентификатор экзона (используется в заголовке subplot)
        :param e_start: начальная геномная координата экзона
        :param e_end: конечная геномная координата экзона
        :param olig_list: список [[ol_id, ol_start, ol_end], ...] олигов данного экзона;
                          может быть пустым
        :param cordx: индекс столбца в сетке subplots
        :param cordy: индекс строки в сетке subplots
        :param axis: 2-D массив Axes из plt.subplots(squeeze=False)
        """
        curax = axis[cordy, cordx]

        if not olig_list:
            curax.bar("Exon", e_end - e_start, width=0.5, bottom=0,
                      color='black', align='center')
            curax.set_ylim([-40, e_end - e_start + 40])
            curax.tick_params(labelrotation=90, length=1, width=1, pad=1)
            curax.tick_params(axis='y', labelrotation=90)
            curax.set_title(
                f"(no oligos)\n {exon_name.split(';')[0]}\n [{e_start}:{e_end}]",
                rotation=90
            )
            return

        bottom = []
        len_olig = []
        olig_id = []

        for i in olig_list:
            bottom.append(i[1] - e_start)
            len_olig.append(i[2] - i[1])
            olig_id.append(i[0])

        self.len_lim = len_olig[0]

        for i in range(len(bottom)):
            curax.bar(olig_id[i], len_olig[i], width=0.2, bottom=bottom[i], snap=False)

        curax.bar("Exon", e_end - e_start, width=0.5, bottom=0, color='black', align='center')
        curax.set_ylim([-40, e_end - e_start + 40])
        curax.tick_params(labelrotation=90, length=1, width=1, pad=1)
        curax.tick_params(axis='y', labelrotation=90)
        curax.set_title(
            f"{olig_id[0]} : {olig_id[-1]}\n {exon_name.split(';')[0]}\n [{e_start}:{e_end}]\n QWC Score = {score}",
            rotation=90
        )

    def quality_subplot_designer(self, col, row, axis, exon_name, olig_list, e_start):
        """
        Рисует профиль качества для одного экзона в виде градиентной линии.

        Линия разбивается на N-1 сегментов, каждый окрашен по среднему значению
        качества двух соседних точек (цветовая карта RdYlGn: красный — низкое,
        зелёный — высокое качество). На линии отмечаются точки начала каждого олига.

        :param col: индекс столбца в сетке subplots
        :param row: индекс строки в сетке subplots
        :param axis: 2-D массив Axes из plt.subplots(squeeze=False)
        :param exon_name: ключ в self.qual_data для данного экзона
        :param olig_list: список [[ol_id, ol_start, ol_end], ...] олигов данного экзона;
                          используется для размещения маркеров на линии качества
        :param e_start: начальная геномная координата экзона; используется для перевода
                        абсолютных координат олигов в относительные позиции по оси Y
        """
        curax = axis[row, col]
        curax.invert_xaxis()

        xs = self.qual_data[exon_name]          # quality values (normalized 0–1)
        ys = list(range(len(xs)))               # positional indices along the exon

        curax.set_ylim([-40, ys[-1] + 40 + self.len_lim])
        curax.set_xlim([1, 0])
        curax.tick_params(labelrotation=90, length=1, width=1, pad=1)
        curax.tick_params(axis='y', labelrotation=90)

        cmap = plt.get_cmap('RdYlGn')
        norm = mcolors.Normalize(vmin=0, vmax=1)

        points = np.array([xs, ys], dtype=float).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        seg_colors = np.array([(xs[i] + xs[i + 1]) / 2 for i in range(len(xs) - 1)])

        lc = mcollections.LineCollection(
            segments,
            cmap=cmap,
            norm=norm,
            linewidth=0.8,
            zorder=2
        )
        lc.set_array(seg_colors)
        curax.add_collection(lc)

        # Thin colorbar on the right of this subplot
        plt.colorbar(lc, ax=curax, fraction=0.08, pad=0.02, label='Quality')

        # Dot markers at each oligo's start position, interpolated onto the quality line.
        xs_arr = np.array(xs, dtype=float)
        ys_arr = np.array(ys, dtype=float)

        for oligo in olig_list:
            ol_start = oligo[1]
            y_pos = ol_start - e_start          # relative position on y-axis

            # Only draw if the position falls within the quality profile range
            if ys_arr[0] <= y_pos <= ys_arr[-1]:
                x_pos = float(np.interp(y_pos, ys_arr, xs_arr))
                dot_color = cmap(norm(x_pos))
                curax.scatter(
                    x_pos, y_pos,
                    color=dot_color,
                    edgecolors='black',
                    linewidths=0.3,
                    s=4,
                    zorder=5
                )

    def exome_container(self, exome_dict, olig_dict):
        """
        Основной цикл отрисовки: для каждого гена создаёт фигуру с сеткой subplots,
        вызывает exon_subplot_designer и (при with_distr=True) quality_subplot_designer
        для каждого экзона, сохраняет PNG и поворачивает его на 270°.
        Экзоны, отсутствующие в olig_dict, пропускаются с предупреждением в консоль.
        :param exome_dict: dict {gene_id: [[ex_id, ex_start, ex_end], ...]} из data_parser()
        :param olig_dict: dict {ex_id: [[ol_id, ol_start, ol_end], ...]} из data_parser()
        """
        for gene, exons in exome_dict.items():
            gene_name = gene
            coef_sp = 2 if self.with_distr else 1

            plot_rows, plot_cols = self.plot_splitter(len(exons), coef_sp=coef_sp)
            figsize = self.compute_figsize(exons, olig_dict)
            fig, ax = plt.subplots(plot_rows, plot_cols, figsize=figsize, squeeze=False)

            for id in range(len(exons)):
                ex_id = exons[id][0]
                if ex_id not in olig_dict:
                    continue

                ex_start = exons[id][1]
                ex_qual = self.qual_data[ex_id]
                qwc_score = self.quality_weighted_coverage(olig_dict[ex_id], ex_qual, ex_start)
                if qwc_score is None:
                    qwc_score = 0.5
                else:
                    qwc_score = round(qwc_score, 2)

                # cordx: column = which exon within the current row-group
                # cordy: row    = exon bar row within its pair (quality row = cordy-1)
                cordx = id % self.subplots_in_row
                cordy = (id // self.subplots_in_row) * coef_sp + 1 if self.with_distr \
                    else id // self.subplots_in_row

                self.exon_subplot_designer(
                    exon_name=ex_id,
                    score=qwc_score,
                    e_start=exons[id][1],
                    e_end=exons[id][2],
                    olig_list=olig_dict[ex_id],
                    cordx=cordx,
                    cordy=cordy,
                    axis=ax
                )

                if self.with_distr and ex_id in self.qual_data:
                    self.quality_subplot_designer(
                        col=cordx,
                        row=cordy - 1,
                        axis=ax,
                        exon_name=ex_id,
                        olig_list=olig_dict[ex_id],
                        e_start=exons[id][1]
                    )

            # hide empty subplots of the last row-group
            n_used_in_last_group = len(exons) % self.subplots_in_row
            if n_used_in_last_group != 0:
                last_group = (len(exons) - 1) // self.subplots_in_row
                for empty_col in range(n_used_in_last_group, self.subplots_in_row):
                    for row_offset in range(coef_sp):
                        ax[last_group * coef_sp + row_offset, empty_col].set_visible(False)

            fig.tight_layout()
            curfile = rf"{self.result_path}/{gene_name}_exon_oligo_distr.png"
            plt.savefig(curfile, dpi=self.dpi, format='png')
            plt.close(fig)
            self.png_rotater(curfile)

    def png_rotater(self, currimg):
        """
        Поворачивает сохранённый PNG на 270° (эквивалентно повороту на 90° по часовой стрелке)
        и перезаписывает файл.
        :param currimg: путь до PNG-файла для поворота
        """
        ori = Image.open(currimg)
        ori = ori.rotate(270, expand=True)
        ori.save(currimg, quality=100)

    def save(self):
        """
        Точка входа: запускает парсинг данных и полный цикл отрисовки для всех генов.
        Результаты сохраняются в self.result_path.
        """
        exs, oligs = self.data_parser()
        self.exome_container(exs, oligs)
