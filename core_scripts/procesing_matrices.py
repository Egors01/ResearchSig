import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt

from constants import MATRICES_PATH, PATH_TO_PLOTS
from core_scripts.matrix_operations import normalize_matrix_to_type_sum, \
    reduce_192_to_96_common_notation, reduce_192_to_96_kelley_notation, \
    create_ratio_table, create_heatmap_df_for_species, \
    normalize_frequencies_192_context, matrix_to_r_output


def plot_heatmaps_from_ratios_list(ratios_table, ratios_to_plot, title='',
                                   kelley_notation=True,if_save=True,filename=''):
    heatmaps_to_plot = [
        create_heatmap_df_for_species(relations_df=ratios_table,
                                             column_name=x,
                                      kelley_notation=kelley_notation)
        for x in ratios_to_plot]

    fig = plt.figure(figsize=(15, 10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    plt.title(title, y=1.08, fontsize=18)

    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    fig.axes[0].get_xaxis().set_visible(False)
    fig.axes[0].get_yaxis().set_visible(False)
    fig.axes[0].set_frame_on(False)
    for i, heatmap in enumerate(heatmaps_to_plot):
        ax = fig.add_subplot(1, len(heatmaps_to_plot), i + 1)
        # ax.text(0.5, 0.5, str((1, 4, i)), fontsize=18, ha='center')
        sns.heatmap(heatmap, annot=True, annot_kws={"size": 9}, linewidths=.6,
                    center=1,vmin=0.5,vmax=1.5, cmap='RdBu_r', ax=ax)
        fig.subplots_adjust(hspace=0.8)
        ax.set_title(ratios_to_plot[i])

        #fig.subplots_adjust(top=0.5)
    #fig.subplots_adjust(top=.5)
    #fig.suptitle(title, size=12)
    #fig.subplots_adjust(hspace=0.4, wspace=0.4)
    #fig.subplots_adjust(top=.5)
    # plt.tight_layout()
    # plt.subplots_adjust(top=1.15)
    if if_save:
        plt.savefig(os.path.join(PATH_TO_PLOTS,filename+'.png'))

    plt.show()

import matplotlib.backends.backend_pdf

def plot_mutational_spectra(matrix, species_names, title='',
                            kelley_notation=False,
                            save_to_pdf=False):
    if save_to_pdf:
        pp = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
    for species_name in species_names:
        if not kelley_notation:
            if len(matrix['SUBS'].values[:]) > 96:
                g = sns.FacetGrid(matrix, col="SUBS", hue="SUBS", col_wrap=6,
                                  col_order=['C>A', 'C>G', 'C>T', 'T>A', "T>C",
                                             "T>G",'G>T', 'G>T', 'G>A', 'A>T', "A>G",
                                             "A>C"])
            else:
                g = sns.FacetGrid(matrix, col="SUBS", hue="SUBS", col_wrap=6,
                                  col_order=['C>A', 'C>G', 'C>T', 'T>A', "T>C",
                                             "T>G"])
        if kelley_notation:
            g = sns.FacetGrid(matrix, col="SUBS", hue="SUBS", col_wrap=6,
                              col_order=['C>A', 'C>G', 'C>T', 'A>G', "A>C", "A>T"])
        g.map(sns.barplot, "Context", species_name)
        g.set_xticklabels(rotation=90)
        g.set_titles("{col_name} ")
        g.fig.subplots_adjust(top=.5)
        g.fig.suptitle(title, size=16)
        if save_to_pdf:
            plt.savefig(pp, format='pdf')
        else:
            plt.show()
    if save_to_pdf:
        pp.close()





def run_process():
    m192_raw_filename = os.path.join(MATRICES_PATH, '192_matrix_29_all_recent.csv')



    #roll up and normalize 96 ONLY by sum as in paper
    m96_kelley_notation_filename = \
        reduce_192_to_96_kelley_notation(path_to_192_matrix=m192_raw_filename)
    m96_kelley_not_norm = pd.read_csv(m96_kelley_notation_filename, sep='\t')
    m96_kelley = normalize_matrix_to_type_sum(m96_kelley_not_norm)


    ##normalize 192 BOTH by context abundance and sum and roll up in PAPER notation
    m192_norm_filename = normalize_frequencies_192_context(
    path_to_192_matrix=m192_raw_filename)
    m96_kelley_norm_filename = reduce_192_to_96_kelley_notation(m192_norm_filename)
    m96_kelley_norm = pd.read_csv(m96_kelley_norm_filename, sep='\t')

    # ##normalize 192 BOTH by context abundance and sum and roll up in COMMON notation
    # #common notation to plot spectra as in cosmic
    # m192 = pd.read_csv(m192_norm_filename, sep='\t')
    # m96_norm_common_filename = \
    #     reduce_192_to_96_common_notation(path_to_192_matrix=m192_norm_filename)
    # m96 = pd.read_csv(m96_norm_common_filename, sep='\t')

    # plotting spectra
    # plot_mutational_spectra(matrix=m192,
    #                         species_names=['hg38_VS_panTro4.r_ponA'],
    #                         title='192 normalized by context and sum ')
    #
    # species_to_plot = ['hg38_VS_panTro4.r_inner.18',
    #                    'panPan1_VS_hg38.r_inner.18',
    #                    'gorGor3_VS_hg38.r_inner.19',
    #                    'ponAbe2_VS_hg38.r_inner.20',
    #                    'panTro4_VS_panPan1.r_inner.17'
    #                    ]
    # plot_mutational_spectra(matrix=m96,
    #                         species_names=species_to_plot,
    #                         title='96 normalized by context and sum',
    #                         save_to_pdf=False)
    #
    #
    # plot_mutational_spectra(matrix=m96_kelley,
    #                         species_names=species_to_plot,
    #                         title='Without context norm ',
    #                         kelley_notation=True,
    #                         save_to_pdf=False)

    # plotting heatmaps

    ratios_df_3 = create_ratio_table(m96=m96_kelley,
                                   hg38_pair='panTro4',hg38_ref='inner.18',
                                   pantro_pair='hg38',pantro_ref='inner.18',
                                   ponabe_pair='hg38',ponabe_ref='inner.20',
                                   gorgor_pair='hg38',gorgor_ref='inner.19',
                                   panpan_pair='panTro4',panpan_ref='inner.17')
    ratios_df = create_ratio_table(m96=m96_kelley,
                                   hg38_pair='panTro4',hg38_ref='inner.18',
                                   pantro_pair='hg38',pantro_ref='inner.18',
                                   ponabe_pair='hg38',ponabe_ref='inner.20',
                                   gorgor_pair='hg38',gorgor_ref='inner.19',
                                   panpan_pair='hg38',panpan_ref='inner.18')
    ratios_df_norm = create_ratio_table(m96=m96_kelley_norm,
                                   hg38_pair='panTro4',hg38_ref='inner.18',
                                   pantro_pair='hg38',pantro_ref='inner.18',
                                   ponabe_pair='hg38',ponabe_ref='inner.20',
                                   gorgor_pair='hg38',gorgor_ref='inner.19',
                                   panpan_pair='hg38',panpan_ref='inner.18')
    ratios_df_oran = create_ratio_table(m96=m96_kelley,
                                   hg38_pair='panTro4',hg38_ref='ponA',
                                   pantro_pair='hg38',pantro_ref='ponA',
                                   ponabe_pair='gorGor3',ponabe_ref='nomL',
                                   gorgor_pair='panPan1',gorgor_ref='ponA',
                                   panpan_pair='hg38',panpan_ref='ponA')

    ratios_groups = [
        ['hg38_VS_panTro4',
         'hg38_VS_ponAbe2',
         'panTro4_VS_ponAbe2',
         'gorGor3_VS_ponAbe2'],
        ['panTro4_VS_hg38',
         "panPan1_VS_hg38",
         'gorGor3_VS_hg38'],
        ['hg38_VS_gorGor3',
         "panTro4_VS_gorGor3",
         "panPan1_VS_gorGor3"],
        ['hg38_VS_panPan1',
         "panTro4_VS_panPan1",
         "gorGor3_VS_panPan1"],
        ['hg38_VS_panTro4',
         "panPan1_VS_panTro4",
         "gorGor3_VS_panTro4"]

    ]
    for i,ratios_list in enumerate(ratios_groups):
        PREFIX = 'fig_'+str(i)+'.'
        plot_heatmaps_from_ratios_list(ratios_table=ratios_df,
                                       ratios_to_plot=ratios_list,
                                       title=PREFIX+'1\n'+'Ancestral alleles from Great Ape genome',
                                       kelley_notation=True,
                                       if_save=True,
                                       filename=PREFIX+'1.rec_mrca')

        plot_heatmaps_from_ratios_list(ratios_table=ratios_df_oran,
                                       ratios_to_plot=ratios_list,
                                       title=PREFIX+'2\n'+'Orangutan as outgroup ',
                                       kelley_notation=True,
                                       if_save=True,
                                       filename=PREFIX + '2.orang_outgroup'
                                       )

        plot_heatmaps_from_ratios_list(ratios_table=ratios_df_norm,
                                       ratios_to_plot=ratios_list,
                                       title=PREFIX+'3\n'+'Ancestral alleles from Great Ape genome (context normalized)',
                                       kelley_notation=True,
                                       if_save=True,
                                       filename=PREFIX + '3.rec_mrca.cnorm'
                                       )



#outgroup Orang
#ancestral alleles from Great Ape diversity panel
run_process()

# matrix_to_r_output(
#     os.path.join(MATRICES_PATH, '96matrix.recent.normalized.csv'),
#     new_order=True)

