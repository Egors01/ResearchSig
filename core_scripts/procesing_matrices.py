import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt

from constants import MATRICES_PATH
from core_scripts.matrix_operations import normalize_96_to_total_sum, \
    reduce_192_to_96_common_notation, reduce_192_to_96_kelley_notation, \
    create_ratio_table, create_heatmap_df_for_species_kelley, \
    normalize_frequencies_192_context, matrix_to_r_output


def plot_heatmaps_from_ratios_list(ratios_table, ratios_to_plot,title=''):
    heatmaps_to_plot = [
        create_heatmap_df_for_species_kelley(relations_df=ratios_table,
                                             column_name=x)
        for x in ratios_to_plot]

    fig = plt.figure(figsize=(20, 10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i, heatmap in enumerate(heatmaps_to_plot):
        ax = fig.add_subplot(1, len(heatmaps_to_plot), i + 1)
        # ax.text(0.5, 0.5, str((1, 4, i)), fontsize=18, ha='center')
        sns.heatmap(heatmap, annot=True, annot_kws={"size": 9}, linewidths=.6,
                    center=1, cmap='RdBu_r', ax=ax)
        fig.subplots_adjust(hspace=0.8)
        ax.set_title(ratios_to_plot[i])
    fig.suptitle(title, size=16)
    fig.subplots_adjust(top=.5)
    plt.show()


def plot_mutational_spectra(matrix, species_name, title='',
                            kelley_notation=False):
    if not kelley_notation:
        if len(matrix['SUBS'].values[:]) > 96:
            g = sns.FacetGrid(matrix, col="SUBS", hue="SUBS", col_wrap=6,
                              col_order=['C>A', 'C>G', 'C>T', 'T>A', "T>C",
                                         "T>G",
                                         'G>T', 'G>T', 'G>A', 'A>T', "A>G",
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
    g.fig.suptitle(title, size=16)
    g.fig.subplots_adjust(top=.5)

    plt.show()


def with_context_norm():
    m192_raw_filename = os.path.join(MATRICES_PATH, '192_matrix.recent.csv')

    # normalize 192 by context counts and contribution sum per species
    m192_norm_filename = normalize_frequencies_192_context(
        path_to_192_matrix=m192_raw_filename)
    m192_norm_sum_filename = normalize_96_to_total_sum(
        path_to_96_matrix=m192_norm_filename)

    # rollup to 96 and norm to contribution
    m96_common_filename = \
        reduce_192_to_96_common_notation(path_to_192_matrix=m192_norm_filename)
    m96_common_normsum_filename = \
        normalize_96_to_total_sum(path_to_96_matrix=m96_common_filename)

    # just roll up and normalize by sum as in paper
    m96_filename_normc_kelley_notation = \
        reduce_192_to_96_kelley_notation(path_to_192_matrix=m192_norm_filename)
    m96_filename_normcs_kelley = \
        normalize_96_to_total_sum(
            path_to_96_matrix=m96_filename_normc_kelley_notation)

    # plotting spectra
    m96 = pd.read_csv(m96_common_normsum_filename, sep='\t')
    m96_norm_kelley = pd.read_csv(m96_filename_normcs_kelley, sep='\t')
    m192 = pd.read_csv(m192_norm_sum_filename, sep='\t')

    plot_mutational_spectra(matrix=m192,
                            species_name='hg38_VS_panTro4.r_ponA',
                            title='192 normalized by context and sum ')

    plot_mutational_spectra(matrix=m96,
                            species_name='hg38_VS_panTro4.r_ponA',
                            title='192 normalized by context and sum, reduced to 96')

    plot_mutational_spectra(matrix=m96_norm_kelley,
                            species_name='hg38_VS_panTro4.r_ponA',
                            title='192 normalized reduced to 96 in kelley not ',
                            kelley_notation=True)

    # plotting heatmaps
    ratios_df = create_ratio_table(
        path_to_matrix_file= m96_filename_normcs_kelley)
    ratios_list = ['hg38_VS_panTro4.r_ponA', "hg38_VS_ponAbe2.r_nomL",
                   "gorGor3_VS_ponAbe2.r_nomL"]
    plot_heatmaps_from_ratios_list(ratios_table=ratios_df,
                                   ratios_to_plot=ratios_list)


with_context_norm()
# matrix_to_r_output(
#     os.path.join(MATRICES_PATH, '96matrix.recent.normalized.csv'),
#     new_order=False)
# matrix_to_r_output(
#     os.path.join(MATRICES_PATH, '96matrix.recent.normalized.csv'),
#     new_order=True)


def kelley_norm():
    m96_filename = \
        reduce_192_to_96_kelley_notation(path_to_192_matrix=
                                         os.path.join(MATRICES_PATH,
                                                      '192_matrix.recent.csv'))
    norm_matrix_kelly_filename = \
        normalize_96_to_total_sum(path_to_96_matrix=m96_filename)
    ratios_df = create_ratio_table(path_to_matrix_file=
                                   norm_matrix_kelly_filename)

    # plotting spectra
    m96 = pd.read_csv(norm_matrix_kelly_filename, sep='\t')
    spectra_to_plot = ['hg38_VS_panTro4.r_ponA', 'panTro4_VS_hg38.r_ponA']
    for colname in spectra_to_plot:
        g = sns.FacetGrid(m96, col="SUBS", hue="SUBS", col_wrap=6)
        g.map(sns.barplot, "Context", colname)
        g.set_xticklabels(rotation=90)
        g.set_titles("{col_name} ")
        plt.show()

    # plotting heatmaps
    relations_to_plot = ['hg38_VS_panTro4.r_ponA', 'panTro4_VS_hg38.r_ponA']
    heatmaps_to_plot = [
        create_heatmap_df_for_species_kelley(relations_df=ratios_df,
                                             column_name=x)
        for x in relations_to_plot]

    fig = plt.figure(figsize=(20, 10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i, heatmap in enumerate(heatmaps_to_plot):
        ax = fig.add_subplot(1, len(heatmaps_to_plot), i + 1)
        # ax.text(0.5, 0.5, str((1, 4, i)), fontsize=18, ha='center')
        sns.heatmap(heatmap, annot=True, annot_kws={"size": 9}, linewidths=.6,
                    center=1, cmap='RdBu_r', ax=ax)
        fig.subplots_adjust(hspace=0.8)
        ax.set_title(relations_to_plot[i])
    plt.show()
