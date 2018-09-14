import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt





# fill heatmap for pairname to plot


def create_relations_table_kelley(path_to_matrix_file):
    m96 = pd.read_csv(path_to_matrix_file, sep='\t')
    relations_df = m96[['SUBS', 'Context']].copy()
    # create df with single column for comparation
    for colname in m96.columns.values:
        if colname != 'SUBS' and colname != 'Context':
            name = colname.split('_')[0]
            pair = colname.split('_')[2].split('.')[0]
            ref = colname.split('.r_')[-1]
            new_col_name = colname
            pair_col_name = pair + '_VS_' + name + '.r_' + ref
            if colname not in relations_df.columns.values:
                relations_df[colname] = m96[colname] / m96[pair_col_name]

    return relations_df


def create_heatmap_df_for_species(relations_df, column_name):
    relations = relations_df

    iterables = [['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'],
                 ["A", "T", "G", "C"]]
    index = pd.MultiIndex.from_product(iterables, names=['subs', "3'- "])
    template_heatmap_df = pd.DataFrame(np.zeros(shape=(24, 4)),
                                       columns={"A", "C", "G", "T"},
                                       index=index)
    template_heatmap_df = template_heatmap_df[["A", "C", "G", "T"]]

    for idx in relations.index:
        subs = relations.loc[idx, 'SUBS']
        three_prime = relations.loc[idx, 'Context'][0]
        five_prime = relations.loc[idx, 'Context'][2]
        value = relations.loc[idx, column_name]
        template_heatmap_df.loc[subs, three_prime][five_prime] = value

    return template_heatmap_df


normalize_96_as_kelley(path_to_96_matrix='96matrix.run_0_without_norm.csv')
relations_df = create_relations_table_kelley(path_to_matrix_file=
                                             '96matrix.run_0_without_norm.normkelley.csv')

relations_to_plot = ['hg38_VS_panTro4.r_ponA', 'panTro4_VS_hg38.r_ponA']
heatmaps_to_plot = [
    create_heatmap_df_for_species(relations_df=relations_df, column_name=x)
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
