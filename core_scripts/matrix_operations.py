from tools.logger import Logger
from tools.names_generator import get_sepecies_names
from tools.context_count import load_context_abundance
import pandas as pd
import constants
import os
import numpy as np

def normalize_96_to_total_sum(path_to_96_matrix):
    all_context_sum = {}
    m96 = pd.read_csv(path_to_96_matrix, sep='\t')
    m96_norm_kelley = m96[['SUBS', 'Context']].copy()
    for colname in m96.columns.values:
        if colname != 'SUBS' and colname != 'Context' and colname != 'Full_Context':
            m96_norm_kelley[colname] = m96[colname] / m96[colname].sum()
            all_context_sum[colname] = m96[colname].sum()
    new_matrix_name = path_to_96_matrix.split('.csv')[0] + '.normsum.csv'
    print(new_matrix_name)
    m96_norm_kelley.to_csv(new_matrix_name, sep='\t', index=False)

    return new_matrix_name



def normalize_frequencies_192_context(path_to_192_matrix, run_id='999',
                                      just_name=False):
    if just_name:
        return path_to_192_matrix.split('.csv')[
                   0] + '.192norm.csv'

    print('normalizing 192 by context counts ')
    # names = get_sepecies_names()
    context_abundance = load_context_abundance()

    m192 = pd.read_csv(path_to_192_matrix, sep='\t')
    # m96 = m96.rename(columns=lambda x: x.split('_')[0])
    colname_to_sp_name = dict(zip(
        [x for x in m192.columns.values if x != 'SUBS' and x != 'Context'],
        [x.split('_')[0] for x in m192.columns.values if
         x != 'SUBS' and x != 'Context']))

    m192['Full_context'] = m192['Context'].str[0] + m192['SUBS'].str[0] + \
                           m192['Context'].str[-1]

    for colname in colname_to_sp_name.keys():
        name = colname_to_sp_name[colname]
        for idx in m192.index:
            full_context = m192.ix[idx, 'Full_context']
            denominator = \
                context_abundance.loc[
                    context_abundance["Context"] == full_context][
                    name].item()
            m192.ix[idx, colname] = m192.ix[idx, colname] / denominator

    new_matrix_name = path_to_192_matrix.split('.csv')[0] + '.normcontext.csv'
    m192.to_csv(new_matrix_name, sep='\t', index=False)

    return new_matrix_name

def get_complementary(nucleotide):

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return "C"
    elif nucleotide == 'C':
        return "G"
    else:
        return '.'


def reduce_192_to_96_common_notation(path_to_192_matrix, run_id='999'):
    print('started matrix rollup. File_path ={} \n'.format(path_to_192_matrix))
    m192 = pd.read_csv(path_to_192_matrix, sep='\t')

    colname_to_sp_name = dict(zip(
        [x for x in m192.columns.values if x != 'SUBS' and x != 'Context'],
        [x.split('_')[0] for x in m192.columns.values if
         x != 'SUBS' and x != 'Context' and x != 'Full_context']))

    for colname in colname_to_sp_name.keys():
        for idx in m192.index:
            sbs = m192.ix[idx, 'SUBS']
            if sbs == 'A>T':
                m192.ix[idx, 'SUBS'] = 'T>A'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'G>C':
                m192.ix[idx, 'SUBS'] = 'C>G'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'G>A':
                m192.ix[idx, 'SUBS'] = 'C>T'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'A>G':
                m192.ix[idx, 'SUBS'] = 'T>C'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'A>C':
                m192.ix[idx, 'SUBS'] = 'T>G'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'G>T':
                m192.ix[idx, 'SUBS'] = 'C>A'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])

    m96 = m192.groupby(['SUBS', 'Context'], as_index=False).sum()

    new_matrix_name = os.path.join(os.path.dirname(path_to_192_matrix),
                                   '96' + path_to_192_matrix.split('192_')[-1])
    m96.to_csv(new_matrix_name, sep='\t', index=False)

    return new_matrix_name

def reduce_192_to_96_kelley_notation(path_to_192_matrix, run_id='999'):
    print('started matrix rollup. File_path ={} \n'.format(path_to_192_matrix))
    m192 = pd.read_csv(path_to_192_matrix, sep='\t')

    colname_to_sp_name = dict(zip(
        [x for x in m192.columns.values if x != 'SUBS' and x != 'Context'],
        [x.split('_')[0] for x in m192.columns.values if
         x != 'SUBS' and x != 'Context' and x != 'Full_context']))

    for colname in colname_to_sp_name.keys():
        for idx in m192.index:
            sbs = m192.ix[idx, 'SUBS']
            if sbs == 'T>A':
                m192.ix[idx, 'SUBS'] = 'A>T'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'G>C':
                m192.ix[idx, 'SUBS'] = 'C>G'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'G>A':
                m192.ix[idx, 'SUBS'] = 'C>T'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'T>C':
                m192.ix[idx, 'SUBS'] = 'A>G'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'T>G':
                m192.ix[idx, 'SUBS'] = 'A>C'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])
            elif sbs == 'G>T':
                m192.ix[idx, 'SUBS'] = 'C>A'
                m192.ix[idx, 'Context'] = "".join(
                    [get_complementary(x) for x in m192.ix[idx, 'Context']])

    m96 = m192.groupby(['SUBS', 'Context'], as_index=False).sum()

    new_matrix_name = os.path.join(os.path.dirname(path_to_192_matrix),
                                   '96' + path_to_192_matrix.split('192_')[-1])

    new_matrix_name =  new_matrix_name.split('.csv')[0] + '.kel.csv'

    m96.to_csv(new_matrix_name, sep='\t', index=False)

    return new_matrix_name


def matrix_to_r_output(path_to_96_matrix,new_order=False, run_id='999'):
    m96 = pd.read_csv(path_to_96_matrix, sep='\t')
    m96['subs'] = m96.apply(
        lambda row: row.Context[0] + '[' + row.SUBS + ']' + row.Context[2],
        axis=1)
    m96['ind'] = 0
    if new_order:
        sorter = constants.R_ORDER
        sorterIndex = dict(zip(sorter, range(len(sorter))))
        m96['ind'] = m96['subs'].map(sorterIndex)
        m96.sort_values(by=['ind'], ascending=[True], inplace=True)
        new_matrix_name = path_to_96_matrix.split('.csv')[0] + '.R.newind.csv'
    else:
        new_matrix_name = path_to_96_matrix.split('.csv')[0] + '.R.csv'

    m96.drop(['ind', 'SUBS', 'Context'], 1, inplace=True)
    m96.to_csv(new_matrix_name, sep='\t', index=False)

    return new_matrix_name

def create_ratio_table(path_to_matrix_file):
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


def create_heatmap_df_for_species_common(relations_df, column_name):
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

def create_heatmap_df_for_species_kelley(relations_df, column_name):
    relations = relations_df

    iterables = [['C>A', 'C>G', 'C>T', 'A>G', 'A>C','A>T'],
                 ["T", "G", "C", "A"]]
    index = pd.MultiIndex.from_product(iterables, names=['subs', "5'- "])
    template_heatmap_df = pd.DataFrame(np.zeros(shape=(24, 4)),
                                       columns={"A", "T", "G", "C"},
                                       index=index)
    template_heatmap_df = template_heatmap_df[["A", "C", "G", "T"]]

    for idx in relations.index:
        subs = relations.loc[idx, 'SUBS']
        three_prime = relations.loc[idx, 'Context'][2]
        five_prime = relations.loc[idx, 'Context'][0]
        value = relations.loc[idx, column_name]
        template_heatmap_df.loc[subs,five_prime][three_prime] = value

    return template_heatmap_df

# def normalize_frequencies_96_context(path_to_96_matrix, run_id='999', just_name=False):
#     if just_name:
#         return path_to_96_matrix.split('.csv')[0] + '.normalized.csv'
#     logger = Logger(source_name='matrix_normilize',
#                     msg='started matrix normilization', run_id=run_id)
#     # names = get_sepecies_names()
#     sp_names = set()
#     context_abundance = load_context_abundance()
#
#     m96 = pd.read_csv(path_to_96_matrix, sep='\t')
#     # m96 = m96.rename(columns=lambda x: x.split('_')[0])
#     for s in m96.columns.values[:-1]:
#         sp_names.add(s.split('_')[0])
#     sp_names = list(sp_names)
#     for name in sp_names:
#         for context in context_abundance["Context"]:
#             mult = \
#                 context_abundance.loc[context_abundance["Context"] == context][
#                     name].item()
#             # print(" sp_name {} context {} mult {}".format(name,context,mult))
#             m96.update(m96.loc[m96["SUBS"].str.contains(
#                 r'' + context[0] + '.....' + context[2]),
#                                m96.filter(regex=r'' + name + '_.*',
#                                           axis=1).columns] / mult)
#     new_matrix_name = path_to_96_matrix.split('.csv')[0] + '.normalized.csv'
#     m96.to_csv(new_matrix_name, sep='\t', index=False)
#     logger.print_end()
#
#     return new_matrix_name