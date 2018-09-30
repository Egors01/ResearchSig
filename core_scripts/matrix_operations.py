from tools.logger import Logger
from tools.names_generator import get_sepecies_names
from tools.context_count import load_context_abundance
import pandas as pd
import constants
import os
import numpy as np

def normalize_matrix_to_type_sum(matrix):
    all_context_sum = {}
    matrix_norm_sum = matrix[['SUBS', 'Context']].copy()
    for colname in matrix.columns.values:
        if colname != 'SUBS' and colname != 'Context' and colname != 'Full_context':
            matrix_norm_sum[colname] = matrix[colname] / matrix[colname].sum()
            all_context_sum[colname] = matrix[colname].sum()

    return matrix_norm_sum


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

    m192_norm_sum = normalize_matrix_to_type_sum(m192)
    new_matrix_name = path_to_192_matrix.split('.csv')[0] + '.norm_sum_cont.csv'
    m192_norm_sum.to_csv(new_matrix_name, sep='\t', index=False)

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
def reverse_context(context):
    return "".join([get_complementary(x) for x in ''.join(reversed(context)) ])

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
            context = m192.ix[idx, 'Context']
            if sbs == 'A>T':
                m192.ix[idx, 'SUBS'] = 'T>A'
                m192.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>C':
                m192.ix[idx, 'SUBS'] = 'C>G'
                m192.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>A':
                m192.ix[idx, 'SUBS'] = 'C>T'
                m192.ix[idx, 'Context'] =  reverse_context(context)
            elif sbs == 'A>G':
                m192.ix[idx, 'SUBS'] = 'T>C'
                m192.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'A>C':
                m192.ix[idx, 'SUBS'] = 'T>G'
                m192.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>T':
                m192.ix[idx, 'SUBS'] = 'C>A'
                m192.ix[idx, 'Context'] = reverse_context(context)

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
            context = m192.ix[idx, 'Context']
            if sbs == 'T>A':
                m192.ix[idx, 'SUBS'] = 'A>T'
                m192.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>C':
                m192.ix[idx, 'SUBS'] = 'C>G'
                m192.ix[idx, 'Context'] =reverse_context(context)
            elif sbs == 'G>A':
                m192.ix[idx, 'SUBS'] = 'C>T'
                m192.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'T>C':
                m192.ix[idx, 'SUBS'] = 'A>G'
                m192.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'T>G':
                m192.ix[idx, 'SUBS'] = 'A>C'
                m192.ix[idx, 'Context'] = reverse_context(context)
            elif sbs == 'G>T':
                m192.ix[idx, 'SUBS'] = 'C>A'
                m192.ix[idx, 'Context'] = reverse_context(context)

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

def create_ratio_table(m96,
                       hg38_pair,hg38_ref,
                       pantro_pair,pantro_ref,
                       ponabe_pair,ponabe_ref,
                       gorgor_pair,gorgor_ref,
                       panpan_pair,panpan_ref):
    relations_df = m96[['SUBS', 'Context']].copy()
    # create df with single column for comparation
    pairs_to_compare = [['panTro4','hg38'],
                        ['hg38','ponAbe2'],
                        ['panTro4','ponAbe2'],
                        ['panPan1','ponAbe2'],
                        ['gorGor3','ponAbe2'],
                        ['hg38','gorGor3'],
                        ["panPan1","gorGor3"],
                        ['panTro4', 'gorGor3'],
                        ['panPan1', 'panTro4'],
                        ['panPan1', 'hg38']
                        ]
    name_to_col ={
        'hg38':'hg38'+'_VS_'+hg38_pair+'.r_'+hg38_ref,
        'panTro4': 'panTro4' + '_VS_' + pantro_pair+'.r_'+pantro_ref,
        'gorGor3': 'gorGor3' + '_VS_' + gorgor_pair+'.r_'+gorgor_ref,
        'panPan1': 'panPan1' + '_VS_' + panpan_pair + '.r_' + panpan_ref,
        'ponAbe2': 'ponAbe2' + '_VS_' + ponabe_pair + '.r_' + ponabe_ref,
    }
    # colnames_of_comp = [x for x in m96.columns.values if
    #  'VS_hg38' in x and 'nomLeu3' not in x] + ['hg38_VS_panTro4.r_ponA']
    # species_in_colnames = [x.split('_')[0] for x in colnames_of_comp]
    # name_to_col =dict(zip(species_in_colnames,colnames_of_comp))
    for i,pair in enumerate(pairs_to_compare):
            first = pair[0]
            second = pair[1]
            colname_first = name_to_col[first]
            colname_second = name_to_col[second]
            relations_df[first+'_VS_'+second] = m96[colname_first] / m96[colname_second]
            relations_df[second+'_VS_'+first] = m96[colname_second] / m96[
                colname_first]
    return relations_df


# def create_heatmap_df_for_species_common(relations_df, column_name):
#     relations = relations_df
#
#     iterables = [['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'],
#                  ["A", "T", "G", "C"]]
#     index = pd.MultiIndex.from_product(iterables, names=['subs', "3'- "])
#     template_heatmap_df = pd.DataFrame(np.zeros(shape=(24, 4)),
#                                        columns={"A", "C", "G", "T"},
#                                        index=index)
#     template_heatmap_df = template_heatmap_df[["A", "C", "G", "T"]]
#
#     for idx in relations.index:
#         subs = relations.loc[idx, 'SUBS']
#         five_prime = relations.loc[idx, 'Context'][0]
#         three_prime = relations.loc[idx, 'Context'][2]
#
#         value = relations.loc[idx, column_name]
#         template_heatmap_df.loc[subs, five_prime][three_prime] = value
#
#     return template_heatmap_df

def create_heatmap_df_for_species(relations_df, column_name,kelley_notation=True):
    relations = relations_df
    if kelley_notation:
        iterables = [['C>A', 'C>G', 'C>T', 'A>G', 'A>C','A>T'],
                     ["T", "G", "C", "A"]]
        index = pd.MultiIndex.from_product(iterables, names=['subs', "5'- "])

    else:
        iterables = [['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'],
                     ["A", "T", "G", "C"]]
        index = pd.MultiIndex.from_product(iterables, names=['subs', "5'- "])

    template_heatmap_df = pd.DataFrame(np.zeros(shape=(24, 4)),
                                       columns={"A", "C", "G", "T"},
                                       index=index)
    template_heatmap_df = template_heatmap_df[["A", "C", "G", "T"]]
    for idx in relations.index:
        subs = relations.loc[idx, 'SUBS']
        five_prime = relations.loc[idx, 'Context'][0]
        three_prime = relations.loc[idx, 'Context'][2]
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