from tools.logger import Logger
from tools.names_generator import get_sepecies_names
from tools.context_count import load_context_abundance
import pandas as pd
import constants
import os


def normalize_96_as_kelley(path_to_96_matrix):
    all_context_sum = {}
    m96 = pd.read_csv(path_to_96_matrix, sep='\t')
    m96_norm_kelley = m96[['SUBS', 'Context']].copy()
    for colname in m96.columns.values:
        if colname != 'SUBS' and colname != 'Context':
            m96_norm_kelley[colname] = m96[colname] / m96[colname].sum()
            all_context_sum[colname] = m96[colname].sum()
    new_matrix_name = path_to_96_matrix.split('.csv')[0] + '.normkelley.csv'
    print(new_matrix_name)
    m96_norm_kelley.to_csv(new_matrix_name, sep='\t', index=False)
    return new_matrix_name


def normalize_frequencies_96(path_to_96_matrix, run_id='999', just_name=False):
    if just_name:
        return path_to_96_matrix.split('.csv')[0] + '.normalized.csv'
    logger = Logger(source_name='matrix_normilize',
                    msg='started matrix normilization', run_id=run_id)
    # names = get_sepecies_names()
    sp_names = set()
    context_abundance = load_context_abundance()

    m96 = pd.read_csv(path_to_96_matrix, sep='\t')
    # m96 = m96.rename(columns=lambda x: x.split('_')[0])
    for s in m96.columns.values[:-1]:
        sp_names.add(s.split('_')[0])
    sp_names = list(sp_names)
    for name in sp_names:
        for context in context_abundance["Context"]:
            mult = \
                context_abundance.loc[context_abundance["Context"] == context][
                    name].item()
            # print(" sp_name {} context {} mult {}".format(name,context,mult))
            m96.update(m96.loc[m96["SUBS"].str.contains(
                r'' + context[0] + '.....' + context[2]),
                               m96.filter(regex=r'' + name + '_.*',
                                          axis=1).columns] / mult)
    new_matrix_name = path_to_96_matrix.split('.csv')[0] + '.normalized.csv'
    m96.to_csv(new_matrix_name, sep='\t', index=False)
    logger.print_end()

    return new_matrix_name


def normalize_frequencies_192(path_to_192_matrix, run_id='999', just_name=False):
    if just_name:
        return path_to_192_matrix.split('.csv')[
                   0] + '.normalized.csv'

    logger = Logger(source_name='matrix_normilize',
                    msg='started matrix normilization', run_id=run_id)
    # names = get_sepecies_names()
    context_abundance = load_context_abundance()

    m192 = pd.read_csv(path_to_192_matrix, sep='\t')
    # m96 = m96.rename(columns=lambda x: x.split('_')[0])
    colname_to_sp_name = dict(zip(
        [x for x in m192.columns.values if x != 'SUBS' and x != 'Context'],
        [x.split('_')[0] for x in m192.columns.values if
         x != 'SUBS' and x != 'Context']))
    sp_names = list(set([x.split('_')[0] for x in m192.columns.values if
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

    new_matrix_name = path_to_192_matrix.split('.csv')[0] + '.normalized.csv'
    m192.to_csv(new_matrix_name, sep='\t', index=False)
    logger.print_end()

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


def reduce_by_strand(path_to_192_matrix, run_id='999'):
    logger = Logger(source_name='matrix rollup ',
                    msg='started', run_id=run_id)
    m192 = pd.read_csv(path_to_192_matrix, sep='\t')
    # m96 = m96.rename(columns=lambda x: x.split('_')[0])
    # dict to select column by species name
    colname_to_sp_name = dict(zip(
        [x for x in m192.columns.values if x != 'SUBS' and x != 'Context'],
        [x.split('_')[0] for x in m192.columns.values if
         x != 'SUBS' and x != 'Context' and x != 'Full_context']))
    # sp_names = list(set([x.split('_')[0] for x in m192.columns.values if
    #                      x != 'SUBS' and x != 'Context']))
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
    pass
    m96 = m192.groupby(['SUBS', 'Context'], as_index=False).sum()
    new_matrix_name = os.path.join(os.path.dirname(path_to_192_matrix),
                                   '96' + path_to_192_matrix.split('192_')[-1])
    m96.to_csv(new_matrix_name, sep='\t', index=False)
    logger.print_end()

    return new_matrix_name


def matrix_to_r_output(path_to_96_matrix, run_id='999'):
    m96 = pd.read_csv(path_to_96_matrix, sep='\t')
    m96['subs'] = m96.apply(
        lambda row: row.Context[0] + '[' + row.SUBS + ']' + row.Context[2],
        axis=1)
    sorter = constants.R_ORDER
    sorterIndex = dict(zip(sorter, range(len(sorter))))
    m96['ind'] = m96['subs'].map(sorterIndex)
    m96.sort_values(by=['ind'], \
                    ascending=[True], inplace=True)
    m96.drop(['ind', 'SUBS', 'Context'], 1, inplace=True)
    # del m96['SUBS']
    # del m96['Context']
    new_matrix_name = path_to_96_matrix.split('.csv')[0] + '.R.csv'

    # sorterIndex = dict(zip(sorter, range(len(sorter))))

    # Generate a rank column that will be used to sort
    # the dataframe numerically
    # df['Tm_Rank'] = df['Tm'].map(sorterIndex)

    m96.to_csv(new_matrix_name, sep='\t', index=False)

    return new_matrix_name
