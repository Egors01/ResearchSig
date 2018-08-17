from constants import MATRIX_NAME
from tools.names_generator import get_sepecies_names
from tools.context_count import load_context_abundance
import pandas as pd

def normilize_frequencies(path_to_96_matrix=MATRIX_NAME):

    #names = get_sepecies_names()
    sp_names=set()
    context_abundance = load_context_abundance()

    m96 = pd.read_csv(path_to_96_matrix, sep='\t')
    #m96 = m96.rename(columns=lambda x: x.split('_')[0])
    for s in m96.columns.values[:-1]:
        sp_names.add(s.split('_')[0])
    sp_names = list(sp_names)
    print(s)
    # print(sp_names)
    # for name in sp_names:
    #     col_to_select = []
    #     a=[colname for colname in m96.columns.values if name + '_VS_' in colname ]
    #     m96.set_index(m96['SUBS'],inplace=True)
    #     for colname in m96.columns.values:
    #         if name + '_VS_' in colname:
    #             col_to_select.append(colname)
    #         for context in context_abundance["Context"].values:
    #             norm_factor =
    #             #context_abundance.loc[context_abundance["Context"] == context][name].item()
    #             m96.loc[m96["SUBS"].str.contains(r'' + context[0] + '.....' + context[2]), col_to_select] =\
    #                 m96.loc[m96["SUBS"].str.contains(r'' + context[0] + '.....' + context[2]), col_to_select] / norm_factor
    for name in sp_names:
        for context in context_abundance["Context"]:
            mult=context_abundance.loc[context_abundance["Context"]==context][name].item()
            m96.update(m96.loc[m96["SUBS"].str.contains(r''+context[0]+'.....'+context[2]),
                    m96.filter(regex=r''+name+'_.*', axis=1).columns]/mult)
    new_matrix_name=path_to_96_matrix.split('.csv')[0]+'.normalized.csv'
    print(new_matrix_name)
    m96.to_csv(new_matrix_name,sep='\t',index=False)
