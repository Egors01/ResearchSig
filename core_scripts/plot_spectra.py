import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set()
m96 = pd.read_csv(
    '/home/egors/Research/ExploreSig/'
    'data/variants/matrix_tables/'
    '96matrix.run_999.normalized.csv',sep='\t')
for colname in ['hg38_VS_panTro4.r_ponA','panTro4_VS_hg38.r_ponA']:
    g = sns.FacetGrid(m96, col="SUBS", hue="SUBS", col_wrap=6)
    g.map(sns.barplot, "Context", colname)
    g.set_xticklabels(rotation=90)
    g.set_titles("{col_name} ")
    plt.show()


