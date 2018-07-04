import os
import re
import pandas as pd
from tools.set_env import Environment

def load_ref_startpos():
    startpos_dict={}
    CHR_NAMES=['chr1', 'chr2', 'chr3', 'chr4',
           'chr5','chr6','chr7','chr8',
           'chr9','chr10','chr11','chr12',
           'chr13','chr14','chr15','chr16',
           'chr17','chr18','chr19','chr20',
           'chr21','chr22','chrX','chrY','example','example2','example3','example4']
    for chrname in CHR_NAMES:
        file_path=os.path.join(os.getcwd(),'maf_templates',chrname+'.maf')
        inp = open(file_path, "r")
        for line in inp:
            if "hg38" in line and line[0] == 's':
                startpos_dict[chrname]=int(re.findall(r'\w+', line)[3])
                break
        inp.close()
    startpos_df =pd.DataFrame(startpos_dict,index=['pos'])
    startpos_df.to_csv(Environment().get_ref_startpospath(), sep='\t', index=False)
    return startpos_df

def get_ref_startpos(chrname):
    startpos_df = pd.read_csv(Environment().get_ref_startpospath(),sep='\t')
    startpos = startpos_df[chrname].values[0]
    return startpos
