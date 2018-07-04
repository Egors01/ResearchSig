import os
import re
import pandas as pd
from tools.set_env import Environment
from tools.logger import Logger

def get_names_from_maf(chrname):
    maf_path=os.path.join(Environment().raw_maf_path,chrname+'.maf')
    inp = open(maf_path, "r")
    names = []
    linecount=0
    for line in inp:
        linecount+=1
        if line[0] == 's':
            newname=re.search(r'\w+[.]\w+',line).group()
            newname=newname[:newname.index('.')]
            if newname not in names:
                names.append(newname)
        if linecount>1000000:
            break
    inp.close()
    return names

def make_chrs_names_file():
    names_dict = {}
    CHR_NAMES = ['chr1','chrY','example', 'example2', 'example3', 'example4']
    for chmane in CHR_NAMES:
        names_chr_list = get_names(chmane)
        while len(names_chr_list) < 20:
            names_chr_list.append('')
        names_dict[chmane] = names_chr_list
        print(chmane,names_chr_list)
    names_df = pd.DataFrame(names_dict)
    names_df.to_csv(Environment().species_names_path,sep='\t',index=False)
    return names_df

def get_names(chrname):
    try:
        names_df =pd.read_csv(Environment().species_names_path,sep='\t')
        sp_names = names_df[chrname].dropna().values
    except:
        log_def = Logger(msg = 'cannot load sp_names from file. trying from maf...')
        try:
            sp_names = get_names_from_maf(chrname)
            log_def.message(msg='loaded from maf ')
        except FileNotFoundError:
            log_def.message(msg ='unable to load from maf')
            raise
    return sp_names