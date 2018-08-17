import os
import re

import pandas as pd

from constants import OUTPUT_FILT_VCF_NAMES_PATH, FILT_VCF_VARIANT_PATH, VCF_VARIANT_PATH, SPECIES_NAMES_PATH
from tools.logger import Logger
from tools.set_env import Environment


def create_vcf_name_from_pair_and_ref(sp_name, sp_pair, refname, filt=False):
    if not filt:
        return os.path.join(VCF_VARIANT_PATH, sp_name + '_VS_' + sp_pair + '.r_' + refname[:4] + '.vcf')

def make_vcf_filt_name(sp_vcf_old_name):
     new_name = os.path.join(FILT_VCF_VARIANT_PATH,
                sp_vcf_old_name.split('/')[-1].split('.vcf')[0] + '.filt.vcf')
     return new_name
def create_output_file_names(PAIRS_AND_REF):
    vcf_variants_filenames_list = []
    filt_vcf_variants_filenames_list = []
    filt = open(OUTPUT_FILT_VCF_NAMES_PATH, 'a')
    for pair_and_ref in PAIRS_AND_REF:
        # TODO change to combination
        for sp, sp_pair in zip([0, 1], [1, 0]):
            vcf_name = create_vcf_name_from_pair_and_ref(sp_name=pair_and_ref[sp],
                                                         sp_pair=pair_and_ref[sp_pair],
                                                         refname=pair_and_ref[2])
            vcf_variants_filenames_list.append(vcf_name)

    return vcf_variants_filenames_list


def get_names_from_maf(chrname):
    maf_path = os.path.join(Environment().raw_maf_path, chrname + '.maf')
    inp = open(maf_path, "r")
    names = []
    linecount = 0
    for line in inp:
        linecount += 1
        if line[0] == 's':
            newname = re.search(r'\w+[.]\w+', line).group()
            newname = newname[:newname.index('.')]
            if newname not in names:
                names.append(newname)
        if linecount > 1000000:
            break
    inp.close()
    return names


def make_chrs_names_file():
    names_dict = {}
    CHR_NAMES = ['chr1', 'chrY', 'example', 'example2', 'example3', 'example4']
    for chmane in CHR_NAMES:
        names_chr_list = get_names(chmane)
        while len(names_chr_list) < 20:
            names_chr_list.append('')
        names_dict[chmane] = names_chr_list
        print(chmane, names_chr_list)
    names_df = pd.DataFrame(names_dict)
    names_df.to_csv(Environment().species_names_path, sep='\t', index=False)
    return names_df

def get_sepecies_names():
    names_df=pd.read_csv(SPECIES_NAMES_PATH, sep='\t')
    sp_names = pd.read_csv(SPECIES_NAMES_PATH, header=None)[0].tolist()
    return sp_names

get_sepecies_names()

def get_names(chrname='chrY'):
    try:
        names_df = pd.read_csv(SPECIES_NAMES_PATH, sep='\t')
        sp_names = names_df[chrname].dropna().values
    except:
        log_def = Logger(msg='cannot load sp_names from file. trying from maf...')
        try:
            sp_names = get_names_from_maf(chrname)
            log_def.message(msg='loaded from maf ')
        except FileNotFoundError:
            log_def.message(msg='unable to load from maf')
            raise
    return sp_names
