import os
import re
import sys
import time
from itertools import product

import pandas as pd

sys.path.append((os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
sys.path.append(os.path.dirname(os.path.realpath(__file__)))


#print(sys.path)
from constants import FASTA_CHRS_PATH, PRIMARY_PATH
from constants import SPECIES_NAMES, CHR_NAMES, RESOURSES_PATH
from logger import Logger
def context_count():
    chr_names = CHR_NAMES[:-2]
    #chr_names = ['chrY']
    sp_names = [ 'chlSab2', 'macFas5',
                 'papAnu2','rheMac3', 'tarSyr2','rhiRox1',]
    #sp_names = ['ponAbe2','hg38']
    contexts = list(product(['A', 'T', 'G', 'C'], repeat=3))
    contexts = [''.join(x) for x in contexts]
    context_abundance = pd.DataFrame()
    context_abundance['Context'] = contexts
    log=Logger(run_id=778, source_name='cc ', filename='context_count', msg='',
                 reset=True)
    for name in sp_names:
        print(name)
        log.tick(timer_name='species')
        log.message(msg='started {}'.format(name))
        context_counts = dict.fromkeys(contexts, 0)
        for CHR in chr_names:# ['chr1']:
            log.message(msg='chr {} started '.format(CHR))
            file_path = os.path.join(FASTA_CHRS_PATH, CHR, name + '.fasta')
            f = open(file_path, 'r')
            f.readline()
            first = f.read(1).upper()
            second = f.read(1).upper()
            third = f.read(1).upper()
            log.tick(timer_name='chr')
            while True:
                if bool(re.search(r'^[ATGC]+$', first + second + third)):
                    context_counts[first + second + third] += 1
                first = second
                second = third
                third = f.read(1).upper()
                if not third:
                    break
            log.message(msg='chr {} finished '.format(CHR),print_time=True,
                        timer_name='chr')
            f.close()
        log.message(msg='finished {} '.format(name),print_time=True,timer_name='species')
        context_abundance[name] = context_abundance['Context'].apply(
            lambda x: context_counts[x])
    context_abundance.to_csv(
        os.path.join(PRIMARY_PATH, 'context_abundance.csv'),
        sep='\t', index=False)
    return

def load_context_abundance():
    return pd.read_csv(os.path.join(RESOURSES_PATH, 'context_abundance.csv'),
                       sep='\t')


if __name__ == '__main__':
    context_count()
