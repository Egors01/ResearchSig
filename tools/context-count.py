from constants import SPECIES_NAMES,CHR_NAMES
import os
import pandas as pd
import time
DIR=os.getcwd()
logfile=DIR+'/log-count.txt'



def context_count():
    chr_names = CHR_NAMES[:2]
    sp_names = SPECIES_NAMES
    contexts=[ c1+'.'+c2 for c1 in ['A','T','G','C']  for c2 in ['A','T','G','C']]
    context_abundance=pd.DataFrame()
    context_abundance['Context']=contexts
    log=open(logfile,'w')
    log.write(str(sp_names)+'\n')
    log.close()

    for name in sp_names:
        log=open(logfile,'a')
        log.write('Start time \n'+str(time.asctime(time.localtime()))+'\n')
        log.write('Start '+name+'\n')
        log.close()
        context_counts=dict.fromkeys(contexts,0)
        for CHR in chr_names: #['chr1']:
            log=open(logfile,'a')
            log.write('chr '+CHR+'\n')
            log.close()
            file_path=CHR+'/'+name+'.fasta'
            f=open(file_path,'r')
            f.readline()
            first = f.read(1).upper()
            second= f.read(1).upper()
            third = f.read(1).upper()
            while third!='' and third!=' ' and third!=None:
              if third.isalpha() and first.isalpha() and ('N' not in third+first) and ('n' not in third+first):
                context_counts[first+'.'+third]+=1
              first =second
              second=third
              third =f.read(1).upper()
            f.close()
        context_abundance[name]=context_counts.values()
    log=open(logfile,'a')
    log.write('End time \n'+str(time.asctime(time.localtime()))+'\n')
    log.close()
    context_abundance.to_csv('context_abundance.csv',sep='\t')