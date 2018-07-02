import os
import re
dirpath = os.getcwd() 
dirpath+='/'

chr_names=['chr1','chr2','chr3','chr4',
           'chr5','chr6','chr7','chr8',
           'chr9','chr10','chr11','chr12',
           'chr13','chr14','chr15','chr16',
           'chr17','chr18','chr19','chr20',
           'chr21','chr22','chrX','chrY']
outp=open('ref_chr_startpos.txt','w')
for CHR in chr_names:
	inp = open(dirpath+CHR+'.maf', "r")
	for line in inp:
	    if "hg38" in line and line[0] == 's':
	        startpos =int(re.findall(r'\w+', line)[3])
	        firstlen =int(re.findall(r'\w+', line)[4])
	        outp.write(CHR+' '+str(startpos)+'\n')
	        break
inp.close()
outp.close()



