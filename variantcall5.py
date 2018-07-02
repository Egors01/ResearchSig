import re
import os
import pandas as pd
import time
import gc
import shutil
#DEFINE OS

WIN=0
CSV=1
CHRNAME="example"
#PROVIDES GLOBAL VARS
dirpath,outputpath= configure_dir()
log=open(dirpath+'log2.txt','w')
log.write('Global outputpath'+outputpath+'\n')
log.write('Global dirpath'+dirpath+'\n')
log.close()
def get_names(chrname=CHRNAME,dirpath=dirpath):
    dirpath+=chrname
    inp = open(dirpath, "r")
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
def get_ref_startpos(chrname=CHRNAME,dirpath=dirpath):
    dirpath+=chrname
    inp = open(dirpath+'.maf', "r")
    for line in inp:
        if "hg38" in line and line[0] == 's':
            startpos =int(re.findall(r'\w+', line)[3])
            firstlen =int(re.findall(r'\w+', line)[4])
            break
    inp.close()
    return startpos

def get_context(f,fileindex):
    f[fileindex].seek(-2,1)
    nuc_prev=f[fileindex].read(1)
    f[fileindex].seek(1,1)
    nuc_next=f[fileindex].read(1)
    f[fileindex].seek(-1,1)
    nuc_next= str(nuc_next,'utf-8')
    nuc_prev= str(nuc_prev,'utf-8')
    return nuc_prev+'.'+nuc_next

def chr_variant_call(CHRNAME,sp_name_1="nomLeu3",sp_name_2="ponAbe2",ref_name='hg38',dirpath=dirpath):
    outputpath=dirpath+CHRNAME+'/'
    print('Variants counting for '+CHRNAME+ ' ' + sp_name_1+' '+sp_name_2)
    sp_names=[ref_name,sp_name_1,sp_name_2]
    f=[open(outputpath+name+".fasta",'rb') for name in sp_names ]
    for i in range(0,3):
        f[i].readline()
    docpos=0
    alt_name=''
    [ALT,REF,SP1,SP2,REFPOS,INFO,CHR,CONTEXT]=[[],[],[],[],[],[],[],[]]
    refpos=get_ref_startpos(CHRNAME)
    #print('IN CALL2')
    while True:
        sp_1 = f[1].read(1)
        sp_2 = f[2].read(1)
        ref  = f[0].read(1)
        sp_1 = str(sp_1,'utf-8').upper()
        sp_2 = str(sp_2,'utf-8').upper()
        ref  = str(ref,'utf-8').upper()
        #print(sp_1,sp_2)
        refpos+=1
        alt='Q'
        docpos+=1
        if sp_1 != sp_2:
            if not (sp_1!=ref and sp_2!=ref):
                if ('-' not in sp_1+sp_2) and ('n' not in sp_1+sp_2) and ('N' not in sp_1+sp_2) and ('*' not in sp_1+sp_2):
                    if sp_1== ref:
                        #print("sp-1")
                        alt = sp_2
                        alt_name = sp_names[2]
                        context=get_context(f,2)
                    elif sp_2 == ref:
                        alt = sp_1
                        alt_name = sp_names[1]
                        context=get_context(f,1)
                    if ('-' not in context and 'N' not in context and  'n' not in context and '*' not in context):
                        #print("app sector")
                        info="dp "+ str(docpos)
                        CONTEXT.append(context.upper())
                        CHR.append(CHRNAME+'_'+alt_name)
                        ALT.append(alt)
                        REF.append(ref)
                        REFPOS.append(refpos)
                        INFO.append(info)
        #print('IN CALL3')
        if not sp_1 or not sp_2:
            print('break at ',refpos,docpos)
            break
    SUBS =[nref.upper()+'>'+nmut.upper() for nref,nmut in zip(REF,ALT)]
    for index,sbs in enumerate(SUBS):
        if sbs == 'A>T':
            SUBS[index]='T>A'
            INFO[index]+=('CPL A>T')
        if sbs == 'G>C':
            SUBS[index]='C>G'
            INFO[index]+=('CPL G>C')
        if sbs == 'G>A':
            SUBS[index]='C>T'
            INFO[index]+=('CPL G>A')
        if sbs == 'A>G':
            SUBS[index]='T>C'
            INFO[index]+=('CPL A>G')
        if sbs == 'A>C':
            SUBS[index]='T>G'
            INFO[index]+=('CPL A>C')
        if sbs == 'G>T':
            SUBS[index]='C>A'
            INFO[index]+=('CPL G>T')
      
    dictdata={'Chr':CHR,'SUBS':SUBS,'Alt' : ALT,'Ref' : REF,'Refpos':REFPOS,'Context':CONTEXT,'Info':INFO}
    chrdata=pd.DataFrame(dictdata,columns=dictdata.keys())
    data_sp_1=pd.DataFrame(chrdata.loc[chrdata['Chr'] == (CHRNAME+'_'+sp_name_1)])
    data_sp_2=pd.DataFrame(chrdata.loc[chrdata['Chr'] == (CHRNAME+'_'+sp_name_2)])
        
    return data_sp_1,data_sp_2

def prep_print(result,sp_name,sp_pair,refname,dirpath=dirpath):

    final_pair_path=dirpath+'/Variants-new/'
    if not os.path.exists(final_pair_path):
        os.makedirs(final_pair_path)
    filename=final_pair_path+sp_name+'_VS_'+sp_pair+'_r.'+refname[:4]+'.vcf'
    print(filename)
    header='##fileformat=VCFv4.1 \n'
    header+=('##contig=<ID='+filename+', length=111111,assembly=hg38> \n##reference='+refname+ '\n')
    header+='#CHROM POS ID REF ALT QUAL FILTER INFO\n'
    f=open(filename,'w')
    f.write(header)
    f.close()

        
    del result['Info']
    result["CHROM"]=result['Chr'].replace(r'_\w+',"", inplace=False,regex=True)
    result['Chr']=result["Chr"]+"_"+result["SUBS"]+"_"+result["Context"]
    result['QUAL']='.'
    result['FILTER']='.'
    result['INFO']='.'
    del result['Context']
    result.rename(index=str, columns={"Chr": "ID", "Refpos": "POS","Ref":"REF","Alt":"ALT"},inplace=True)
    result=result[["CHROM", "POS", "ID", "REF", "ALT","QUAL", "FILTER","INFO"]]
    with open(filename, 'a') as f:
        result.to_csv(f, header=False,sep='\t',index=False)
    return #result
            
PAIRS=[["nomLeu3","ponAbe2","hg38"]]
       
chr_names=['chr1','chr2','chr3','chr4',
           'chr5','chr6','chr7','chr8',
           'chr9','chr10','chr11','chr12',
           'chr13','chr14','chr15','chr16',
           'chr17','chr18','chr19','chr20',
           'chr21','chr22']

start_a=time.time()

log=open(dirpath+'log2.txt','a')
log.write('Start time \n'+str(time.asctime(time.localtime()))+'\n')
log.close()
for pair in PAIRS:
    log=open(dirpath+'log2.txt','a')
    log.write('\n**********\nStart pair \n'+str(pair)+'\n')
    log.close()
    start_c=time.time()
    v_df1=pd.DataFrame()
    v_df2=pd.DataFrame()
    for CHRNAME in chr_names:
        
        log=open(dirpath+'log2.txt','a')
        log.write('Start chr '+CHRNAME+'\n')
        start_b=time.time()
        log.close()
        
        chrv1,chrv2=chr_variant_call(CHRNAME,pair[1],pair[2],pair[0])
        
        log=open(dirpath+'log2.txt','a')
        log.write("Finished "+ str(round((time.time() - start_b),3) )+' seconds \n')
        log.close()
        
        v_df1 = pd.concat([ chrv1, v_df1], ignore_index=True)
        v_df2 = pd.concat([ chrv2, v_df2], ignore_index=True)
    
    log=open(dirpath+'log2.txt','a')
    log.write('Start print 1 \n')
    start_b=time.time()
    prep_print(v_df1,pair[1],pair[2],pair[0],dirpath)
    log.write("print finished "+ str(round((time.time() - start_b),3) )+' seconds \n')
    log.write("Time per pair "+ str(round((time.time() - start_c),3) )+' seconds \n')
    log.close()
    
    log=open(dirpath+'log2.txt','a')
    log.write('Start print 2 \n')
    start_b=time.time()
    prep_print(v_df2,pair[2],pair[1],pair[0],dirpath)
    log.write("print finished "+ str(round((time.time() - start_b),3) )+' seconds \n')
    log.write('End time \n'+str(time.asctime(time.localtime()))+'\n')
    log.close()

print('end of run')
