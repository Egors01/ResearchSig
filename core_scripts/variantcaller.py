import re
import os
import pandas as pd
import time
import gc
import shutil
from tools.set_env import Environment
from tools.logger import Logger
from tools.find_refstartpos import get_ref_startpos
from tools.list_names import get_names
import constants

#wrapper to extract SNP's     
def variant_call_pairs(pairs, chr_to_select):
    mainlog = Logger(reset=True)
    chr_names = chr_to_select
    for pair in pairs:
        mainlog.tick(timer_name='pair_timer')
        mainlog.message(msg='Started pair ' + str(pair))
        v_df1, v_df2 = [pd.DataFrame(), pd.DataFrame()]
        for chrname in chr_names:
            mainlog.tick(timer_name='chr_timer')
            mainlog.message(msg='Start chr ' + chrname )
            chrv1, chrv2 = chr_variant_call(chrname, pair[1], pair[2], pair[0])
            mainlog.message(msg="finished chr "+chrname, print_time=True, timer_name='chr_timer')
            v_df1 = pd.concat([chrv1, v_df1], ignore_index=True)
            v_df2 = pd.concat([chrv2, v_df2], ignore_index=True)
        mainlog.message(source_name='variant_call_pairs', msg='Ended pair ' + str(pair), print_time=True,
                        timer_name='pair_timer')
        df_to_vcf(v_df1, pair[1], pair[2], pair[0])
        df_to_vcf(v_df2, pair[2], pair[1], pair[0])
    mainlog.print_end()

def get_context(f, fileindex):
    f[fileindex].seek(-2, 1)
    nuc_prev = f[fileindex].read(1)
    f[fileindex].seek(1, 1)
    nuc_next = f[fileindex].read(1)
    f[fileindex].seek(-1, 1)
    nuc_next = str(nuc_next, 'utf-8')
    nuc_prev = str(nuc_prev, 'utf-8')
    return nuc_prev.upper() + '.' + nuc_next.upper()


def chr_variant_call(chr_name, sp_name_1="", sp_name_2="", ref_name=''):
    outputpath = os.path.join(Environment().fasta_chrs_path, chr_name)
    variant_call_logger = Logger(source_name='variant_call')
    #                             msg='variants counting for ' + chr_name + ' ' + sp_name_1 + ' ' + sp_name_2)
    sp_names = [ref_name, sp_name_1, sp_name_2]
    f = [open(os.path.join(outputpath, name + ".fasta"), 'rb') for name in sp_names]
    for i in range(0, 3):
        f[i].readline()
    docpos = 0
    alt_name = ''
    [ALT, REF, SP1, SP2, REFPOS, INFO, CHR, CONTEXT] = [[], [], [], [], [], [], [], []]
    refpos = get_ref_startpos(chr_name)
    # print('IN CALL2')
    while True:
        sp_1 = f[1].read(1)
        sp_2 = f[2].read(1)
        ref = f[0].read(1)
        sp_1 = str(sp_1, 'utf-8').upper()
        sp_2 = str(sp_2, 'utf-8').upper()
        ref = str(ref, 'utf-8').upper()
        # print(sp_1,sp_2)
        refpos += 1
        alt = 'Q'
        docpos += 1
        if sp_1 != sp_2:
            if not (sp_1 != ref and sp_2 != ref):
                if bool(re.match(r'[ATGCatgc]', sp_1)) and bool(re.match(r'[ATGCatgc]', sp_2)):
                    # if (sp_1 in constants.ALLOWED_NUCLEOTIDE_SYMBOLS) and (sp_2 in constants.ALLOWED_NUCLEOTIDE_SYMBOLS):
                    # if ('-' not in sp_1+sp_2) and ('n' not in sp_1+sp_2) and ('N' not in sp_1+sp_2) and ('*' not in sp_1+sp_2):
                    if sp_1 == ref:
                        # print("sp-1")
                        alt = sp_2
                        alt_name = sp_names[2]
                        context = get_context(f, 2)
                    elif sp_2 == ref:
                        alt = sp_1
                        alt_name = sp_names[1]
                        context = get_context(f, 1)
                    else:
                        variant_call_logger.message(msg='Error when comparing ')
                    if bool(re.match(r'[ATGCatgc].[ATGCatgc]', context)):
                        # if ('-' not in context and 'N' not in context and  'n' not in context and '*' not in context):
                        # print("app sector")
                        info = "dp " + str(docpos)
                        CONTEXT.append(context.upper())
                        CHR.append(chr_name + '_' + alt_name)
                        ALT.append(alt)
                        REF.append(ref)
                        REFPOS.append(refpos)
                        INFO.append(info)
        # print('IN CALL3')
        if not sp_1 or not sp_2:
            break
    SUBS = [nref.upper() + '>' + nmut.upper() for nref, nmut in zip(REF, ALT)]
    for index, sbs in enumerate(SUBS):
        if sbs == 'A>T':
            SUBS[index] = 'T>A'
            INFO[index] += ('CPL A>T')
        if sbs == 'G>C':
            SUBS[index] = 'C>G'
            INFO[index] += ('CPL G>C')
        if sbs == 'G>A':
            SUBS[index] = 'C>T'
            INFO[index] += ('CPL G>A')
        if sbs == 'A>G':
            SUBS[index] = 'T>C'
            INFO[index] += ('CPL A>G')
        if sbs == 'A>C':
            SUBS[index] = 'T>G'
            INFO[index] += ('CPL A>C')
        if sbs == 'G>T':
            SUBS[index] = 'C>A'
            INFO[index] += ('CPL G>T')

    dictdata = {'Chr': CHR, 'SUBS': SUBS, 'Alt': ALT, 'Ref': REF, 'Refpos': REFPOS, 'Context': CONTEXT, 'Info': INFO}
    chrdata = pd.DataFrame(dictdata, columns=dictdata.keys())
    data_sp_1 = pd.DataFrame(chrdata.loc[chrdata['Chr'] == (chr_name + '_' + sp_name_1)])
    data_sp_2 = pd.DataFrame(chrdata.loc[chrdata['Chr'] == (chr_name + '_' + sp_name_2)])

    return data_sp_1, data_sp_2


def df_to_vcf(result, sp_name, sp_pair, refname):
    final_pair_path = os.path.join(Environment().variant_path,'raw_variants')
    if not os.path.exists(final_pair_path):
        os.makedirs(final_pair_path)
    filename = os.path.join(final_pair_path, sp_name + '_VS_' + sp_pair + '_r.' + refname[:4] + '.vcf')
    header = '##fileformat=VCFv4.1 \n'
    header += ('##contig=<ID=' + filename + ', length=111111,assembly=hg38> \n##reference=' + refname + '\n')
    header += '#CHROM POS ID REF ALT QUAL FILTER INFO\n'
    f = open(filename, 'w')
    f.write(header)
    f.close()

    del result['Info']
    result["CHROM"] = result['Chr'].replace(r'_\w+', "", inplace=False, regex=True)
    result['Chr'] = result["Chr"] + "_" + result["SUBS"] + "_" + result["Context"]
    result['QUAL'] = '.'
    result['FILTER'] = '.'
    result['INFO'] = '.'
    del result['Context']
    result.rename(index=str, columns={"Chr": "ID", "Refpos": "POS", "Ref": "REF", "Alt": "ALT"}, inplace=True)
    result = result[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]]
    f=open(filename,'r')
    HH = []
    for line in f:
        HH.append(line.strip())
    print(HH)
    f.close()
    with open(filename, 'a') as f:
        result.to_csv(f, header=False, sep='\t', index=False,mode='a')
    Logger(source_name='to_vcf', msg='finished vcf')
    return  # result
