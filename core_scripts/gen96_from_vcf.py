import sys
import os
import pandas as pd
import argparse

from constants import MATRIX_NAME
from tools.logger import Logger
from tools.names_generator import create_vcf_name_from_pair_and_ref

n_files=len(sys.argv[1:])
result=pd.DataFrame()

def old():
    log=open('genlog.txt','w')
    log.write('Accepted args '+str(sys.argv)+ '\nStart matrix gen \n')
    log.write('name example '+sys.argv[2].split(".filt.vcf")[0].split("/")[1]+ ' \n')
    log.close()

    for filename in sys.argv[1:]:
        print('here')
        vcf_data=pd.read_csv(filename, sep='\t', skiprows=4,header=None,names=["CHROM","POS", "ID", "REF", "ALT","QUAL", "FILTER","INFO"])
        vcf_data['Sample']=vcf_data['ID'].str.split('_',expand=True)[1]
        vcf_data['Context']=vcf_data['ID'].str.split('_',expand=True)[3]
        vcf_data['SUBS']=vcf_data['ID'].str.split('_',expand=True)[2]
        print('For filename ',filename)
        name=filename.split(".filt.vcf")[0].split("/")[1]
        df=vcf_data.groupby(['SUBS','Context']).size().reset_index(name=name)
        df.sort_values(by=['SUBS'],inplace=True)
        result[name]=df[name]
        log=open('genlog.txt','a')
        log.write('Done for '+name+ ' \n')
        log.close()

    #append subs coulum by taking from the last df
    df.sort_values(by=['SUBS'],inplace=True)
    df['SUBS'].replace(r'\B','\1>', inplace=True,regex=True)
    df['SUBS'] = df.apply(lambda row: row.Context[0]+'['+row.SUBS +']'+row.Context[2] ,axis=1)
    result['SUBS']=df['SUBS']
    result.to_csv('96-Matrix-recent.csv', sep='\t',index=False)

def gen96_matrix_from_filtered_vcf(PAIRS_AND_REF):
    logger=Logger(source_name='table_generator',msg='Started table creation')
    filt_vcf_names_list=[]
    for pair_and_ref in PAIRS_AND_REF:
        filt_vcf_name=create_vcf_name_from_pair_and_ref(sp_name=pair_and_ref[0],
                                                         sp_pair=pair_and_ref[1],
                                                         refname=pair_and_ref[2])
        filt_vcf_names_list.append(filt_vcf_name)
        filt_vcf_name=create_vcf_name_from_pair_and_ref(sp_name=pair_and_ref[1],
                                                         sp_pair=pair_and_ref[0],
                                                         refname=pair_and_ref[2])
        filt_vcf_names_list.append(filt_vcf_name)
    logger.message(msg='Accepted args '+str(filt_vcf_names_list))
    print('sp names in 96m',filt_vcf_names_list)
#    log.write('name example '+sys.argv[2].split(".filt.vcf")[0].split("/")[1]+ ' \n')
#    log.close()

    for filename in filt_vcf_names_list:

        vcf_data=pd.read_csv(filename, sep='\t', skiprows=4,header=None,names=["CHROM","POS", "ID", "REF", "ALT","QUAL", "FILTER","INFO"])
        vcf_data['Sample']=vcf_data['ID'].str.split('_',expand=True)[1]
        vcf_data['Context']=vcf_data['ID'].str.split('_',expand=True)[3]
        vcf_data['SUBS']=vcf_data['ID'].str.split('_',expand=True)[2]
        print('For filename ',filename)
        #name = filename.split(".filt.vcf")[0].split("/")[1]
        name = filename.split("/")[-1].split('.vcf')[0]
        print(name)
        df=vcf_data.groupby(['SUBS','Context']).size().reset_index(name=name)
        df.sort_values(by=['SUBS'],inplace=True)
        result[name]=df[name]
        logger.message(msg='Done for '+name+ ' \n')

    #append subs coulum by taking from the last df
    df.sort_values(by=['SUBS'],inplace=True)
    df['SUBS'].replace(r'\B','\1>', inplace=True,regex=True)
    df['SUBS'] = df.apply(lambda row: row.Context[0]+'['+row.SUBS +']'+row.Context[2] ,axis=1)
    result['SUBS']=df['SUBS']
    result.to_csv(MATRIX_NAME, sep='\t',index=False)


def gen96_matrix_from_vcf_script():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fileslist",
        help="vcf files list",
        nargs='+',
        type=str,
        required=True,
        action="store"
    )
    args = parser.parse_args()

    input_file_list = args.fileslist
    print(input_file_list)

if __name__ == "__main__":
    gen96_matrix_from_vcf_script()



