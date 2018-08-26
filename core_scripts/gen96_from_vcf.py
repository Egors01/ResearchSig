import sys
import os
import pandas as pd
import argparse

from constants import MATRICES_PATH
from tools.logger import Logger


def gen96_matrix_from_filtered_vcf(filt_vcf_names_list,run_id=None):
    result = pd.DataFrame()
    MATRIX_NAME = os.path.join(MATRICES_PATH,
                               '96_matrix.run_' + run_id + '.csv')
    logger=Logger(source_name='table_generator',msg='Launched table creation',run_id=run_id)
    for filename in filt_vcf_names_list:
        logger.message(msg='started for {}'.format(os.path.basename(filename)))
        vcf_data=pd.read_csv(filename, sep='\t', skiprows=0,header=None,names=["CHROM","POS", "ID", "REF", "ALT","QUAL", "FILTER","INFO"])
        vcf_data['Sample']=vcf_data['ID'].str.split('_',expand=True)[1]
        vcf_data['Context']=vcf_data['ID'].str.split('_',expand=True)[3]
        vcf_data['SUBS']=vcf_data['ID'].str.split('_',expand=True)[2]
        #name = filename.split(".filt.vcf")[0].split("/")[1]
        name = filename.split("/")[-1].split('.filt.vcf')[0]
        print(name)
        df = vcf_data.groupby(['SUBS','Context']).size().reset_index(name=name)
        df.sort_values(by=['SUBS'],inplace=True)
        result[name]=df[name]
        logger.message(msg='Done for '+name+ ' \n')
    #append subs coulum by taking from the last df
    df.sort_values(by=['SUBS'],inplace=True)
    df['SUBS'].replace(r'\B','\1>', inplace=True,regex=True)
    df['SUBS'] = df.apply(lambda row: row.Context[0]+'['+row.SUBS +']'+row.Context[2] ,axis=1)
    result['SUBS']=df['SUBS']
    result.to_csv(MATRIX_NAME, sep='\t',index=False)
    return MATRIX_NAME

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




