# #!/usr/bin/bash
# from constants import EXOME_FILE
# from tools.logger import Logger
# from tools.names_generator import make_vcf_filt_name
#
# echo Launched terminal part
# for file in ../data/variants/vcf_variants/*.vcf
# do
#     OUT=../data/variants/filt_variants/`basename -s .vcf $file`.filt.vcf
#     echo file=$file
# 	echo out=$OUT
# 	bedtools intersect -wao -v -a $file  -b  ../data/exome_beds/exomeD.bed >  $OUT
#
# done
import os
import subprocess

from constants import EXOME_FILE
from tools.logger import Logger
from tools.names_generator import make_vcf_filt_name


def filtering_vcf_files(processed_files):
    logger = Logger(source_name='secondary_proc',
                    msg='launched filtering part')
    filtered_filenames_list = []
    for filename in processed_files:
        logger.message('filtering pair {}'.format(os.path.basename(filename)))
        filt_vcf_filename = make_vcf_filt_name(sp_vcf_old_name=filename)
        filtered_filenames_list.append(filt_vcf_filename)
        command = "bedtools intersect -wao -v -a " + \
                  filename + " -b  " + EXOME_FILE + \
                  " > " + filt_vcf_filename
        subprocess.call(command,shell=True)
    logger.print_end()
    return filtered_filenames_list
