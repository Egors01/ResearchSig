import os
import subprocess

import constants
from tools.logger import Logger

from tools.names_generator import make_vcf_filt_name

def filtering_vcf_files(processed_files,run_id=None):
    logger = Logger(source_name='filtering ',
                    msg='launched filtering part',run_id=run_id)
    filtered_filenames_list = []
    for filename in processed_files:
        logger.message('filtering pair {}'.format(os.path.basename(filename)))
        filt_vcf_filename = make_vcf_filt_name(sp_vcf_old_name=filename)
        filtered_filenames_list.append(filt_vcf_filename)
        try:

            command = "bedtools intersect -v -wao -a " + \
                      filename + " -b  " + constants.EXOME_FILE + \
                      " > " + filt_vcf_filename
            subprocess.call(command,shell=True)
        except:
            logger.message('Error in filtering ')
            raise Exception('Error in filetering ')

    logger.print_end()
    return filtered_filenames_list
