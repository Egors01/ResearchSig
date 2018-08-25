import pybedtools

from constants import EXOME_FILE
from tools.logger import Logger
from tools.names_generator import make_vcf_filt_name


def filter_exome_from_vcf(vcf_filenames_full):
    logger = Logger(source_name='exome_filter', msg='launched exome filter')
    filtered_vcf_filenames = []
    for sp_vcf_name in vcf_filenames_full:
        logger.tick(timer_name='filterer')

        bed = pybedtools.BedTool(EXOME_FILE)
        vcf = pybedtools.BedTool(sp_vcf_name)
        inters = vcf.intersect(bed, v=True, wao=True)
        new_vcf_name = make_vcf_filt_name(sp_vcf_old_name = sp_vcf_name)
        inters.moveto(new_vcf_name)

        logger.message(msg='filtered pair {}'.format(new_vcf_name))
        logger.message(msg='filtering time for chr ',print_time=True,timer_name='filterer')

        filtered_vcf_filenames.append(new_vcf_name)
    return filtered_vcf_filenames