import pybedtools
from tools.logger import Logger

def filter_exome_from_vcf(filenames_full):
    logger = Logger(source_name='exome_filter', msg='Launched exome filter')
    for sp_name in filenames_full:
        vcf = pybedtools.BedTool(EXOME_FILE)
        bed = pybedtools.BedTool(sp_name)
        inters = vcf.intersect(bed, v=True, wao=True)