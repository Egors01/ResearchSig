from core_scripts.variantcaller import variant_call_pairs
from core_scripts.gen96_from_vcf import gen96_matrix_from_filtered_vcf
from constants import RUN_ID
from tools.names_generator import create_output_file_names

PAIRS = [["nomLeu3", "ponAbe2", "hg38"]]
chr_names = ['chrY']
create_output_file_names(PAIRS)
#variant_call_pairs(PAIRS,chr_names)
#gen96_matrix_from_filtered_vcf(PAIRS)

    # env = Environment()
    # var_log = Logger('variant_call','vclog')
    # for i in range(1, 20000000):
    #     a = 8 * i
    # var_log.tick()
    # var_log.message('in the middle ', print_time=True)
    # for i in range(1, 20000000):
    #     a = 8 * i
    # var_log.message('at the end ', print_time=True)
    # var_log.print_end()

