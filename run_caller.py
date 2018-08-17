from core_scripts.variantcaller import variant_call_pairs
from core_scripts.gen96_from_vcf import gen96_matrix_from_filtered_vcf
from constants import RUN_ID
from tools.names_generator import create_output_file_names
from core_scripts.exome_filtering import filter_exome_from_vcf
from core_scripts.normalize_table import normilize_frequencies
PAIRS = [["nomLeu3", "ponAbe2", "hg38"]]
chr_names = ['chrY']
processed_files = variant_call_pairs(PAIRS,chr_names)
print(processed_files)
processed_files = create_output_file_names(PAIRS)
print(processed_files)
filtered_vcfs = filter_exome_from_vcf(processed_files)
gen96_matrix_from_filtered_vcf(filtered_vcfs)
normilize_frequencies()
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

