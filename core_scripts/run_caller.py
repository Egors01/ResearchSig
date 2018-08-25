from constants import RUN_ID
from core_scripts.exome_filtering import filtering_vcf_files
from core_scripts.gen96_from_vcf import gen96_matrix_from_filtered_vcf
from core_scripts.normalize_table import normilize_frequencies
from parameters import input_data
from core_scripts.variantcaller import variant_call_pairs

chr_names = ['chrY']

pairs_list = input_data.get_pairs_by_run_id(RUN_ID)

processed_files = variant_call_pairs(pairs_list, chr_names)
filtered_vcfs = filtering_vcf_files(processed_files)

matrix_name = gen96_matrix_from_filtered_vcf(filtered_vcfs)
normilize_frequencies(matrix_name)
