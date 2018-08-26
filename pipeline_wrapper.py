import argparse

from core_scripts.exome_filtering import filtering_vcf_files
from core_scripts.gen96_from_vcf import gen96_matrix_from_filtered_vcf
from core_scripts.normalize_table import normilize_frequencies
from parameters import input_data
from core_scripts.variantcaller import variant_call_pairs

def pipeline(run_id):
    chr_names = ['chrY']

    pairs_list = input_data.get_pairs_by_run_id(run_id=run_id)

    processed_files = variant_call_pairs(pairs_list, chr_names,run_id=run_id)
    filtered_vcfs = filtering_vcf_files(processed_files,run_id=run_id)

    matrix_name = gen96_matrix_from_filtered_vcf(filt_vcf_names_list=filtered_vcfs,run_id=run_id)
    normilize_frequencies(path_to_96_matrix=matrix_name,run_id=run_id)
    print('completed')
    return

def run_wrapper():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--run_id",
        help="id of run",
        type=int,
        required=True,
        action="store"
    )
    args = parser.parse_args()
    RUN_ID = str(args.run_id)
    return RUN_ID

if __name__=='__main__':
    RUN_ID = run_wrapper()
    pipeline(RUN_ID)
    pass
