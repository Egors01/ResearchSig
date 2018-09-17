import argparse

from core_scripts.exome_filtering import filtering_vcf_files
from core_scripts.gen96_from_vcf import  \
    gen192_matrix_from_filtered_vcf
from core_scripts.matrix_operations import normalize_frequencies_192_context, \
    reduce_192_to_96_common_notation, matrix_to_r_output
from parameters import input_data
from core_scripts.variantcaller import variant_call_pairs
from tools.names_generator import create_output_file_names, make_vcf_filt_names


def pipeline(run_id):
    chr_names = ['chr1', 'chr2', 'chr3', 'chr4',
             'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12',
             'chr13', 'chr14', 'chr15', 'chr16',
             'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22']
    #chr_names=['chrX']#
    for runid in [0]:
        run_id=str(runid)
        pairs_list = input_data.get_pairs_by_run_id(run_id=run_id)

        processed_files = variant_call_pairs(pairs_list, chr_names, run_id=run_id)
        processed_files = create_output_file_names(pairs_list)

        filtered_vcfs = filtering_vcf_files(processed_files, run_id=run_id)
        filtered_vcfs = make_vcf_filt_names(processed_files)
        print(filtered_vcfs)

        matrix_name = gen192_matrix_from_filtered_vcf(
            filt_vcf_names_list=filtered_vcfs, run_id=run_id)

        norm_matrix_name = normalize_frequencies_192_context(path_to_192_matrix=matrix_name, run_id=run_id)

        m96_matrix_name = reduce_192_to_96_common_notation(path_to_192_matrix=norm_matrix_name, run_id=run_id)

        matrix_to_r_output(path_to_96_matrix= m96_matrix_name, run_id='999')

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


if __name__ == '__main__':
    RUN_ID = run_wrapper()
    pipeline(RUN_ID)
    pass
