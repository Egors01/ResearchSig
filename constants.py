import os

OS_TYPE = 'LINUX'
RUN_ID = '5'
CHR_NAMES = ['chr1', 'chr2', 'chr3', 'chr4',
             'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12',
             'chr13', 'chr14', 'chr15', 'chr16',
             'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22', 'chrX', 'chrY']
SPECIES_NAMES = ['calJac3', 'canFam3', 'chlSab2', 'gorGor3', 'hg38', 'macFas5', \
                 'micMur1', 'nasLar1', 'nomLeu3', 'panPan1', 'panTro4', 'papAnu2', \
                 'ponAbe2', 'rheMac3', 'rhiRox1', 'saiBol1', 'tarSyr2']

PRIMARY_PATH = os.path.dirname(os.path.realpath(__file__))
ALLOWED_NUCLEOTIDE_SYMBOLS = ['A', 'T', 'G', 'C', 'a', 't', 'g', 'c']

VARIANT_PATH = os.path.join(PRIMARY_PATH, 'data', 'variants')
VARIANT_PATH_VCF = os.path.join(VARIANT_PATH, 'vcf_variants')
RAW_MAF_PATH = os.path.join(PRIMARY_PATH, 'data', 'maf_data')
FASTA_CHRS_PATH = os.path.join(PRIMARY_PATH, 'data', 'chr_raw_fasta_data')
LOG_PATH = os.path.join(PRIMARY_PATH, 'logs')
RESOURSES_PATH = os.path.join(PRIMARY_PATH, 'resourses')
MATRICES_PATH = os.path.join(VARIANT_PATH, 'matrix_tables')
EXOME_PATH=os.path.join(PRIMARY_PATH, 'data', 'exome_beds')
for path in [
    VARIANT_PATH,
    VARIANT_PATH_VCF,
    RAW_MAF_PATH,
    FASTA_CHRS_PATH,
    LOG_PATH,
    RESOURSES_PATH,
    EXOME_PATH]:
    if not os.path.exists(path):
        os.makedirs(path)

MATRIX_NAME = os.path.join(MATRICES_PATH, '96_matrix.run_' + RUN_ID + '.csv')
REF_STARTPOSPATH = os.path.join(RESOURSES_PATH, 'start_positions.csv')
SPECIES_NAMES_PATH = os.path.join(RESOURSES_PATH, 'species_names.csv')
OUTPUT_FILT_VCF_NAMES_PATH = os.path.join(RESOURSES_PATH, 'output_vcf_list.csv')
OUTPUT_VCF_NAMES_PATH = os.path.join(RESOURSES_PATH, 'output_filtered_vcf_list.csv')
EXOME_FILE = os.path.join(EXOME_PATH, 'ExomeD.bed')