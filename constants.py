import os

CHR_NAMES = ['chr1', 'chr2', 'chr3', 'chr4',
             'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12',
             'chr13', 'chr14', 'chr15', 'chr16',
             'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22', 'chrX', 'chrY']
SPECIES_NAMES = ['calJac3', 'canFam3', 'chlSab2', 'gorGor3', 'hg38', 'macFas5', \
                 'micMur1', 'nasLar1', 'nomLeu3', 'panPan1', 'panTro4',
                 'papAnu2', \
                 'ponAbe2', 'rheMac3', 'rhiRox1', 'saiBol1', 'tarSyr2']
REF_NAMES = ['leaf.9', 'leaf.10', 'leaf.11', 'inner.12', 'inner.13',
             'inner.14', 'inner.15', 'inner.16', 'inner.17', 'inner.18',
             'inner.19', 'inner.20']
ANCESTRAL_COLNAMES= ['chromosome', 'position', 'leaf.1', 'leaf.2', 'leaf.3',
                      'leaf.4', 'leaf.5', 'leaf.6', 'leaf.7', 'leaf.8',
                      'leaf.9', 'leaf.10', 'leaf.11', 'inner.12', 'inner.13',
                      'inner.14', 'inner.15', 'inner.16', 'inner.17',
                      'inner.18', 'inner.19', 'inner.20', 'hg18', 'rheMac2']
PRIMARY_PATH = os.path.dirname(os.path.realpath(__file__))
ALLOWED_NUCLEOTIDE_SYMBOLS = ['A', 'T', 'G', 'C', 'a', 't', 'g', 'c']

VARIANT_PATH = os.path.join(PRIMARY_PATH, 'data', 'variants')
VCF_VARIANT_PATH = os.path.join(VARIANT_PATH, 'vcf_variants')
FILT_VCF_VARIANT_PATH = os.path.join(VARIANT_PATH, 'filt_variants')
RAW_MAF_PATH = os.path.join(PRIMARY_PATH, 'data', 'maf_data')
FASTA_CHRS_PATH = os.path.join(PRIMARY_PATH, 'data', 'chr_raw_fasta_data')
LOG_PATH = os.path.join(PRIMARY_PATH, 'logs')
RESOURSES_PATH = os.path.join(PRIMARY_PATH, 'resourses')
MATRICES_PATH = os.path.join(VARIANT_PATH, 'matrix_tables')
EXOME_PATH = os.path.join(PRIMARY_PATH, 'data', 'exome_beds')
CORE_SCRIPTS_PATH = os.path.join(PRIMARY_PATH, 'core_scripts')
for path in [
    VARIANT_PATH,
    VCF_VARIANT_PATH,
    FILT_VCF_VARIANT_PATH,
    RAW_MAF_PATH,
    FASTA_CHRS_PATH,
    LOG_PATH,
    RESOURSES_PATH,
    EXOME_PATH]:
    if not os.path.exists(path):
        os.makedirs(path)

REF_STARTPOSPATH = os.path.join(RESOURSES_PATH, 'start_positions.csv')
SPECIES_NAMES_PATH = os.path.join(RESOURSES_PATH, 'species_names.csv')
OUTPUT_VCF_NAMES_PATH = os.path.join(RESOURSES_PATH, 'output_vcf_list.csv')
OUTPUT_FILT_VCF_NAMES_PATH = os.path.join(RESOURSES_PATH,
                                          'output_filtered_vcf_list.csv')
EXOME_FILE = os.path.join(EXOME_PATH, 'exome_regions.bed')
R_ORDER = ["A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
           "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
           "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
           "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
           "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
           "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
           "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
           "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
           "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
           "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
           "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
           "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
           "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
           "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
           "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
           "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"]
