# General
FASTA_FILE = 'tests/test_data/input.fasta'

ISOFORM_HELPER_DICT = {}
GROUPS_CSV = 'tests/test_data/groups_protein_pilot.csv'

# this is the default path where the tool will safe the alignment
# just change if you want to supply your own alignment
# CAUTION: the alignment must match with the fasta file
ALIGNED_FASTA_FILE = "tests/test_data/aligned.fasta"

# Protein Pilot
# choose between local and global
PROTEIN_PILOT_INPUT_DIR = 'tests/test_data/protein_pilot/'
FDR_GLOBAL = 'global'
CONFIDENCE_THRESHOLD = 0.1
# choose between 'all' and 'protein'
RELEVANT_MODS = 'all'
