# General
FASTA_FILE = 'tests/test_data/input.fasta'

ISOFORM_HELPER_DICT = { "0N3R": "P10636-2",
                        "1N3R": "P10636-4",
                        "2N3R": "P10636-5",
                        "0N4R": "P10636-6",
                        "1N4R": "P10636-7",
                        "2N4R": "P10636-8",}
GROUPS_CSV = 'tests/test_data/groups.csv'

# this is the default path where the tool will safe the alignment
# just change if you want to supply your own alignment
# CAUTION: the alignment must match with the fasta file
ALIGNED_FASTA_FILE = "tests/test_data/aligned.fasta"

# Mascot
MASCOT_INPUT_DIR = 'tests/test_data/mascot/'

# Protein Pilot
# choose between local and global
PROTEIN_PILOT_INPUT_DIR = 'tests/test_data/protein_pilot/'
FDR_GLOBAL = 'global'
CONFIDENCE_THRESHOLD = 0.1
# choose between 'all' and 'protein'
RELEVANT_MODS = 'all'

# MS Fragger
MS_FRAGGER_FILE = 'tests/test_data/ms_fragger/combined_modified_peptide.tsv'
MS_FRAGGER_MODS = {"42.0106": "Acetyl",
                   "79.9663": "Phospho",
                   "114.0429": "GG",
                   "14.0157": "Methyl",
                   "0.9840": "Citrullination",}

# MaxQuant
MAX_QUANT_FILE = 'tests/test_data/evidence.txt'
THRESHOLD = 0.01
