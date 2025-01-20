# General
FASTA_FILE = 'tests/test_data/input.fasta'

ISOFORM_HELPER_DICT = {}
GROUPS_CSV = 'tests/test_data/groups_ms_fragger.csv'

# this is the default path where the tool will safe the alignment
# just change if you want to supply your own alignment
# CAUTION: the alignment must match with the fasta file
ALIGNED_FASTA_FILE = "tests/test_data/aligned.fasta"

# MS Fragger
MS_FRAGGER_FILE = 'tests/test_data/ms_fragger/combined_modified_peptide.tsv'
MS_FRAGGER_MODS = {"42.0106": "Acetyl",
                   "79.9663": "Phospho",
                   "114.0429": "GG",
                   "14.0157": "Methyl",
                   "0.9840": "Citrullination",}
