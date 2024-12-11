# General
#FASTA_FILE = 'data/uniprot_data/P14136_a_e.fasta'
FASTA_FILE = '/home/talnawa/Desktop/protein_sequencing/data/experiment/test.fasta'
ISOFORM_HELPER_DICT = { "0N3R": "P10636-2",
                        "1N3R": "P10636-4",
                        "2N3R": "P10636-5",
                        "0N4R": "P10636-6",
                        "1N4R": "P10636-7",
                        "2N4R": "P10636-8",}
#GROUPS_CSV = 'data/groups.csv'
GROUPS_CSV = '/home/talnawa/Desktop/protein_sequencing/data/experiment/groups.csv'
# this is the default path where the tool will safe the alignment
# just change if you want to supply your own alignment
# CAUTION: the alignment must match with the fasta file
ALIGNED_FASTA_FILE = "data/uniprot_data/aligned.fasta"

# Mascot
MASCOT_INPUT_DIR = 'data/mascot/'

# Protein Pilot
# choose between local and global
PROTEIN_PILOT_INPUT_DIR = 'data/protein_pilot_gfap/P4/'
FDR_GLOBAL = 'global'
CONFIDENCE_THRESHOLD = 0.1
# choose between 'all' and 'protein'
RELEVANT_MODS = 'all'

# MS Fragger
MS_FRAGGER_FILE = 'data/ms_fragger/combined_modified_peptide.tsv'
MS_FRAGGER_MODS = {"42.0106": "Acetyl",
                   "79.9663": "Phospho",
                   "114.0429": "GG",
                   "14.0157": "Methyl",
                   "0.9840": "Citrullination",}

# MaxQuant
#MAX_QUANT_FILE = 'data/experiment/evidence.txt'
MAX_QUANT_FILE = '/home/talnawa/Desktop/protein_sequencing/data/experiment/evidence_output.txt'
THRESHOLD = 0.01