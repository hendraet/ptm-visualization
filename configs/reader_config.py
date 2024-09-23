# General
FASTA_FILE = 'data/uniprot_data/tau_isoforms2N4R.fasta'
ALIGNED_FASTA_FILE = 'data/uniprot_data/tau_aligned.fasta'
INPUT_DIR = 'data/protein_pilot/'

# Mascot
isoform_helper_dict = { "0N3R": "P10636-2",
                        "1N3R": "P10636-4",
                        "2N3R": "P10636-5",
                        "0N4R": "P10636-6",
                        "1N4R": "P10636-7",
                        "2N4R": "P10636-8",}

# Protein Pilot
# choose between local and global
FDR_GLOBAL = 'global'
CONFIDENCE_THRESHOLD = 0.01