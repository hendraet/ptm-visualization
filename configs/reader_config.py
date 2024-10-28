# General
FASTA_FILE = 'data/uniprot_data/tau_isoforms2N4R.fasta'
INPUT_DIR = 'data/max_quant/'
ISOFORM_HELPER_DICT = { "0N3R": "P10636-2",
                        "1N3R": "P10636-4",
                        "2N3R": "P10636-5",
                        "0N4R": "P10636-6",
                        "1N4R": "P10636-7",
                        "2N4R": "P10636-8",}
ALIGNED_FASTA_FILE = "data/uniprot_data/aligned.fasta"
ISOFORM_TRANSPOSE_DICT = {
    "P10636-2": "P10636-8",
    "P10636-4": "P10636-8",
    "P10636-5": "P10636-8",
    "P10636-6": "P10636-8",
    "P10636-7": "P10636-8",
}

# Mascot

# Protein Pilot
# choose between local and global
FDR_GLOBAL = 'global'
CONFIDENCE_THRESHOLD = 0.01

# MS Fragger
MS_FRAGGER_FILE = 'data/ms_fragger/combined_modified_peptide.tsv'
MS_FRAGGER_MODS = {"42.0106": "Acetyl",
                   "79.9663": "Phospho",
                   "114.0429": "GG",
                   "14.0157": "Methyl",
                   "0.9840": "Citrullination",}

# MaxQuant
MAX_QUANT_FILE = 'data/max_quant/evidence.txt'
THRESHOLD = 0.01