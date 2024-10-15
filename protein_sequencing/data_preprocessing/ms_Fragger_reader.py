import importlib
import os
import pandas as pd
from protein_sequencing.data_preprocessing import reader_helper
from typing import Tuple

# TODO make changable from main input script

CONFIG = importlib.import_module('configs.default_config', 'configs')
READER_CONFIG = importlib.import_module('configs.reader_config', 'configs')

fasta_file = READER_CONFIG.FASTA_FILE
aligned_fasta_file = READER_CONFIG.ALIGNED_FASTA_FILE
input_dir = READER_CONFIG.INPUT_DIR

groups_df = pd.read_csv(f"{os.path.dirname(__file__)}/groups.csv")

sorted_isoform_headers = reader_helper.process_tau_file(fasta_file, aligned_fasta_file)

def get_accession(accession: str, peptide: str) -> Tuple[str, str, int]:
    offset = 0
    for header in sorted_isoform_headers:
        if peptide in header[1]:
            isoform = header[0].strip()
            sequence = header[1]
            offset = sequence.index(peptide)
            break
    if sequence is not None:
        return isoform, sequence, offset
    else:
        raise ValueError(f"Peptide {peptide} with accession {accession} not found in fasta file")

def check_N_term_cleavage(sequence: str, peptide: str, accession: str) -> str:
    index = sequence.index(peptide)
    amino_acid_first = peptide[0]
    amino_acid_before = ""
    if index > 0:
        amino_acid_before = sequence[index - 1]
    if amino_acid_before != "K" and amino_acid_before != "R":
        return f"N-term({amino_acid_first}){index+1}"

    return ""

def check_C_term_cleavage(sequence: str, peptide: str, accession: str) -> str:
    index = sequence.index(peptide) + len(peptide)
    amino_acid_last = peptide[-1]
    if amino_acid_last != "K" and amino_acid_last != "R":
        return f"C-term({amino_acid_last}){index}"

    return ""

def process_ms_fragger_file(file: str):
    pep_seq_idx = -1
    pep_mod_seq_idx = -1
    prot_accession_idx = -1
    pep_exp_z_idx = -1

    exp_idx = []
    exp_names = []
    with open(input_dir + file, 'r') as f:
        while line := f.readline():
            if line.startswith('Peptide Sequence'):
                row = line.split('\t')
                for i, field in enumerate(row):
                    if field == "Protein ID":
                        prot_accession_idx = i
                    elif field == "Peptide Sequence":
                        pep_seq_idx = i
                    elif field == "Modified Sequence":
                        pep_mod_seq_idx = i
                    elif field == "Charges":
                        pep_exp_z_idx = i
                    elif "Intensity" in field:
                        if not "MaxLFQ Intensity" in field:
                            exp_idx.append(i)
                            exp_names.append(field.replace(" Intensity", ""))
            else:
                row = line.split('\t')
                nterm_cleav = ""
                cterm_cleav = ""
                all_mods = []
                mods_for_exp = {}

                if len(row) > 0:
                    if row[prot_accession_idx] in READER_CONFIG.ISOFORM_HELPER_DICT:
                        row[prot_accession_idx] = READER_CONFIG.ISOFORM_HELPER_DICT[row[prot_accession_idx]]
                    if row[prot_accession_idx] not in sorted_isoform_headers:
                        continue
                    else:
                        isoform, sequence, offset = get_accession(row[prot_accession_idx], row[pep_seq_idx])
                        
                        nterm_cleav = check_N_term_cleavage(sequence, row[pep_seq_idx], row[prot_accession_idx])
                        cterm_cleav = check_C_term_cleavage(sequence, row[pep_seq_idx], row[prot_accession_idx])
                        cleavage = ""
                        if nterm_cleav != "" and cterm_cleav != "":
                            cleavage = nterm_cleav + "; " + cterm_cleav
                        elif nterm_cleav != "":
                            cleavage = nterm_cleav
                        elif cterm_cleav != "":
                            cleavage = cterm_cleav



def process_ms_fragger_dir():
    for file in os.listdir(input_dir):
        if file.endswith(".tsv"):
            process_ms_fragger_file(file)
            #TODO: check if we should process multiple files for ms fragger
            break

process_ms_fragger_dir()