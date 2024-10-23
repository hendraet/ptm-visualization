import csv
import importlib
import os
import pandas as pd
from protein_sequencing import uniprot_align
from protein_sequencing.data_preprocessing import reader_helper
from typing import Tuple
import re

# TODO make changable from main input script

CONFIG = importlib.import_module('configs.default_config', 'configs')
READER_CONFIG = importlib.import_module('configs.reader_config', 'configs')

fasta_file = READER_CONFIG.FASTA_FILE
aligned_fasta_file = READER_CONFIG.ALIGNED_FASTA_FILE
input_file = READER_CONFIG.MS_FRAGGER_FILE

groups_df = pd.read_csv(f"{os.path.dirname(__file__)}/groups.csv")

sorted_isoform_headers = reader_helper.process_tau_file(fasta_file, aligned_fasta_file)

ms_fragger_mods = READER_CONFIG.MS_FRAGGER_MODS

def get_accession(accession: str, peptide: str) -> Tuple[str, str, int]:
    offset = 0
    sequence = None
    for header in sorted_isoform_headers:
        if peptide in header[1]:
            isoform = header[0]
            sequence = header[1]
            offset = sequence.index(peptide)
            break
    if sequence is not None:
        return isoform, sequence, offset
    else:
        raise ValueError(f"Peptide {peptide} with accession {accession} not found in fasta file")

def check_N_term_cleavage(peptide: str, accession: str) -> str:
    isoform, sequence, offset = get_accession(accession, peptide)
    if isoform in READER_CONFIG.ISOFORM_TRANSPOSE_DICT:
        isoform = READER_CONFIG.ISOFORM_TRANSPOSE_DICT[isoform]
    amino_acid_first = peptide[0]
    amino_acid_before = ""
    if offset > 0:
        amino_acid_before = sequence[offset - 1]
    if amino_acid_before != "K" and amino_acid_before != "R":
        return f"{amino_acid_first}@{offset+1}_{isoform}"

    return ""

def check_C_term_cleavage(peptide: str, accession: str) -> str:
    isoform, _, offset = get_accession(accession, peptide)
    if isoform in READER_CONFIG.ISOFORM_TRANSPOSE_DICT:
        isoform = READER_CONFIG.ISOFORM_TRANSPOSE_DICT[isoform]
    amino_acid_last = peptide[-1]
    if amino_acid_last not in ["K", "R"]:
        return f"{amino_acid_last}@{offset+len(peptide)}_{isoform}"

    return ""

def check_modification_present(mod_sequence: str) -> bool:
    relevant_mod_present = False
    for mod_idx in ms_fragger_mods:
        if mod_idx in mod_sequence:
            relevant_mod_present = True
            break
    return relevant_mod_present

def get_exact_indexes(mod_sequence: str) -> list:
    indexes = []
    current_index = 1
    inside_brackets = False
    for i, char in enumerate(mod_sequence):
        if char == '[':
            inside_brackets = True
        elif char == ']':
            inside_brackets = False
        elif not inside_brackets and char.isalpha():
            if i + 1 < len(mod_sequence) and mod_sequence[i + 1] == '[':
                indexes.append(current_index)
        if not inside_brackets:
            current_index += 1

    return indexes

def process_modifications(mod_sequence: str, offset: int, isoform: str) -> list:
    all_mods = []
    matches = re.findall(r'\[(\d+\.\d+?)\]', mod_sequence)
    mod_indexes = get_exact_indexes(mod_sequence)
    for i, match in enumerate(matches):
        if match in ms_fragger_mods and ms_fragger_mods[match] in CONFIG.MODIFICATIONS:
            modified_aa = mod_sequence[mod_sequence.index(match)-2]
            if modified_aa in CONFIG.EXCLUDED_MODIFICATIONS:
                if CONFIG.EXCLUDED_MODIFICATIONS[modified_aa] is not None and ms_fragger_mods[match] in CONFIG.EXCLUDED_MODIFICATIONS[modified_aa]:
                    continue
            offset = mod_indexes[i] + offset
            mod_string = f"{ms_fragger_mods[match]}({modified_aa})@{offset}_{isoform}"
            all_mods.append(mod_string)
    return all_mods

def write_results(all_mods, mods_for_exp, cleavages_with_ranges, cleavages_for_exp):
    with open(f"{CONFIG.OUTPUT_FOLDER}/result_ms_fragger_mods.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + all_mods)
        writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mods])
        writer.writerow(['', ''] + [reader_helper.extract_mod_location(mod) for mod in all_mods])
        for key, value in mods_for_exp.items():
            row = [1 if mod in value else 0 for mod in all_mods]
            group = groups_df.loc[groups_df['file_name'] == key]['group_name'].values[0]
            writer.writerow([key, group] + row)

    with open(f"{CONFIG.OUTPUT_FOLDER}/result_ms_fragger_cleavages.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + cleavages_with_ranges)
        writer.writerow(['', ''] + ['Non-Tryptic' for _ in cleavages_with_ranges])
        writer.writerow(['', ''] + [cleavage for cleavage in cleavages_with_ranges])
        ranges = reader_helper.parse_ranges(cleavages_with_ranges)
        for key, value in cleavages_for_exp.items():
            indexes = [reader_helper.extract_index(cleavage) for cleavage in value]
            row = reader_helper.cleavage_score(ranges, indexes)
            group = groups_df.loc[groups_df['file_name'] == key]['group_name'].values[0]
            writer.writerow([key, group] + row)

def process_ms_fragger_file(file: str):
    pep_seq_idx = -1
    pep_mod_seq_idx = -1
    prot_accession_idx = -1

    exp_idx = []
    exp_names = []

    all_mods = []
    mods_for_exp = {}
    all_cleavages = []
    cleavages_for_exp = {}
    for key in groups_df['file_name']:
        mods_for_exp[key] = []
        cleavages_for_exp[key] = []

    with open(file, 'r') as f:
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
                    elif "Intensity" in field:
                        if not "MaxLFQ Intensity" in field:
                            exp_idx.append(i)
                            exp_names.append(field.replace(" Intensity", "").strip())
            else:
                row = line.split('\t')
                row = [field.strip() for field in row]
                nterm_cleav = ""
                cterm_cleav = ""

                if len(row) > 0:
                    if row[prot_accession_idx] in READER_CONFIG.ISOFORM_HELPER_DICT:
                        row[prot_accession_idx] = READER_CONFIG.ISOFORM_HELPER_DICT[row[prot_accession_idx]]
                    try:
                        isoform, _, offset = get_accession(row[prot_accession_idx], row[pep_seq_idx])
                        if isoform in READER_CONFIG.ISOFORM_TRANSPOSE_DICT:
                            isoform = READER_CONFIG.ISOFORM_TRANSPOSE_DICT[isoform]
                    except ValueError:
                        continue                       
                        
                    nterm_cleav = check_N_term_cleavage(row[pep_seq_idx], row[prot_accession_idx])
                    cterm_cleav = check_C_term_cleavage(row[pep_seq_idx], row[prot_accession_idx])
                    cleavage = ""
                    if nterm_cleav != "" and cterm_cleav != "":
                        all_cleavages.append(nterm_cleav)
                        all_cleavages.append(cterm_cleav)
                        cleavage = nterm_cleav + "; " + cterm_cleav
                    elif nterm_cleav != "":
                        all_cleavages.append(nterm_cleav)
                        cleavage = nterm_cleav
                    elif cterm_cleav != "":
                        all_cleavages.append(cterm_cleav)
                        cleavage = cterm_cleav

                    if check_modification_present(row[pep_mod_seq_idx]):
                        mods_for_peptide = process_modifications(row[pep_mod_seq_idx], offset, isoform)
                        all_mods.extend(mods_for_peptide)
                        for i, idx in enumerate(exp_idx):
                            if row[idx] != "0.0":
                                mods_for_exp[exp_names[i]].extend(mods_for_peptide)

                                if cleavage != "":
                                    cleavages_for_exp[exp_names[i]].append(cleavage)

    all_mods = sorted(set(all_mods), key=reader_helper.extract_index)
    for key in mods_for_exp:
        mods_for_exp[key] = sorted(set(mods_for_exp[key]), key=reader_helper.extract_index)

    all_cleavages = sorted(set(all_cleavages), key=reader_helper.extract_cleavage_location)
    cleavages_with_ranges = reader_helper.extract_cleavages_ranges(all_cleavages)
    write_results(all_mods, mods_for_exp, cleavages_with_ranges, cleavages_for_exp)

uniprot_align.get_alignment(fasta_file)
process_ms_fragger_file(input_file)