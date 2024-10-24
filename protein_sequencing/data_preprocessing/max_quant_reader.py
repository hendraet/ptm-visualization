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
input_file = READER_CONFIG.MAX_QUANT_FILE

groups_df = pd.read_csv(f"{os.path.dirname(__file__)}/groups_max_quant.csv")

def get_accession(accession: str, peptide: str) -> Tuple[str, str, int, str]:
    offset = 0
    sequence = None
    for header in sorted_isoform_headers:
        if peptide in header[1]:
            isoform = header[0]
            sequence = header[1]
            offset = sequence.index(peptide)
            # in case we didn't match the longest sequence we have to adapt the offset
            offset += sequence[:offset].count('-')
            break
    if sequence is not None:
        return isoform, sequence, offset, header[2]
    else:
        raise ValueError(f"Peptide {peptide} with accession {accession} not found in fasta file")

def count_missing_amino_acids(peptide: str, aligned_sequence: str, offset: int) -> int:
    missing = 0
    j = 0
    for i in range(offset, len(aligned_sequence)):
        if aligned_sequence[i] == '-':
            missing += 1
        elif peptide[j] == aligned_sequence[i]:
            j+=1
        if j == len(peptide):
            break
    return missing

def check_N_term_cleavage(peptide: str, accession: str) -> str:
    isoform, sequence, offset, _ = get_accession(accession, peptide)
    amino_acid_first = peptide[0]
    amino_acid_before = ""
    if offset > 0:
        amino_acid_before = sequence[offset - 1]
    if amino_acid_before != "K" and amino_acid_before != "R":
        return f"{amino_acid_first}@{offset+1}_{isoform}"

    return ""

def check_C_term_cleavage(peptide: str, accession: str) -> str:
    isoform, sequence, offset, aligned_sequence = get_accession(accession, peptide)
    missing_aa = 0
    if len(sequence) != len(aligned_sequence):
        missing_aa = count_missing_amino_acids(peptide, aligned_sequence, offset)
    amino_acid_last = peptide[-1]
    if amino_acid_last not in ["K", "R"]:
        return f"{amino_acid_last}@{offset+len(peptide)+missing_aa}_{isoform}"

    return ""

def get_exact_indexes(mod_sequence: str) -> list:
    indexes = []
    current_index = 1
    inside_brackets = False
    inside_second_brackets = False
    mod_sequence = mod_sequence.replace("_", "")
    for i, char in enumerate(mod_sequence):
        if char == '(':
            if inside_brackets:
                inside_second_brackets = True
            else:
                inside_brackets = True
        elif char == ')':
            if inside_second_brackets:
                inside_second_brackets = False
            else:
                inside_brackets = False
                continue
        elif not inside_second_brackets and char.isalpha():
            if i + 1 < len(mod_sequence) and mod_sequence[i + 1] == '(':
                indexes.append(current_index)
        if not inside_brackets:
            current_index += 1

    return indexes

def reformat_mod(modified_peptide: str, peptide: str, peptide_offset: int, sequence: str, isoform: str, aligned_sequence: str) -> list[str]:
    mod_strings = []
    
    pattern = r"\((\w+)\s*\(([^)]+)\)\)"
    matches = re.findall(pattern, modified_peptide)
    indexes = get_exact_indexes(modified_peptide)
    counter = 0
    for mod_type, mod_position in matches:
        if mod_position == "Protein N-term":
            mod_location = sequence[peptide_offset-1]
            aa_offset = 0
        elif mod_position == "Protein C-term":
            mod_location = peptide[-1]
            aa_offset = len(peptide)
        else:
            mod_location = peptide[indexes[counter]-1]
            aa_offset = indexes[counter]

        if mod_location in CONFIG.EXCLUDED_MODIFICATIONS:
            if CONFIG.EXCLUDED_MODIFICATIONS[mod_location] is not None and mod_type in CONFIG.EXCLUDED_MODIFICATIONS[mod_location]:
                continue

        if sequence[aa_offset+peptide_offset-1] != mod_location:
            raise ValueError(f"AA don't match for {mod_location} for peptide {peptide} in sequence {sequence} with offset {peptide_offset+aa_offset}")
        
        missing_aa = 0
        if len(sequence) != len(aligned_sequence):
            missing_aa = count_missing_amino_acids(peptide[:aa_offset], aligned_sequence, peptide_offset)
        mod_strings.append(f"{mod_type}({mod_location})@{aa_offset+peptide_offset+missing_aa}_{isoform}")
        counter += 1
    return mod_strings

def process_max_quant_file(input_file: str):
    pep_seq_idx = -1
    pep_mod_seq_idx = -1
    prot_accession_idx = -1
    mods_idx = -1
    exp_idx = -1
    pep_score_idx = -1

    all_mods = []
    mods_for_exp = {}
    all_cleavages = []
    cleavages_for_exp = {}

    for key in groups_df['file_name']:
        mods_for_exp[key] = []
        cleavages_for_exp[key] = []

    with open(input_file, 'r') as f:
        while line := f.readline():
            if line.startswith("Sequence"):
                header = line.split("\t")
                for i, field in enumerate(header):
                    if field == "Sequence":
                        pep_seq_idx = i
                    elif field == "Modified sequence":
                        pep_mod_seq_idx = i
                    elif field == "Modifications":
                        mods_idx = i
                    elif field == "Proteins":
                        prot_accession_idx = i
                    elif field == "PEP":
                        pep_score_idx = i
                    elif field.startswith("Experiment"):
                        exp_idx = i
            else:
                fields = line.split("\t")
                if fields[prot_accession_idx] in READER_CONFIG.ISOFORM_HELPER_DICT:
                    fields[prot_accession_idx] = READER_CONFIG.ISOFORM_HELPER_DICT[fields[prot_accession_idx]]
                try:
                    isoform, sequence, offset, aligned_sequence = get_accession(fields[prot_accession_idx], fields[pep_seq_idx])
                except ValueError:
                    continue

                cleavage = check_N_term_cleavage(fields[pep_seq_idx], fields[prot_accession_idx])
                if cleavage != "":
                    all_cleavages.append(cleavage)
                    if fields[exp_idx] not in cleavages_for_exp:
                        cleavages_for_exp[fields[exp_idx]] = []
                    cleavages_for_exp[fields[exp_idx]].append(cleavage)
                cleavage = check_C_term_cleavage(fields[pep_seq_idx], fields[prot_accession_idx])
                if cleavage != "":
                    all_cleavages.append(cleavage)
                    cleavages_for_exp[fields[exp_idx]].append(cleavage)
                
                if float(fields[pep_score_idx]) < READER_CONFIG.THRESHOLD:
                    if fields[mods_idx] != "Unmodified":
                        mods = reformat_mod(fields[pep_mod_seq_idx], fields[pep_seq_idx], offset, sequence, isoform, aligned_sequence)
                        all_mods.extend(mods)
                        mods_for_exp[fields[exp_idx]].extend(mods)

    all_mods = sorted(set(all_mods), key=reader_helper.extract_index)
    for key in mods_for_exp:
        mods_for_exp[key] = sorted(set(mods_for_exp[key]), key=reader_helper.extract_index)

    all_cleavages = sorted(set(all_cleavages), key=reader_helper.extract_cleavage_location)
    cleavages_with_ranges = reader_helper.extract_cleavages_ranges(all_cleavages)
    write_results(all_mods, mods_for_exp, cleavages_with_ranges, cleavages_for_exp)
                
def write_results(all_mods, mods_for_exp, cleavages_with_ranges, cleavages_for_exp):
    with open(f"{CONFIG.OUTPUT_FOLDER}/result_max_quant_mods.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + all_mods)
        writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mods])
        writer.writerow(['', ''] + [reader_helper.extract_mod_location(mod) for mod in all_mods])
        writer.writerow(['', ''] + [mod.split('_')[1] for mod in all_mods])
        for key, value in mods_for_exp.items():
            row = [1 if mod in value else 0 for mod in all_mods]
            group = groups_df.loc[groups_df['file_name'] == key]['group_name'].values[0]
            writer.writerow([key, group] + row)

    with open(f"{CONFIG.OUTPUT_FOLDER}/result_max_quant_cleavages.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + cleavages_with_ranges)
        writer.writerow(['', ''] + ['Non-Tryptic' for _ in cleavages_with_ranges])
        writer.writerow(['', ''] + [cleavage.split('_')[0] for cleavage in cleavages_with_ranges])
        writer.writerow(['', ''] + [cleavage.split('_')[1] for cleavage in cleavages_with_ranges])
        ranges = reader_helper.parse_ranges(cleavages_with_ranges)
        for key, value in cleavages_for_exp.items():
            indexes = [reader_helper.extract_index(cleavage) for cleavage in value]
            row = reader_helper.cleavage_score(ranges, indexes)
            group = groups_df.loc[groups_df['file_name'] == key]['group_name'].values[0]
            writer.writerow([key, group] + row)           
        
uniprot_align.get_alignment(fasta_file)
sorted_isoform_headers = reader_helper.process_tau_file(fasta_file, aligned_fasta_file)

process_max_quant_file(input_file)