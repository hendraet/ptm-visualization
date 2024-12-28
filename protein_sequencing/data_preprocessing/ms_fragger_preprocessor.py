"""MS Fragger Preprocessor Module. Extracts modifications and cleavages from MS Fragger output file."""
import csv
import importlib
import re
import pandas as pd
from protein_sequencing import exon_helper, uniprot_align
from protein_sequencing.data_preprocessing import preprocessor_helper

CONFIG = importlib.import_module('configs.default_config', 'configs')
PREPROCESSOR_CONFIG = importlib.import_module('configs.preprocessor_config', 'configs')

fasta_file = PREPROCESSOR_CONFIG.FASTA_FILE
aligned_fasta_file = PREPROCESSOR_CONFIG.ALIGNED_FASTA_FILE
input_file = PREPROCESSOR_CONFIG.MS_FRAGGER_FILE

groups_df = pd.read_csv(PREPROCESSOR_CONFIG.GROUPS_CSV)
exon_found, exon_start_index, exon_end_index, exon_length, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length, exon_none_isoforms, max_sequence_length = exon_helper.retrieve_exon(fasta_file, CONFIG.MIN_EXON_LENGTH)

def check_modification_present(mod_sequence: str) -> bool:
    """Check if relevant modification is present in the sequence."""
    relevant_mod_present = False
    for mod_idx in PREPROCESSOR_CONFIG.MS_FRAGGER_MODS:
        if mod_idx in mod_sequence:
            relevant_mod_present = True
            break
    return relevant_mod_present

def get_exact_indexes(mod_sequence: str) -> list:
    """Get exact indexes of the modifications in the sequence."""
    indexes = []
    current_index = 1
    inside_brackets = False
    for i, char in enumerate(mod_sequence):
        if char == '[':
            inside_brackets = True
        elif char == ']':
            inside_brackets = False
            continue
        elif not inside_brackets and char.isalpha():
            if i + 1 < len(mod_sequence) and mod_sequence[i + 1] == '[':
                indexes.append(current_index)
        if not inside_brackets:
            current_index += 1

    return indexes

def process_modifications(mod_sequence: str, peptide_offset: int, isoform: str, sequence: str, aligned_sequence: str) -> list:
    """Process modifications in the sequence."""
    all_mods = []
    matches = re.findall(r'\[(\d+\.\d+?)\]', mod_sequence)
    peptide = re.sub(r'\[(\d+\.\d+?)\]', '', mod_sequence)
    aa_offsets = get_exact_indexes(mod_sequence)
    for i, match in enumerate(matches):
        if match in PREPROCESSOR_CONFIG.MS_FRAGGER_MODS and PREPROCESSOR_CONFIG.MS_FRAGGER_MODS[match] in CONFIG.INCLUDED_MODIFICATIONS:
            modified_aa = peptide[aa_offsets[i]-1]
            if CONFIG.INCLUDED_MODIFICATIONS.get(match):
                if modified_aa not in CONFIG.INCLUDED_MODIFICATIONS[match]:
                    continue
                if modified_aa == 'R' and match == 'Deamidated':
                    match = 'Citrullination'
            missing_aa = 0
            if len(sequence) != len(aligned_sequence):
                missing_aa = preprocessor_helper.count_missing_amino_acids(peptide[:aa_offsets[i]], aligned_sequence, peptide_offset, exon_start_index, exon_end_index)
            offset = preprocessor_helper.calculate_exon_offset(aa_offsets[i] + peptide_offset+missing_aa, isoform, exon_found, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
            if aligned_sequence[offset-1] != modified_aa:
                raise ValueError(f"AA don't match for {modified_aa} for peptide {peptide} in sequence {sequence} with offset {offset}")
            iso = preprocessor_helper.get_isoform_for_offset(isoform, offset, exon_start_index, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length)
            mod_string = f"{PREPROCESSOR_CONFIG.MS_FRAGGER_MODS[match]}({modified_aa})@{offset}_{iso}"
            all_mods.append(mod_string)
    return all_mods

def process_ms_fragger_file(file: str):
    """Process MS Fragger output file."""
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

    with open(file, 'r', encoding="utf-8") as f:
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
                    if row[prot_accession_idx] in PREPROCESSOR_CONFIG.ISOFORM_HELPER_DICT:
                        row[prot_accession_idx] = PREPROCESSOR_CONFIG.ISOFORM_HELPER_DICT[row[prot_accession_idx]]
                    try:
                        isoform, sequence, offset, aligned_sequence = preprocessor_helper.get_accession(row[prot_accession_idx], row[pep_seq_idx], sorted_isoform_headers)
                    except ValueError:
                        continue

                    nterm_cleav = preprocessor_helper.check_N_term_cleavage(row[pep_seq_idx], row[prot_accession_idx], sorted_isoform_headers, exon_found, exon_start_index, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
                    cterm_cleav = preprocessor_helper.check_C_term_cleavage(row[pep_seq_idx], row[prot_accession_idx], sorted_isoform_headers, exon_found, exon_start_index, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
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
                        mods_for_peptide = process_modifications(row[pep_mod_seq_idx], offset, isoform, sequence, aligned_sequence)
                        all_mods.extend(mods_for_peptide)
                        for i, idx in enumerate(exp_idx):
                            if row[idx] != "0.0":
                                mods_for_exp[exp_names[i]].extend(mods_for_peptide)

                                if cleavage != "":
                                    cleavages_for_exp[exp_names[i]].append(cleavage)

    all_mods = sorted(set(all_mods), key=preprocessor_helper.extract_index)
    all_mods = preprocessor_helper.sort_by_index_and_exons(all_mods)
    for key in mods_for_exp:
        mods_for_exp[key] = sorted(set(mods_for_exp[key]), key=preprocessor_helper.extract_index)

    all_cleavages = sorted(set(all_cleavages), key=preprocessor_helper.extract_cleavage_location)
    all_cleavages = preprocessor_helper.sort_by_index_and_exons(all_cleavages)
    cleavages_with_ranges = preprocessor_helper.extract_cleavages_ranges(all_cleavages)
    preprocessor_helper.write_results(all_mods, mods_for_exp, cleavages_with_ranges, cleavages_for_exp, CONFIG.OUTPUT_FOLDER, groups_df)

uniprot_align.get_alignment(fasta_file)
sorted_isoform_headers = preprocessor_helper.process_tau_file(fasta_file, aligned_fasta_file)

process_ms_fragger_file(input_file)
