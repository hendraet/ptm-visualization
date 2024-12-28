"""MaxQuant preprocessor module. Extracts modifications and cleavages from MaxQuant output file."""
import importlib
import re
import pandas as pd
from protein_sequencing import exon_helper, uniprot_align
from protein_sequencing.data_preprocessing import preprocessor_helper

CONFIG = importlib.import_module('configs.default_config', 'configs')
PREPROCESSOR_CONFIG = importlib.import_module('configs.preprocessor_config', 'configs')

fasta_file = PREPROCESSOR_CONFIG.FASTA_FILE
aligned_fasta_file = PREPROCESSOR_CONFIG.ALIGNED_FASTA_FILE
input_file = PREPROCESSOR_CONFIG.MAX_QUANT_FILE

groups_df = pd.read_csv(PREPROCESSOR_CONFIG.GROUPS_CSV)
exon_found, exon_start_index, exon_end_index, exon_length, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length, exon_none_isoforms, max_sequence_length = exon_helper.retrieve_exon(fasta_file, CONFIG.MIN_EXON_LENGTH)

def get_exact_indexes(mod_sequence: str) -> list:
    """Get exact indexes of the modifications in the sequence."""
    indexes = []
    current_index = 0
    inside_brackets = False
    inside_second_brackets = False
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
        elif not inside_second_brackets and char.isalpha() or char == '_':
            if i + 1 < len(mod_sequence) and mod_sequence[i + 1] == '(':
                indexes.append(current_index)
        if not inside_brackets:
            current_index += 1

    return indexes

def reformat_mod(modified_peptide: str, peptide: str, peptide_offset: int, sequence: str, isoform: str, aligned_sequence: str) -> list[str]:
    """Reformat the modification string."""
    mod_strings = []

    pattern = r"\(([^()]*)\)"
    matches = re.findall(pattern, modified_peptide)
    indexes = get_exact_indexes(modified_peptide)
    counter = 0
    for mod_position in matches:
        if mod_position == "Protein N-term":
            aa = sequence[peptide_offset-1]
            aa_offset = 0
        elif mod_position == "Protein C-term":
            aa = peptide[-1]
            aa_offset = len(peptide)
        else:
            aa = peptide[indexes[counter]-1]
            aa_offset = indexes[counter]
        if mod_position == 'ph':
            mod_type = 'Phospho'
        elif mod_position == 'ac':
            mod_type = 'Acetyl'
        elif mod_position == 'gg':
            mod_type = 'GG'
        elif mod_position == 'me':
            mod_type = 'Methyl'
        elif mod_position == 'ci':
            mod_type = 'Citrullination'
        elif mod_position == 'de':
            mod_type = 'Deamidated'
        else:
            mod_type = mod_position.split(' ')[0]
        if CONFIG.INCLUDED_MODIFICATIONS.get(mod_type):
            if aa not in CONFIG.INCLUDED_MODIFICATIONS[mod_type]:
                continue
            if aa == 'R' and mod_type == 'Deamidated':
                mod_type = 'Citrullination'
        else:
            continue
        missing_aa = 0
        if len(sequence) != len(aligned_sequence):
            missing_aa = preprocessor_helper.count_missing_amino_acids(peptide[:aa_offset], aligned_sequence, peptide_offset, exon_start_index, exon_end_index)
        offset = preprocessor_helper.calculate_exon_offset(aa_offset+peptide_offset+missing_aa, isoform, exon_found, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
        if aligned_sequence[offset-1] != aa:
            raise ValueError(f"AA don't match for {aa} for peptide {peptide} in sequence {sequence} with offset {offset}")
        iso = preprocessor_helper.get_isoform_for_offset(isoform, offset, exon_start_index, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length)
        mod_strings.append(f"{mod_type}({aa})@{offset}_{iso}")
        counter += 1
    return mod_strings

def process_max_quant_file(evidence_file: str):
    """Process MaxQuant file."""
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

    with open(evidence_file, 'r', encoding="utf-8") as f:
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
                if fields[prot_accession_idx] in PREPROCESSOR_CONFIG.ISOFORM_HELPER_DICT:
                    fields[prot_accession_idx] = PREPROCESSOR_CONFIG.ISOFORM_HELPER_DICT[fields[prot_accession_idx]]
                try:
                    isoform, sequence, peptide_offset, aligned_sequence = preprocessor_helper.get_accession(fields[prot_accession_idx], fields[pep_seq_idx], sorted_isoform_headers)
                except ValueError:
                    continue

                cleavage = preprocessor_helper.check_N_term_cleavage(fields[pep_seq_idx], fields[prot_accession_idx], sorted_isoform_headers, exon_found, exon_start_index, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
                if cleavage != "":
                    all_cleavages.append(cleavage)
                    if fields[exp_idx] not in cleavages_for_exp:
                        cleavages_for_exp[fields[exp_idx]] = []
                    cleavages_for_exp[fields[exp_idx]].append(cleavage)
                cleavage = preprocessor_helper.check_C_term_cleavage(fields[pep_seq_idx], fields[prot_accession_idx], sorted_isoform_headers, exon_found, exon_start_index, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
                if cleavage != "":
                    all_cleavages.append(cleavage)
                    if fields[exp_idx] not in cleavages_for_exp:
                        cleavages_for_exp[fields[exp_idx]] = []
                    cleavages_for_exp[fields[exp_idx]].append(cleavage)

                if float(fields[pep_score_idx]) < PREPROCESSOR_CONFIG.THRESHOLD:
                    if fields[mods_idx] != "Unmodified":
                        mods = reformat_mod(fields[pep_mod_seq_idx], fields[pep_seq_idx], peptide_offset, sequence, isoform, aligned_sequence)
                        all_mods.extend(mods)
                        if fields[exp_idx] in mods_for_exp:
                            mods_for_exp[fields[exp_idx]].extend(mods)

    all_mods = sorted(set(all_mods), key=preprocessor_helper.extract_index)
    all_mods = preprocessor_helper.sort_by_index_and_exons(all_mods)
    for key in mods_for_exp:
        mods_for_exp[key] = sorted(set(mods_for_exp[key]), key=preprocessor_helper.extract_index)
        mods_for_exp[key] = preprocessor_helper.sort_by_index_and_exons(mods_for_exp[key])

    all_cleavages = sorted(set(all_cleavages), key=preprocessor_helper.extract_cleavage_location)
    all_cleavages = preprocessor_helper.sort_by_index_and_exons(all_cleavages)
    cleavages_with_ranges = preprocessor_helper.extract_cleavages_ranges(all_cleavages)
    preprocessor_helper.write_results(all_mods, mods_for_exp, cleavages_with_ranges, cleavages_for_exp, CONFIG.OUTPUT_FOLDER, groups_df)

uniprot_align.get_alignment(fasta_file)
sorted_isoform_headers = preprocessor_helper.process_tau_file(fasta_file, aligned_fasta_file)

process_max_quant_file(input_file)
