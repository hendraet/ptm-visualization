import importlib
import os
import csv
import re
import pandas as pd
from protein_sequencing import exon_helper, uniprot_align
from protein_sequencing.data_preprocessing import reader_helper


CONFIG = importlib.import_module('configs.default_config', 'configs')
READER_CONFIG = importlib.import_module('configs.reader_config', 'configs')

fasta_file = READER_CONFIG.FASTA_FILE
aligned_fasta_file = READER_CONFIG.ALIGNED_FASTA_FILE
input_dir = READER_CONFIG.MASCOT_INPUT_DIR

groups_df = pd.read_csv(READER_CONFIG.GROUPS_CSV)
exon_found, exon_start_index, exon_end_index, exon_length, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length, exon_none_isoforms, max_sequence_length = exon_helper.retrieve_exon(fasta_file, CONFIG.MIN_EXON_LENGTH)

def reformmods(mods, sites, peptide, variable_mods, isoform, sequence, aligned_sequence):
    modstrings = []
    sites = sites.split('.')[1]
    for i, site in enumerate(sites):
        site = int(site)
        if site != 0:
            mod = variable_mods[site]
            aa = peptide[i]
            if CONFIG.INCLUDED_MODIFICATIONS.get(mod):
                if aa not in CONFIG.INCLUDED_MODIFICATIONS[mod]:
                    continue
                if aa == 'R' and mod == 'Deamidated':
                    mod = 'Citrullination'
            peptide_offset = sequence.index(peptide)+1
            missing_aa = 0
            if len(sequence) != len(aligned_sequence):
                missing_aa = reader_helper.count_missing_amino_acids(peptide[:i], aligned_sequence, peptide_offset, exon_start_index, exon_end_index)
            offset = reader_helper.calculate_exon_offset(peptide_offset + i + missing_aa, isoform, exon_found, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
            if aligned_sequence[offset-1] != aa:
                raise ValueError(f"AA don't match for {aa} for peptide {peptide} in sequence {sequence} with offset {offset}")
            iso = reader_helper.get_isoform_for_offset(isoform, offset, exon_start_index, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length)
            modstring = f"{mod}({aa})@{offset}_{iso}"
            modstrings.append(modstring)
    return modstrings

def process_mascot_file(file, fasta_headers):

    def replace_comma_in_quotes(line):
        pattern = r'\"(.*?)\"'
        def replace_comma(match):
            return match.group(0).replace(',', '-')
        modified_line = re.sub(pattern, replace_comma, line)
        return modified_line

    pep_seq_idx = -1
    pep_var_mod_idx = -1
    pep_var_mod_pos_idx = -1
    pep_accession_idx = -1

    header_length = -1

    varmods = False
    table = False
    emptyrowpars = False

    variable_mods = {}
    all_mod_strings = []
    with open(input_dir + file, 'r') as f:
        while line := f.readline():
            if line.startswith('\"Variable modifications'):
                varmods = True
            elif line.startswith('\"Protein hits'):
                table = True
            elif varmods == True and line == '\n':
                emptyrowpars = True
            elif table == True and line == '\n':
                emptyrowpars = True
            elif '--------------------------------------------------------' in line:
                varmods = False
                table = False
                emptyrowpars = False
            elif varmods == True and emptyrowpars == True:
                if not line.startswith('\"Identifier') and line != '\n':
                    mod_string = line.split(',')[1]
                    mod = mod_string.split(' ')[0][1:]
                    variable_mods[int(line[0])] = mod
            elif table == True and emptyrowpars == True:
                if line.startswith('prot_hit_num'):
                    headers = line.split(',')
                    header_length = len(headers)
                    for i, header in enumerate(headers):
                        if header == 'pep_seq':
                            pep_seq_idx = i
                        elif header == 'pep_var_mod':
                            pep_var_mod_idx = i
                        elif header == 'pep_var_mod_pos':
                            pep_var_mod_pos_idx = i
                        elif header == 'prot_acc':
                            pep_accession_idx = i
                else:
                    line = line.replace(', ', '-')
                    line = replace_comma_in_quotes(line)
                    row = line.split(',')
                    def replace_comma_in_quotes(line):
                        pattern = r'\"(.*?)\"'
                        def replace_comma(match):
                            return match.group(0).replace(',', '-')
                        modified_line = re.sub(pattern, replace_comma, line)
                        return modified_line
                    header_found = False
                    search_header = row[pep_accession_idx][1:-1]
                    if search_header in READER_CONFIG.ISOFORM_HELPER_DICT:
                        search_header = READER_CONFIG.ISOFORM_HELPER_DICT[search_header]
                    for fasta_header in fasta_headers:
                        if search_header in fasta_header[0]:
                            header_found = True
                            break
                    if not header_found:
                        continue
                    if len(row) != header_length:
                        pattern = r"Protein: \w{6,} has \d+ cached-but \d+ from the re-load*"
                        if re.search(pattern, line):
                            next_line = f.readline()
                            line = re.sub(pattern, '', line).strip()
                            line = line + next_line
                            if len(line.split(',')) != header_length:
                                print(f"Failed automatically fixing line: {line}")
                                continue
                            else:
                                row = line.split(',')
                        else:
                            print(f"Error with line: {line}")
                            continue
                    if row[pep_var_mod_idx] == '""':
                        continue
                    isoform = None
                    sequence = None
                    aligned_sequence = None
                    for fasta_header in fasta_headers:
                        if row[pep_seq_idx] in fasta_header[1]:
                                isoform = fasta_header[0]
                                sequence = fasta_header[1]
                                aligned_sequence = fasta_header[2]
                                break
                    if isoform is None:
                        continue
                    reform_mod_strings = reformmods(row[pep_var_mod_idx], row[pep_var_mod_pos_idx], row[pep_seq_idx], variable_mods, isoform, sequence, aligned_sequence)
                    all_mod_strings.extend(reform_mod_strings)
    return all_mod_strings

def process_results(all_mod_strings, mod_strings_for_files):    
    all_mod_strings = sorted(set(all_mod_strings), key=reader_helper.extract_index)
    all_mods = reader_helper.sort_by_index_and_exons(all_mod_strings)
    with open(f"{CONFIG.OUTPUT_FOLDER}/result_mascot.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Group'] + all_mods)
        writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mods])
        writer.writerow(['', ''] + [reader_helper.extract_mod_location(mod) for mod in all_mods])
        writer.writerow(['', ''] + [mod.split('_')[1] for mod in all_mods])
        for file, mods in mod_strings_for_files.items():
            row = [1 if mod in mods else 0 for mod in all_mods]
            group = groups_df.loc[groups_df['file_name'] == file]['group_name'].values[0]
            writer.writerow([file, group] + row)


def process_mascot_dir(input_dir, tau_headers):
    files = os.listdir(input_dir)
    all_mod_strings = []
    mod_strings_for_files = {}
    for file in files:
        result= process_mascot_file(file, tau_headers)
        all_mod_strings.extend(result)
        mod_strings_for_files[file] = result

    process_results(all_mod_strings, mod_strings_for_files)

uniprot_align.get_alignment(fasta_file)
fasta_headers = reader_helper.process_tau_file(fasta_file, aligned_fasta_file)

process_mascot_dir(input_dir, fasta_headers)