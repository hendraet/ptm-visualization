import importlib
import os
import csv
import re
import pandas as pd
from protein_sequencing.data_preprocessing import reader_helper


CONFIG = importlib.import_module('configs.default_config', 'configs')
READER_CONFIG = importlib.import_module('configs.reader_config', 'configs')

fasta_file = READER_CONFIG.FASTA_FILE
aligned_fasta_file = READER_CONFIG.ALIGNED_FASTA_FILE
input_dir = READER_CONFIG.INPUT_DIR

current_dir = os.path.dirname(__file__)
groups_df = pd.read_csv(f"{current_dir}/groups.csv")

isoform_helper_dict = READER_CONFIG.isoform_helper_dict

def reformmods(mods, sites, short_sequence, variable_mods, isoform, sequence, aligned_sequence):
    modstrings = []
    sites = sites.split('.')[1]
    sequence_split_offset = 0
    for i, site in enumerate(sites):
        site = int(site)
        if site != 0:
            mod = variable_mods[site]
            aa = short_sequence[i]
            pos = sequence.index(short_sequence) + sequence_split_offset
            modstring = f"{mod}({aa})@{pos}_{isoform}"
            modstrings.append(modstring)
        sequence_split_offset += 1
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
                    if search_header in isoform_helper_dict:
                        search_header = isoform_helper_dict[search_header]
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
    with open(f"{CONFIG.OUTPUT_FOLDER}/result_mascot.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + all_mod_strings)
        writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mod_strings])
        writer.writerow(['', ''] + [reader_helper.extract_mod_location(mod) for mod in all_mod_strings])
        for file, mods in mod_strings_for_files.items():
            row = [1 if mod in mods else 0 for mod in all_mod_strings]
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

fasta_headers = reader_helper.process_tau_file(fasta_file, aligned_fasta_file)
process_mascot_dir(input_dir, fasta_headers)