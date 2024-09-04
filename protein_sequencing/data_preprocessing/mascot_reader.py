import importlib
import os
import csv

fasta_file = 'data/test_data/uniprot/test.fasta'
aligned_fasta_file = 'data/test_data/uniprot/test_aligned.fasta'
input_dir = 'data/test_data/mascot/'

# TODO make changable from main input script
CONFIG = importlib.import_module('configs.default_config', 'configs')

def process_tau_file(fasta_file, aligned_fasta_file):
    headers = []
    aligned_sequences = {}
    with open(aligned_fasta_file, 'r') as file:
        lines = file.readlines()
        header = ''
        aligned_seq = ''
        for line in lines:
            if line.startswith('>'):
                if header != '':
                    aligned_sequences[header] = aligned_seq
                header = line.split('|')[1]
                aligned_seq = ''
            else:
                aligned_seq += line.strip()
        if header:
            aligned_sequences[header] = aligned_seq

    with open(fasta_file, 'r') as file:
        lines = file.readlines()
        seq = ''
        header = ''
        for line in lines:
            if line[0] == '>':
                if header != '':
                    headers.append((header, seq, aligned_sequences[header]))
                    seq = ''
                    header = ''
                header = line.split('|')[1]
            else:
                seq += line.strip()
        headers.append((header, seq, aligned_sequences[header]))
    sorted_headers = sorted(headers, key=lambda x: -len(x[1]))
    return sorted_headers

def reformmods(mods, sites, short_sequence, pep_start, variable_mods, isoform, aligned_sequence):
    # get new start index for aligned isoform
    aligned_pep_start = 0
    break_counter = int(pep_start)
    for i in range(len(aligned_sequence)):
        aligned_pep_start += 1
        if aligned_sequence[i] != '-':
            break_counter -= 1
            if break_counter == 0:
                break

    modstrings = []

    # TODO do we acutally need this?
    mods = mods[1:-1].split(';')
    single_mods = []
    for mod in mods:
        if mod[0].isdigit():
            single_mod = mod[2:]
            for i in range(int(mod[0])):
                single_mods.append(single_mod)
        else:
            single_mods.append(mod)
    mods = [mod.strip() for mod in single_mods]

    # "2 Acetyl (K); GG (K); Oxidation (M); 2 Phospho (ST)"
    # 0.0400030011042.0
    #   SSQLQMGQKKNSK
    # abcdefabc
    # ABC----------DEFABC------
    # pep_start = 6
    # sequence = "abc"
    # mods = "Oxidation (B)"
    # sites = "0.030.0"
    sites = sites.split('.')[1]
    sequence_split_offset = 0
    for i, site in enumerate(sites):
        site = int(site)
        if site != 0:
            mod = variable_mods[site]
            aa = short_sequence[i]
            # -1 because the sequence is 1 indexed
            while aligned_sequence[(aligned_pep_start-1) + sequence_split_offset] == '-':
                sequence_split_offset += 1
            siteidx = aligned_pep_start + sequence_split_offset
                
            modstring = f"{mod}({aa})@{siteidx}_{isoform}"
            modstrings.append(modstring)
        sequence_split_offset += 1
    return modstrings

def process_mascot_file(file, fasta_headers):
    # TODO do we need all these variables?
    prot_acc_idx = -1
    pep_seq_idx = -1
    pep_var_mod_idx = -1
    pep_var_mod_pos_idx = -1
    pep_score_idx = -1
    pep_delta_idx = -1
    pep_exp_mr_idx = -1
    pep_exp_mz_idx = -1
    pep_calc_mr_idx = -1
    pep_exp_z_idx = -1
    pep_start_idx = -1
    pep_end_idx = -1

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
                    for i, header in enumerate(headers):
                        if header == 'prot_acc':
                            prot_acc_idx = i
                        elif header == 'pep_seq':
                            pep_seq_idx = i
                        elif header == 'pep_var_mod':
                            pep_var_mod_idx = i
                        elif header == 'pep_var_mod_pos':
                            pep_var_mod_pos_idx = i
                        elif header == 'pep_score':
                            pep_score_idx = i
                        elif header == 'pep_delta':
                            pep_delta_idx = i
                        elif header == 'pep_exp_mr':
                            pep_exp_mr_idx = i
                        elif header == 'pep_exp_mz':
                            pep_exp_mz_idx = i
                        elif header == 'pep_calc_mr':
                            pep_calc_mr_idx = i
                        elif header == 'pep_exp_z':
                            pep_exp_z_idx = i
                        elif header == 'pep_start':
                            pep_start_idx = i
                        elif header == 'pep_end':
                            pep_end_idx = i
                else:
                    line = line.replace(', ', '-')
                    row = line.split(',')
                    if row[pep_var_mod_idx] == '""':
                        continue
                    isoform = None
                    sequence = None
                    aligned_sequence = None
                    for fasta_header in fasta_headers:
                        if row[pep_seq_idx] in fasta_header[1]:
                            start = int(row[pep_start_idx])
                            end = int(row[pep_end_idx])
                            wrong_sequence = False
                            for i in range(end - start + 1):
                                if fasta_header[1][start + i -1] != row[pep_seq_idx][i]:
                                    wrong_sequence = True
                                    break
                            if not wrong_sequence:
                                isoform = fasta_header[0]
                                sequence = fasta_header[1]
                                aligned_sequence = fasta_header[2]
                                break
                    if isoform is None:
                        print(f"Could not find isoform for sequence {row[pep_seq_idx]}")
                        continue
                    reform_mod_strings = reformmods(row[pep_var_mod_idx], row[pep_var_mod_pos_idx], row[pep_seq_idx], row[pep_start_idx], variable_mods, isoform, aligned_sequence)
                    all_mod_strings.extend(reform_mod_strings)
    return all_mod_strings

def process_results(all_mod_strings, mod_strings_for_files):
    def extract_index(mod_string):
        return int(mod_string.split('@')[1].split('_')[0])
    
    def extract_location(mod_string):
        return mod_string.split('(')[1].split(')')[0]+mod_string.split('@')[1].split('_')[0]
    
    all_mod_strings = sorted(set(all_mod_strings), key=extract_index)
    with open(f"{CONFIG.OUTPUT_FOLDER}/result_mascot.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + all_mod_strings)
        writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mod_strings])
        writer.writerow(['', ''] + [extract_location(mod) for mod in all_mod_strings])
        writer.writerow(['', ''] + [mod.split('_')[1] for mod in all_mod_strings])
        for file, mods in mod_strings_for_files.items():
            row = [1 if mod in mods else 0 for mod in all_mod_strings]
            writer.writerow([file, ''] + row)

def process_mascot_dir(input_dir, tau_headers):
    files = os.listdir(input_dir)
    all_mod_strings = []
    mod_strings_for_files = {}
    for file in files:
        result= process_mascot_file(file, tau_headers)
        all_mod_strings.extend(result)
        mod_strings_for_files[file] = result

    process_results(all_mod_strings, mod_strings_for_files)

fasta_headers = process_tau_file(fasta_file, aligned_fasta_file)
process_mascot_dir(input_dir, fasta_headers)