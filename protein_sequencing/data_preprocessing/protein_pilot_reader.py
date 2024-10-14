import importlib
import os
import csv
import re
import pandas as pd
from python_calamine import CalamineWorkbook
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

def extract_confidence_score(calamine_sheet):
    global_column = READER_CONFIG.FDR_GLOBAL == 'global'
    threshold = READER_CONFIG.CONFIDENCE_THRESHOLD
    fdr_column, threshold_column = None, None
    if global_column:
        fdr_column = 'Fit Global FDR'
    else:
        fdr_column = 'Fit Local FDR'
    threshold_column = 'Fit Confidence Threshold'

    rows = iter(calamine_sheet.to_python())
    for row in rows:
        if threshold_column in row:
            threshold_index = row.index(threshold_column)
            fdr_index = row.index(fdr_column)
        else:
            if row[fdr_index] >= threshold:
                return row[threshold_index]
    raise Exception("No confidence score found.")

def split_mod(mod, seq):
    if '(' not in mod:
        mod_name = mod.split('@')[0]
        mod_pos = mod.split('@')[1]
        amino_acid = ''
    else:
        mod_name = mod.split('(')[0]
        amino_acid = mod.split('(')[1][0]
        mod_pos = mod.split('@')[1]
    if mod_pos == 'N-term':
        mod_pos = 0
        amino_acid = seq[0]
    if mod_pos == 'C-term':
        mod_pos = len(seq)+1
        amino_acid = seq[len(seq)-1]
    return mod_name.strip(), amino_acid.strip(), int(mod_pos)

def get_accession(row, accession_index, seq_index) -> Tuple[str, int] | Tuple[None, None]:
    search_header = row[accession_index].split('|')[1]
    if search_header in READER_CONFIG.ISOFORM_HELPER_DICT:
        search_header = READER_CONFIG.ISOFORM_HELPER_DICT[search_header]
    if search_header not in [header[0] for header in sorted_isoform_headers]:
        return None, None
    index_offset = None
    for header in sorted_isoform_headers:
        if row[seq_index] in header[1]:
            isoform = header[0].strip()
            index_offset = header[1].index(row[seq_index])
            break
    if index_offset is not None:
        return isoform, index_offset
    else:
        return None, None

def extract_mods_from_rows(rows, protein_mod_index, mod_index, seq_index, accession_index) -> list:
    mods = []
    for row in rows:
        isoform, index_offset = get_accession(row, accession_index, seq_index)
        if isoform is None:
            continue

        if index_offset is None:
            print(f"Error: sequence {row[seq_index]} with assumed accession: {row[accession_index]} not found in fasta file.")
            continue
        relevant_mods = row[protein_mod_index].split(';')
        all_mods = row[mod_index].split(';')
        seq = row[seq_index]
        for rel_mod in relevant_mods:
            matched_mod = None
            rel_mod_name, rel_amino_acid, rel_mod_pos = split_mod(rel_mod, seq)
            if rel_mod_name not in CONFIG.MODIFICATIONS:
                continue
            for mod in all_mods:
                mod_name, amino_acid, mod_pos = split_mod(mod, seq)
                if rel_mod_name == mod_name and rel_amino_acid == amino_acid:
                    matched_mod = mod
                    break
            if matched_mod is not None:
                mod_name, amino_acid, mod_pos = split_mod(matched_mod, seq)
                if CONFIG.EXCLUDED_MODIFICATIONS.get(amino_acid) is not None and mod_name in CONFIG.EXCLUDED_MODIFICATIONS[amino_acid]:
                    continue
                if amino_acid != seq[int(mod_pos)-1]:
                    print(f"Error: amino acid mismatch at position {mod_pos} for sequence {seq} and mod {matched_mod}")
                modstring = f"{mod_name}({amino_acid})@{mod_pos+index_offset}"
                mods.append(modstring)

    return mods

def extract_cleavages_from_rows(rows, cleavage_index, seq_index, accession_index) -> list:
    cleavages = []
    for row in rows:
        isoform, index_offset = get_accession(row, accession_index, seq_index)
        if index_offset is None:
            continue
        cleaved_sites = row[cleavage_index].split(';')
        for site in cleaved_sites:
            if not 'cleaved' in site:
                continue
            
            amino_acid = site.split('@')[0][-1]
            if 'N-term' in site:
                site_index = 0
            elif 'C-term' in site:
                site_index = len(row[seq_index])+1
            
            cleavages.append(f"{amino_acid}@{index_offset+site_index}")

    return cleavages

def extract_data_with_threshold(calamine_sheet, threshold):
    rows = iter(calamine_sheet.to_python())
    relevant_mod_rows = []
    relavent_cleavage_rows = []
    for row in rows:
        if 'ProteinModifications' in row and 'Conf' in row:
            protein_mod_index = row.index('ProteinModifications')
            mod_index = row.index('Modifications')
            seq_index = row.index('Sequence')
            conf_index = row.index('Conf')
            cleavage_index = row.index('Cleavages')
            accession_index = row.index('Accessions')
        else:
            if row[conf_index] >= threshold*100:
                if row[protein_mod_index] is not None and row[protein_mod_index] != '':
                    relevant_mod_rows.append(row)
                if row[cleavage_index] is not None and row[cleavage_index] != '' and 'cleaved' in row[cleavage_index]:
                    relavent_cleavage_rows.append(row)

    mods = extract_mods_from_rows(relevant_mod_rows, protein_mod_index, mod_index, seq_index, accession_index)
    cleavages = extract_cleavages_from_rows(relavent_cleavage_rows, cleavage_index, seq_index, accession_index)

    return mods, cleavages

def process_protein_pilot_xlsx_file(file) -> Tuple[list, list]:
    workbook = CalamineWorkbook.from_path(file)

    distinct_peptide_level_data = workbook.get_sheet_by_name('Distinct Peptide Level Data')
    fdr_threshold = extract_confidence_score(distinct_peptide_level_data)

    peptide_summary = workbook.get_sheet_by_name('Peptide Summary')
    mods, cleavages = extract_data_with_threshold(peptide_summary, fdr_threshold)

    return mods, cleavages

def parse_ranges(ranges_list):
    ranges = []
    for part in ranges_list:
        if '-' in part:
            start, end = map(int, part.split('-'))
            ranges.append(list(range(start, end + 1)))
        else:
            ranges.append([int(part)])
    return ranges

def cleavage_score(ranges, cleavages):
    results = []
    for r in ranges:
        range_length = len(r)
        cleavage_hits = len([x for x in r if x in cleavages])
        
        if cleavage_hits == 0:
            results.append(0)
        elif cleavage_hits == range_length:
            results.append(1)
        else:
            results.append(cleavage_hits / range_length)
    return results

def process_protein_pilot_dir():
    all_mods = []
    all_cleavages = []
    mods_per_file = {}
    cleavages_per_file = {}
    for file in os.listdir(input_dir):
        if file.endswith('.xlsx'):
            mods_for_file, cleavages_for_file = process_protein_pilot_xlsx_file(input_dir+file)
            mods_for_file = set(mods_for_file)
            cleavages_for_file = set(cleavages_for_file)
            all_mods.extend(mods_for_file)
            all_cleavages.extend(cleavages_for_file)
            mods_per_file[file] = mods_for_file
            cleavages_per_file[file] = cleavages_for_file
    all_mods = sorted(set(all_mods), key=reader_helper.extract_index)
    all_cleavages = sorted(set(all_cleavages), key=reader_helper.extract_index)

    for file, mods in mods_per_file.items():
        if 'Acetyl' not in file:
            file_name = '_'.join(file.split('_')[0:3] + ['Acetyl'] + file.split('_')[3:])
            mods_per_file[file] = mods.union(mods_per_file[file_name])

    for file in list(mods_per_file.keys()):
        if 'Acetyl' in file:
            del mods_per_file[file]

    for file, cleavages in cleavages_per_file.items():
        if 'Acetyl' not in file:
            file_name = '_'.join(file.split('_')[0:3] + ['Acetyl'] + file.split('_')[3:])
            cleavages_per_file[file] = cleavages.union(cleavages_per_file[file_name])

    for file in list(cleavages_per_file.keys()):
        if 'Acetyl' in file:
            del cleavages_per_file[file]

    with open(f"{CONFIG.OUTPUT_FOLDER}/result_protein_pilot_mods.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + all_mods)
        writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mods])
        writer.writerow(['', ''] + [reader_helper.extract_mod_location(mod) for mod in all_mods])
        for file, mods in mods_per_file.items():
            row = [1 if mod in mods else 0 for mod in all_mods]
            group = groups_df.loc[groups_df['file_name'] == file]['group_name'].values[0]
            writer.writerow([file[:-10], group] + row)

    cleavages_with_ranges = []
    i = 0
    while i < len(all_cleavages):
        first_cleavage = all_cleavages[i]
        cleavage = all_cleavages[i]
        while i < len(all_cleavages)-1 and int(reader_helper.extract_cleavage_location(cleavage))+1 == int(reader_helper.extract_cleavage_location(all_cleavages[i+1])):
            i += 1
            cleavage = all_cleavages[i]
        if int(reader_helper.extract_cleavage_location(first_cleavage)) != int(reader_helper.extract_cleavage_location(cleavage)):
            cleavage_range = reader_helper.extract_cleavage_location(first_cleavage) + '-' + reader_helper.extract_cleavage_location(cleavage)
            cleavages_with_ranges.append(cleavage_range)
        else:
            cleavages_with_ranges.append(reader_helper.extract_cleavage_location(first_cleavage))
        i += 1
        
    with open(f"{CONFIG.OUTPUT_FOLDER}/result_protein_pilot_cleavages.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + cleavages_with_ranges)
        writer.writerow(['', ''] + ['Non-Tryptic' for _ in cleavages_with_ranges])
        writer.writerow(['', ''] + [cleavage for cleavage in cleavages_with_ranges])
        ranges = parse_ranges(cleavages_with_ranges)
        for file, cleavages in cleavages_per_file.items():
            indexes = [reader_helper.extract_index(cleavage) for cleavage in cleavages]
            row = cleavage_score(ranges, indexes)
            group = groups_df.loc[groups_df['file_name'] == file]['group_name'].values[0]
            writer.writerow([file[:-10], group] + row)

    return all_mods, all_cleavages


process_protein_pilot_dir()
