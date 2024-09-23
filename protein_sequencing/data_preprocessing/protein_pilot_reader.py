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

current_dir = os.path.dirname(__file__)
groups_df = pd.read_csv(f"{current_dir}/groups.csv")

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
    for index, row in enumerate(rows):
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
    if row[accession_index].split('|')[1] not in [header[0] for header in sorted_isoform_headers]:
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

def extract_mods_from_rows(rows, protein_mod_index, mod_index, seq_index, accession_index):
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
                modstring = f"{mod_name}({amino_acid})@{mod_pos+index_offset}_{isoform}"
                mods.append(modstring)

    return mods

def extract_cleavages_from_rows(rows, cleavage_index, seq_index, accession_index):
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
            
            cleavages.append(f"{amino_acid}@{index_offset+site_index}_{isoform}")

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

        


def process_protein_pilot_xlsx_file(file):
    workbook = CalamineWorkbook.from_path(file)
    if 'Distinct Peptide Level Data' in workbook.sheet_names:
        distinct_peptide_level_data = workbook.get_sheet_by_name('Distinct Peptide Level Data')
        fdr_threshold = extract_confidence_score(distinct_peptide_level_data)
    else:
        print("No 'Distinct Peptide Level Data' worksheet found in the Excel file.")
    if 'Peptide Summary' in workbook.sheet_names:
        peptide_summary = workbook.get_sheet_by_name('Peptide Summary')
        mods, cleavages = extract_data_with_threshold(peptide_summary, fdr_threshold)

        return mods, cleavages

def process_protein_pilot_dir():
    all_mods = []
    all_cleavages = []
    mods_per_file = {}
    cleavages_per_file = {}
    for file in os.listdir(input_dir):
        if file.endswith('.xlsx'):
            mods_for_file, cleavages_for_file = process_protein_pilot_xlsx_file(input_dir+file)
            all_mods.extend(mods_for_file)
            all_cleavages.extend(cleavages_for_file)
            mods_per_file[file] = mods_for_file
            cleavages_per_file[file] = cleavages_for_file

    all_mods = sorted(set(all_mods), key=reader_helper.extract_index)
    all_cleavages = sorted(set(all_cleavages), key=reader_helper.extract_index)

    with open(f"{CONFIG.OUTPUT_FOLDER}/result_protein_pilot_mods.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Neuropathology'] + all_mods)
        writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mods])
    return all_mods, all_cleavages


process_protein_pilot_dir()
