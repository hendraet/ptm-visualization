"""Mascot Preprocessor Module. Extracting modifications from Mascot results and creating a CSV file with the results."""
import os
import csv
import re
from pathlib import Path

import pandas as pd
from protein_sequencing import exon_helper, uniprot_align
from protein_sequencing.data_preprocessing import preprocessor_helper

class MascotPreprocessor:
    """Mascot Preprocessor Class."""

    def __init__(self, config, preprocessor_config) -> None:
        self.CONFIG = config
        self.PREPROCESSOR_CONFIG = preprocessor_config

        self.fasta_file = self.PREPROCESSOR_CONFIG.FASTA_FILE
        self.aligned_fasta_file = self.PREPROCESSOR_CONFIG.ALIGNED_FASTA_FILE
        self.input_dir = self.PREPROCESSOR_CONFIG.MASCOT_INPUT_DIR

        uniprot_align.get_alignment(self.fasta_file)
        self.fasta_headers = preprocessor_helper.process_tau_file(self.fasta_file, self.aligned_fasta_file)

        self.groups_df = pd.read_csv(self.PREPROCESSOR_CONFIG.GROUPS_CSV)
        (
            self.exon_found,
            self.exon_start_index,
            self.exon_end_index,
            self.exon_length,
            self.exon_1_isoforms,
            self.exon_1_length,
            self.exon_2_isoforms,
            self.exon_2_length,
            self.exon_none_isoforms,
            self.max_sequence_length
        ) = exon_helper.retrieve_exon(self.fasta_file, self.CONFIG.MIN_EXON_LENGTH)

        self.process_mascot_dir()

    def reformmods(self, sites, peptide, variable_mods, isoform, sequence, aligned_sequence):
        """Reformat the modification string to have the correct indexes and exon information."""
        modstrings = []
        sites = sites.split('.')[1]
        for i, site in enumerate(sites):
            site = int(site)
            if site != 0:
                mod = variable_mods[site]
                aa = peptide[i]
                if self.CONFIG.INCLUDED_MODIFICATIONS.get(mod):
                    if aa not in self.CONFIG.INCLUDED_MODIFICATIONS[mod]:
                        continue
                    if aa == 'R' and mod == 'Deamidated':
                        mod = 'Citrullination'
                peptide_offset = sequence.index(peptide)+1
                missing_aa = 0
                if len(sequence) != len(aligned_sequence):
                    missing_aa = preprocessor_helper.count_missing_amino_acids(peptide[:i+1], aligned_sequence, peptide_offset, self.exon_start_index, self.exon_end_index)
                offset = preprocessor_helper.calculate_exon_offset(peptide_offset + i + missing_aa, isoform, self.exon_found, self.exon_end_index, self.exon_1_isoforms, self.exon_2_isoforms, self.exon_1_length, self.exon_2_length, self.exon_length)
                if aligned_sequence[offset-1] != aa:
                    raise ValueError(f"AA don't match for {aa} for peptide {peptide} in sequence {sequence} with offset {offset}")
                iso = preprocessor_helper.get_isoform_for_offset(isoform, offset, self.exon_start_index, self.exon_1_isoforms, self.exon_1_length, self.exon_2_isoforms, self.exon_2_length)
                modstring = f"{mod}({aa})@{offset}_{iso}"
                modstrings.append(modstring)
        return modstrings

    def process_mascot_file(self, file, fasta_dict):
        """Process a Mascot file and extract the modifications."""
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
        with open(self.input_dir + file, 'r', encoding="utf-8") as f:
            while line := f.readline():
                if line.startswith('\"Variable modifications'):
                    varmods = True
                elif line.startswith('\"Protein hits'):
                    table = True
                elif varmods and line == '\n':
                    emptyrowpars = True
                elif table and line == '\n':
                    emptyrowpars = True
                elif '--------------------------------------------------------' in line:
                    varmods = False
                    table = False
                    emptyrowpars = False
                elif varmods and emptyrowpars:
                    if not line.startswith('\"Identifier') and line != '\n':
                        mod_string = line.split(',')[1]
                        mod = mod_string.split(' ')[0][1:]
                        variable_mods[int(line[0])] = mod
                elif table and emptyrowpars:
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
                        header_found = False
                        search_header = row[pep_accession_idx].strip('"\'')
                        if search_header in self.PREPROCESSOR_CONFIG.ISOFORM_HELPER_DICT:
                            search_header = self.PREPROCESSOR_CONFIG.ISOFORM_HELPER_DICT[search_header]
                        for fasta_header in fasta_dict:
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
                                try:
                                    assert len(line.split(',')) == header_length
                                except AssertionError:
                                    print(f"Failed automatically fixing line: {line}")
                                    continue
                                row = line.split(',')
                            else:
                                print(f"Error with line: {line}")
                                continue
                        if row[pep_var_mod_idx] == '""':
                            continue
                        isoform = None
                        sequence = None
                        aligned_sequence = None
                        for fasta_header in fasta_dict:
                            if row[pep_seq_idx] in fasta_header[1]:
                                isoform = fasta_header[0]
                                sequence = fasta_header[1]
                                aligned_sequence = fasta_header[2]
                                break
                        if isoform is None:
                            continue
                        reform_mod_strings = self.reformmods(row[pep_var_mod_pos_idx], row[pep_seq_idx], variable_mods, isoform, sequence, aligned_sequence)
                        all_mod_strings.extend(reform_mod_strings)
        return all_mod_strings

    def process_results(self, all_mod_strings, mod_strings_for_files):
        """Process the results and write it to a CSV file."""
        all_mod_strings = sorted(set(all_mod_strings), key=preprocessor_helper.extract_index)
        all_mods = preprocessor_helper.sort_by_index_and_exons(all_mod_strings)
        out_dir = Path(self.CONFIG.OUTPUT_FOLDER)
        if not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)
        with (out_dir / "result_mascot.csv").open('w', newline='', encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(['ID', 'Group'] + all_mods)
            writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mods])
            writer.writerow(['', ''] + [preprocessor_helper.extract_mod_location(mod) for mod in all_mods])
            writer.writerow(['', ''] + [mod.split('_')[1] for mod in all_mods])
            for file, mods in mod_strings_for_files.items():
                row = [1 if mod in mods else 0 for mod in all_mods]
                if file not in self.groups_df['file_name'].values:
                    raise KeyError(f"File {file} not found in groups CSV")
                group = self.groups_df.loc[self.groups_df['file_name'] == file]['group_name'].values[0]
                writer.writerow([file, group] + row)


    def process_mascot_dir(self):
        """Process all Mascot files in a directory."""
        files = os.listdir(self.input_dir)
        all_mod_strings = []
        mod_strings_for_files = {}
        for file in files:
            result= self.process_mascot_file(file, self.fasta_headers)
            all_mod_strings.extend(result)
            mod_strings_for_files[file] = result

        self.process_results(all_mod_strings, mod_strings_for_files)
