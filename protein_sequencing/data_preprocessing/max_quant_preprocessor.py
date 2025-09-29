"""MaxQuant preprocessor module. Extracts modifications and cleavages from MaxQuant output file."""
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from protein_sequencing import exon_helper, uniprot_align
from protein_sequencing.data_preprocessing import preprocessor_helper


class MaxQuantPreprocessor:
    """MaxQuant Preprocessor."""

    def __init__(self, config, preprocessor_config, evidence_df: pd.DataFrame) -> None:
        self.CONFIG = config
        self.PREPROCESSOR_CONFIG = preprocessor_config

        self.fasta_file = self.PREPROCESSOR_CONFIG.FASTA_FILE
        self.aligned_fasta_file = self.PREPROCESSOR_CONFIG.ALIGNED_FASTA_FILE
        self.out_dir = config.OUTPUT_FOLDER

        uniprot_align.get_alignment(Path(self.fasta_file), Path(self.out_dir))
        self.sorted_isoform_headers = preprocessor_helper.process_fasta_files(self.fasta_file, self.aligned_fasta_file)

        group_file_keys = {'file_name', 'group_name', 'replicate'}
        self.groups_df = pd.read_csv(gf) if (gf := self.PREPROCESSOR_CONFIG.GROUPS_CSV) is not None else pd.DataFrame(
            {k: [] for k in group_file_keys}
        )
        assert group_file_keys.issubset(set(self.groups_df.columns)), (f"Groups file must contain the columns: "
                                                                       f"{group_file_keys}")
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
            self.max_sequence_leng
        ) = exon_helper.retrieve_exon(Path(self.fasta_file), self.CONFIG.MIN_EXON_LENGTH, Path(self.out_dir))

        self.process_max_quant_file(evidence_df)

    def get_exact_indexes(self, mod_sequence: str) -> list:
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

    def reformat_mod(self, modified_peptide: str, peptide: str, peptide_offset: int, sequence: str, isoform: str, aligned_sequence: str) -> list[str]:
        """Reformat the modification string."""
        mod_strings = []

        pattern = r"\(([^()]*(\([^()]*\))*)\)"
        matches = re.findall(pattern, modified_peptide)
        indexes = self.get_exact_indexes(modified_peptide)
        counter = 0
        for mod_position, _ in matches:
            if ' ' in mod_position:
                mod_type = mod_position.split(' ')[0]
                mod_position = mod_position.split('(')[1][:-1]
            else:
                mod_type = mod_position.split(' ')[0]
                if mod_type == 'ph':
                    mod_type = 'Phospho'
                elif mod_type == 'ac':
                    mod_type = 'Acetyl'
                elif mod_type == 'gg':
                    mod_type = 'GG'
                elif mod_type == 'me':
                    mod_type = 'Methyl'
                elif mod_type == 'ci':
                    mod_type = 'Citrullination'
                elif mod_type == 'de':
                    mod_type = 'Deamidated'
            if mod_position == "Protein N-term":
                aa = sequence[peptide_offset-1]
                aa_offset = 0
            elif mod_position == "Protein C-term":
                aa = peptide[-1]
                aa_offset = len(peptide)
            else:
                aa = peptide[indexes[counter]-1]
                aa_offset = indexes[counter]
            if self.CONFIG.INCLUDED_MODIFICATIONS.get(mod_type):
                if aa not in self.CONFIG.INCLUDED_MODIFICATIONS[mod_type]:
                    continue
                if aa == 'R' and mod_type == 'Deamidated':
                    mod_type = 'Citrullination'
            else:
                continue
            missing_aa = 0
            if len(sequence) != len(aligned_sequence):
                missing_aa = preprocessor_helper.count_missing_amino_acids(peptide[:aa_offset], aligned_sequence, peptide_offset, self.exon_start_index, self.exon_end_index)
            offset = preprocessor_helper.calculate_exon_offset(aa_offset+peptide_offset+missing_aa, isoform, self.exon_found, self.exon_end_index, self.exon_1_isoforms, self.exon_2_isoforms, self.exon_1_length, self.exon_2_length, self.exon_length)
            if aligned_sequence[offset-1] != aa:
                raise ValueError(f"AA don't match for {aa} for peptide {peptide} in sequence {sequence} with offset {offset}")
            iso = preprocessor_helper.get_isoform_for_offset(isoform, offset, self.exon_start_index, self.exon_1_isoforms, self.exon_1_length, self.exon_2_isoforms, self.exon_2_length)
            mod_strings.append(f"{mod_type}({aa})@{offset}_{iso}")
            counter += 1
        return mod_strings

    def process_max_quant_file(self, evidence_df: pd.DataFrame):
        """Process MaxQuant file."""
        all_mods = []
        mods_for_exp = defaultdict(list)
        all_cleavages = []
        cleavages_for_exp = defaultdict(list)
        group_names = self.groups_df['file_name'].values

        isoforms_from_fasta = [i for (i, _, _) in self.sorted_isoform_headers]
        filtered_evidence_df = evidence_df[evidence_df['Protein ID'].apply(lambda x: any(
            isoform in x.split(';') for isoform in isoforms_from_fasta
        ))]

        if len(filtered_evidence_df) == 0:
            raise ValueError("No matching isoform IDs found between the uploaded evidence file and the fasta file.")

        for _, row in filtered_evidence_df.iterrows():
            try:
                isoform, sequence, peptide_offset, aligned_sequence = preprocessor_helper.get_accession(
                    row['Protein ID'],
                    row['Sequence'],
                    self.sorted_isoform_headers
                )  # Should be enough for now
            except ValueError:
                continue

            cleavage = preprocessor_helper.check_N_term_cleavage(
                row['Sequence'],
                row['Protein ID'],
                self.sorted_isoform_headers,
                self.exon_found,
                self.exon_start_index,
                self.exon_end_index,
                self.exon_1_isoforms,
                self.exon_2_isoforms,
                self.exon_1_length,
                self.exon_2_length,
                self.exon_length
            )
            if cleavage != "":
                all_cleavages.append(cleavage)
                if len(self.groups_df) > 0:
                    cleavages_for_exp[row['Sample']].append(cleavage)

            cleavage = preprocessor_helper.check_C_term_cleavage(row['Sequence'], row['Protein ID'],
                                                                 self.sorted_isoform_headers, self.exon_found,
                                                                 self.exon_start_index, self.exon_end_index,
                                                                 self.exon_1_isoforms, self.exon_2_isoforms,
                                                                 self.exon_1_length, self.exon_2_length,
                                                                 self.exon_length)
            if cleavage != "":
                all_cleavages.append(cleavage)
                if len(self.groups_df) > 0:
                    cleavages_for_exp[row['Sample']].append(cleavage)

            if float(row['PEP']) < self.PREPROCESSOR_CONFIG.THRESHOLD:
                if row['Modifications'] != "Unmodified":
                    mods = self.reformat_mod(
                        row['Modified sequence'],
                        row['Sequence'],
                        peptide_offset,
                        sequence,
                        isoform,
                        aligned_sequence
                    )
                    all_mods.extend(mods)
                    if row['Sample'] in group_names:
                        mods_for_exp[row['Sample']].extend(mods)

        all_mods = sorted(set(all_mods), key=preprocessor_helper.extract_index)
        all_mods = preprocessor_helper.sort_by_index_and_exons(all_mods)
        for key in mods_for_exp:
            mods_for_exp[key] = sorted(set(mods_for_exp[key]), key=preprocessor_helper.extract_index)
            mods_for_exp[key] = preprocessor_helper.sort_by_index_and_exons(mods_for_exp[key])

        all_cleavages = sorted(set(all_cleavages), key=preprocessor_helper.extract_cleavage_location)
        all_cleavages = preprocessor_helper.sort_by_index_and_exons(all_cleavages)
        cleavages_with_ranges = preprocessor_helper.extract_cleavages_ranges(all_cleavages)
        preprocessor_helper.write_results(
            all_mods,
            mods_for_exp,
            cleavages_with_ranges,
            cleavages_for_exp,
            f"{self.CONFIG.OUTPUT_FOLDER}/result_max_quant",
            self.groups_df
        )
