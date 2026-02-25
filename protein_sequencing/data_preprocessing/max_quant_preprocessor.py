"""MaxQuant preprocessor module. Extracts modifications and cleavages from MaxQuant output file."""

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from protein_sequencing import exon_helper, uniprot_align
from protein_sequencing.data_preprocessing import preprocessor_helper


class MaxQuantPreprocessor:
    """MaxQuant Preprocessor."""

    def __init__(
        self,
        config,
        preprocessor_config,
        evidence_df: pd.DataFrame,
        metadata_df: pd.DataFrame | None = None,
        metadata_column: str | None = None,
    ) -> None:
        self.CONFIG = config
        self.PREPROCESSOR_CONFIG = preprocessor_config

        self.fasta_file = self.PREPROCESSOR_CONFIG.FASTA_FILE
        self.aligned_fasta_file = self.PREPROCESSOR_CONFIG.ALIGNED_FASTA_FILE
        self.out_dir = config.OUTPUT_FOLDER
        if metadata_df is not None:
            if metadata_column not in metadata_df.columns:
                raise ValueError(
                    f"Metadata column '{metadata_column}' not found in metadata DataFrame columns: {metadata_df.columns}"
                )
            metadata_samples = set(metadata_df["Sample"].unique())
            evidence_samples = set(evidence_df["Sample"].unique())
            if not evidence_samples.issubset(metadata_samples):
                missing_groups = evidence_samples - metadata_samples
                raise ValueError(
                    f"The following samples from the evidence file are missing in the metadata file: {missing_groups}"
                )

            self.groups_df = metadata_df[["Sample", metadata_column]].rename(
                columns={"Sample": "file_name", metadata_column: "group_name"}
            )
        else:
            self.groups_df = pd.DataFrame(columns=["file_name", "group_name"])

        # Makes sure that groups are correctly extracted from the metadata file if provided
        assert (
            metadata_column is None
            and metadata_df is None
            or not self.groups_df.empty
        ), (
            "Metadata was provided but no valid groups could be extracted. Please check the metadata file and the "
            "specified metadata column."
        )

        aligned_sequence = uniprot_align.get_alignment(
            Path(self.fasta_file), Path(self.out_dir)
        )
        longest_original_sequence = max(
            [len(i) for i in aligned_sequence.alignment.inverse_indices]
        )
        max_configured_region = max(region[1] for region in config.REGIONS)
        if longest_original_sequence != max_configured_region:
            raise ValueError(
                f"The longest original sequence has a length of {longest_original_sequence}, but the regions file only "
                f"contains regions up to position {max_configured_region}. Please adjust the regions in the "
                "corresponding CSV file to match the sequence length."
            )
        self.sorted_isoform_headers = preprocessor_helper.process_fasta_files(
            self.fasta_file, self.aligned_fasta_file
        )

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
            self.max_sequence_leng,
        ) = exon_helper.retrieve_exon(
            Path(self.fasta_file), self.CONFIG.MIN_EXON_LENGTH, Path(self.out_dir)
        )

        remapped_evidence_df = self._remap_isoforms_in_evidence(evidence_df)
        self.process_max_quant_file(remapped_evidence_df)

    def get_exact_indices(self, mod_sequence: str) -> list:
        """Get exact indices of the modifications in the sequence."""
        indices = []
        current_index = 0
        inside_brackets = False
        inside_second_brackets = False
        for i, char in enumerate(mod_sequence):
            if char == "(":
                if inside_brackets:
                    inside_second_brackets = True
                else:
                    inside_brackets = True
            elif char == ")":
                if inside_second_brackets:
                    inside_second_brackets = False
                else:
                    inside_brackets = False
                    continue
            elif not inside_second_brackets and char.isalpha() or char == "_":
                if i + 1 < len(mod_sequence) and mod_sequence[i + 1] == "(":
                    indices.append(current_index)
            if not inside_brackets:
                current_index += 1

        return indices

    @staticmethod
    def _remap_isoforms_in_evidence(evidence_df: pd.DataFrame) -> pd.DataFrame:
        """Remap isoform IDs in the evidence DataFrame to match those in the fasta file."""
        def remap_protein_id(protein_id: str) -> str:
            protein_ids = protein_id.split(";")
            remapped_ids = []
            for protein_id in protein_ids:
                new_protein_id = f"{protein_id}-1" if "-" not in protein_id else protein_id
                remapped_ids.append(new_protein_id)

            return ";".join(remapped_ids)

        evidence_df["Protein ID"] = evidence_df["Protein ID"].apply(remap_protein_id)
        return evidence_df

    def reformat_mod(
        self,
        modified_peptide: str,
        peptide: str,
        peptide_offset: int,
        sequence: str,
        isoform: str,
        aligned_sequence: str,
    ) -> list[str]:
        """Reformat the modification string."""
        mod_strings = []

        pattern = r"\(([^()]*(\([^()]*\))*)\)"
        matches = re.findall(pattern, modified_peptide)
        indices = self.get_exact_indices(modified_peptide)
        for counter, (mod_site, _) in enumerate(matches):
            if " " in mod_site:
                mod_type = mod_site.split(" ")[0]
                mod_site = mod_site.split("(")[1][:-1]
            else:
                mod_type = mod_site.split(" ")[0]
                if mod_type == "ph":
                    mod_type = "Phospho"
                elif mod_type == "ac":
                    mod_type = "Acetyl"
                elif mod_type == "gg":
                    mod_type = "GG"
                elif mod_type == "me":
                    mod_type = "Methyl"
                elif mod_type == "ci":
                    mod_type = "Citrullination"
                elif mod_type == "de":
                    mod_type = "Deamidated"

            if mod_site == "Protein N-term" or mod_site == "Protein C-term" or indices[counter] == 0:
                # We skip these modifications because they are not as relevant
                continue

            assert indices[counter] > 0, "Modification site index should be greater than 0"
            mod_site_index = indices[counter] - 1
            aa = peptide[mod_site_index]
            aa_offset = indices[counter]

            if aa == "R" and mod_type == "Deamidated":
                mod_type = "Citrullination"

            missing_aa = 0
            if len(sequence) != len(aligned_sequence):
                missing_aa = preprocessor_helper.count_missing_amino_acids(
                    peptide[:aa_offset],
                    aligned_sequence,
                    peptide_offset,
                    self.exon_start_index,
                    self.exon_end_index,
                )
            offset = preprocessor_helper.calculate_exon_offset(
                aa_offset + peptide_offset + missing_aa,
                isoform,
                self.exon_found,
                self.exon_end_index,
                self.exon_1_isoforms,
                self.exon_2_isoforms,
                self.exon_1_length,
                self.exon_2_length,
                self.exon_length,
            )
            assert offset < len(aligned_sequence)
            if aligned_sequence[offset - 1] != aa:
                raise ValueError(
                    f"Amino acid doesn't match for {aa} for peptide {peptide} in sequence {sequence} with offset {offset}"
                )
            iso = preprocessor_helper.get_isoform_for_offset(
                isoform,
                offset,
                self.exon_start_index,
                self.exon_1_isoforms,
                self.exon_1_length,
                self.exon_2_isoforms,
                self.exon_2_length,
            )
            mod_strings.append(f"{mod_type}({aa})@{offset}_{iso}")
        return mod_strings

    def process_max_quant_file(self, evidence_df: pd.DataFrame):
        """Process MaxQuant file."""
        all_mods = []
        mods_for_exp = defaultdict(list)
        all_cleavages = []
        cleavages_for_exp = defaultdict(list)
        group_names = self.groups_df["file_name"].values

        isoforms_from_fasta = [i for (i, _, _) in self.sorted_isoform_headers]
        filtered_evidence_df = evidence_df[
            evidence_df["Protein ID"].apply(
                lambda x: any(
                    isoform in x.split(";") for isoform in isoforms_from_fasta
                )
            )
        ]

        if len(filtered_evidence_df) == 0:
            raise ValueError(
                "No matching isoform IDs found between the uploaded evidence file and the fasta file."
            )

        for _, row in filtered_evidence_df.iterrows():
            try:
                isoform, sequence, peptide_offset, aligned_sequence = (
                    preprocessor_helper.get_accession(
                        row["Protein ID"], row["Sequence"], self.sorted_isoform_headers
                    )
                )  # Should be enough for now
            except ValueError:
                continue

            cleavage = preprocessor_helper.check_N_term_cleavage(
                row["Sequence"],
                row["Protein ID"],
                self.sorted_isoform_headers,
                self.exon_found,
                self.exon_start_index,
                self.exon_end_index,
                self.exon_1_isoforms,
                self.exon_2_isoforms,
                self.exon_1_length,
                self.exon_2_length,
                self.exon_length,
            )
            if cleavage != "":
                all_cleavages.append(cleavage)
                if len(self.groups_df) > 0:
                    cleavages_for_exp[row["Sample"]].append(cleavage)

            cleavage = preprocessor_helper.check_C_term_cleavage(
                row["Sequence"],
                row["Protein ID"],
                self.sorted_isoform_headers,
                self.exon_found,
                self.exon_start_index,
                self.exon_end_index,
                self.exon_1_isoforms,
                self.exon_2_isoforms,
                self.exon_1_length,
                self.exon_2_length,
                self.exon_length,
            )
            if cleavage != "":
                all_cleavages.append(cleavage)
                if len(self.groups_df) > 0:
                    cleavages_for_exp[row["Sample"]].append(cleavage)

            if float(row["PEP"]) < self.PREPROCESSOR_CONFIG.THRESHOLD:
                if row["Modifications"] != "Unmodified":
                    mods = self.reformat_mod(
                        row["Modified sequence"],
                        row["Sequence"],
                        peptide_offset,
                        sequence,
                        isoform,
                        aligned_sequence,
                    )
                    all_mods.extend(mods)
                    if row["Sample"] in group_names:
                        mods_for_exp[row["Sample"]].extend(mods)

        all_mods = sorted(set(all_mods), key=preprocessor_helper.extract_index)
        all_mods = preprocessor_helper.sort_by_index_and_exons(all_mods)
        for key in mods_for_exp:
            mods_for_exp[key] = sorted(
                set(mods_for_exp[key]), key=preprocessor_helper.extract_index
            )
            mods_for_exp[key] = preprocessor_helper.sort_by_index_and_exons(
                mods_for_exp[key]
            )

        all_cleavages = sorted(
            set(all_cleavages), key=preprocessor_helper.extract_cleavage_location
        )
        all_cleavages = preprocessor_helper.sort_by_index_and_exons(all_cleavages)
        cleavages_with_ranges = preprocessor_helper.extract_cleavages_ranges(
            all_cleavages
        )
        preprocessor_helper.write_results(
            all_mods,
            mods_for_exp,
            cleavages_with_ranges,
            cleavages_for_exp,
            f"{self.CONFIG.OUTPUT_FOLDER}/result_max_quant",
            self.groups_df,
        )
