"""Helper functions for exon retrieval"""

from pathlib import Path

import pandas as pd

from protein_sequencing import uniprot_align


def levenshtein_distance(str1: str, str2: str, min_exon_length: int) -> bool:
    """Calculate the Levenshtein distance between two strings.
    Args:
        str1 (str): First exon.
        str2 (str): Second exon.
        min_exon_length (int): Minimum exon length.
    Returns:
        bool: True if the Levenshtein distance is less than or equal to the minimum exon length.
    """
    if abs(len(str1) - len(str2)) > 1:
        return False
    matrix = [[0] * (len(str2) + 1) for _ in range(len(str1) + 1)]
    for i in range(len(str1) + 1):
        matrix[i][0] = i
    for j in range(len(str2) + 1):
        matrix[0][j] = j
    for i in range(1, len(str1) + 1):
        for j in range(1, len(str2) + 1):
            if str1[i - 1] == str2[j - 1]:
                matrix[i][j] = matrix[i - 1][j - 1]
            else:
                matrix[i][j] = 1 + min(
                    matrix[i - 1][j],  # Deletion
                    matrix[i][j - 1],  # Insertion
                    matrix[i - 1][j - 1],  # Substitution
                )

    return matrix[-1][-1] <= min_exon_length


def add_exon_info_to_regions(
    fasta_file: Path, regions_df: pd.DataFrame, min_exon_length: int, out_dir: Path
) -> tuple[pd.DataFrame, int | None, int | None, int | None, int | None, int | None, int]:
    """Retrieve exon from protein sequence."""
    alignments = list(uniprot_align.get_alignment(fasta_file, out_dir))
    max_sequence_length = 0
    for alignment in alignments:
        max_sequence_length = max(max_sequence_length, len(alignment.seq))

    assert all(len(alignment.seq) == max_sequence_length for alignment in alignments)

    # TODO: would be nice to replace the -1 with None
    different_possibilities = [-1] * max_sequence_length
    for i in range(max_sequence_length):
        amino_acids = dict(list())
        for alignment in alignments:
            amino_acid = alignment.seq[i]
            if amino_acid not in amino_acids:
                amino_acids[amino_acid] = [alignment.id.split("|")[1]]
            else:
                amino_acids[amino_acid].append(alignment.id.split("|")[1])

        if "-" in amino_acids:
            if len(amino_acids) == 2:
                different_possibilities[i] = -1
            if len(amino_acids) > 2:
                different_possibilities[i] = len(amino_acids) - 1
        else:
            if len(amino_acids) > 1:
                different_possibilities[i] = len(amino_acids)
            else:
                different_possibilities[i] = 1

    i = 0
    exon_start_index = None
    max_exon_length = 0
    exon_found = False
    while i < len(different_possibilities):
        # 2 different possibilities means that there are 2 different amino acids at this position, which is the minimum
        # for an exon.
        if different_possibilities[i] == 2:
            if exon_found:
                raise ValueError(
                    "There are multiple exons in the sequence, currently the tool just supports 1 different exon."
                )

            current_position = i
            exon_start_index = i

            # Go left and right until we find a position where all sequences are the same, which should be the end of
            #  the exon.
            while i > 0:
                i -= 1
                if different_possibilities[i] == -1:
                    exon_start_index = i
                else:
                    break
            i = current_position

            # Go right until we find a position where all sequences are the same again, which should be the other
            #  end of the exon.
            while i < len(different_possibilities):
                if different_possibilities[i] == -1 or different_possibilities[i] == 2:
                    i += 1
                else:
                    break
            exon_end_index = i

            # We have the start and end of the exon, but we need to check if there are minor differences in the
            # sequence which should not count as a different exon. We do this by checking the number of gaps in the
            # alignment for each sequence in the exon region and taking the minimum number of gaps as the number of
            # minor differences. We then subtract this from the exon length to get the maximum exon length.
            min_alignment_offset = None
            for alignment in alignments:
                alignment_offset = alignment.seq[exon_start_index:exon_end_index].count(
                    "-"
                )
                if (
                    min_alignment_offset is None
                    or alignment_offset < min_alignment_offset
                ):
                    min_alignment_offset = alignment_offset

            max_exon_length = exon_end_index - exon_start_index - min_alignment_offset
            # filters minor differences in aa sequence which should not count as a different exon
            if max_exon_length < min_exon_length:
                max_exon_length = 0
                exon_start_index = -1
            else:
                exon_found = True
                exon_1 = None
                exon_2 = None
                exon_1_isoforms = []
                exon_2_isoforms = []
                exon_none_isoforms = []
                exon_1_length = 0
                exon_2_length = 0
                # TODO: could also be calculated outside
                max_sequence_length = (
                    max_sequence_length
                    - (exon_end_index - exon_start_index)
                    + max_exon_length
                )
                for alignment in alignments:
                    exon = alignment.seq[exon_start_index:exon_end_index].replace(
                        "-", ""
                    )
                    isoform = alignment.id.split("|")[1]
                    if exon != "" and len(exon) > min_exon_length:
                        if exon_1 is None:
                            exon_1 = exon
                            exon_1_isoforms.append(isoform)
                            exon_1_length = len(exon)
                        elif levenshtein_distance(exon_1, exon, min_exon_length):
                            exon_1_isoforms.append(isoform)
                        elif exon_2 is None:
                            exon_2 = exon
                            exon_2_isoforms.append(isoform)
                            exon_2_length = len(exon)
                        elif levenshtein_distance(exon_1, exon, min_exon_length):
                            exon_2_isoforms.append(isoform)
                        else:
                            raise ValueError(
                                "There are more than 2 different exons in the sequence."
                            )
                    else:
                        exon_none_isoforms.append(isoform)
        i += 1

    new_regions_df = regions_df.copy()
    new_regions_df["exon_id"] = None
    # TODO: maybe change this to sth. else? Empty lists?
    new_regions_df["isoforms"] = None

    # the additional length checks ensure that we don't count small differences in the sequence as an exon, which can
    # happen for example with cassette exons where one isoform has a single amino acid more than the other isoform.
    if exon_found and exon_1_length > 0 and exon_2_length > 0:
        # TODO: is there a problem with how (casette) exons with single substituion are handled, e.g. one exon length
        #  is 0 while there is one substitution
        # TODO: maybe verify this by mocking an exon at the beginning of the sequence
        # exon_start_index points to the first amino acid of the exon and exon_end_index points to the last amino acid
        # of the exon

        exon_1_end_index = exon_start_index + exon_1_length
        exon_1 = new_regions_df[
            # TODO: we have to find a way to deal with these fucking +1
            # TODO: is exon_1_end_index correct here? be also +1? - double check against data later
            (new_regions_df["region_start"] == exon_start_index + 1)
            & (new_regions_df["region_end"] == exon_1_end_index)
        ]
        if len(exon_1) != 1:
            # TODO: write/adapt test
            raise ValueError(
                "Could not find a unique exon 1 in regions file. Expected to find exactly one exon with start "
                f"{exon_start_index + 1} and end {exon_1_end_index} but found {len(exon_1)} exon(s)"
            )
        exon_2_end_index = exon_start_index + exon_2_length
        exon_2 = new_regions_df[
            (new_regions_df["region_start"] == exon_start_index + 1)
            & (new_regions_df["region_end"] == exon_2_end_index)
        ]
        if len(exon_2) != 1:
            # TODO: write/adapt test
            raise ValueError(
                "Could not find a unique exon 2 in regions file. Expected to find exactly one exon with start "
                f"{exon_start_index + 1} and end {exon_2_end_index} but found {len(exon_1)} exon(s)"
            )

        new_regions_df.loc[exon_1.index, "exon_id"] = 1
        new_regions_df.loc[exon_2.index, "exon_id"] = 2

        # TODO: maybe we have to add none_isoforms here
        new_regions_df.loc[exon_1.index, "isoforms"] = exon_1_isoforms
        new_regions_df.loc[exon_2.index, "isoforms"] = exon_2_isoforms

        # TODO: check what isoforms are actually used for (outside exons) - might not be relevant at all
        # new_regions_df["isoforms"] = [
        #     alignment.id.split("|")[1] for alignment in alignments
        # ]

        # Sanity checks
        for i, region in new_regions_df.iterrows():
            current_region_end = region["region_end"]
            if current_region_end == exon_start_index:
                # We first check that there are actually two exons after the current region to account for
                # potential exons at the end of the sequence
                if exon_1_length > 0 and exon_2_length > 0 and len(new_regions_df) <= i + 2:
                    raise ValueError(
                        f"Exon start {exon_start_index} matches a region end for region {region}, but "
                        "there are not enough regions after it, please check your supplied region "
                        "list."
                    )

                raise ValueError(
                    f"Exon start {exon_start_index} does not match any region end, please check your supplied region "
                    "list - maybe it is missing some regions or it doesn't match the provided fasta sequence."
                )
    else:
        exon_start_index = None
        exon_end_index = None
        exon_1_length = None
        exon_2_length = None
        max_exon_length = None


    return (
        new_regions_df,
        exon_start_index,
        exon_end_index,
        exon_1_length,
        exon_2_length,
        max_exon_length,
        max_sequence_length,
    )

    # TODO: remove
    # return (
    #     True, # can probably also be checked by looking at regions
    #     exon_start_index,  #
    #     exon_end_index,
    #     max_exon_length,  # calc
    #     exon_1_isoforms,  #
    #     exon_1_length,   # calc
    #     exon_2_isoforms,  #
    #     exon_2_length,  # calc
    #     exon_none_isoforms,  # return?
    #     max_sequence_length,  # return
    # )

    # TODO: all these -1 here and also above is just desaster waiting to happen. Ideally, these should be None and
    #  there would be code in place that checks for it.
    # TODO: remove
    # return (
    #     False,
    #     -1,
    #     -1,
    #     -1,
    #     [],
    #     -1,
    #     [],
    #     -1,
    #     [alignment.id.split("|")[1] for alignment in alignments],
    #     max_sequence_length,
    # )
