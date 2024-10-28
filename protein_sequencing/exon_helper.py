import os
from protein_sequencing import uniprot_align, utils

def levenshtein_distance(str1, str2, MIN_EXON_LENGTH):
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
                    matrix[i - 1][j],       # Deletion
                    matrix[i][j - 1],       # Insertion
                    matrix[i - 1][j - 1]    # Substitution
                )

    return matrix[-1][-1] <= MIN_EXON_LENGTH

def retrieve_exon(input_file: str | os.PathLike, MIN_EXON_LENGTH: int):
    alignments = list(uniprot_align.get_alignment(input_file))
    max_sequence_length = 0
    for alignment in alignments:
        if len(alignment.seq) > max_sequence_length:
            max_sequence_length = len(alignment.seq)

    assert all(len(alignment.seq) == max_sequence_length for alignment in alignments)
    for alignment in alignments:
        utils.ISOFORM_IDS.append(alignment.id.split('|')[1])

    different_possibilities = [-1]*max_sequence_length
    for i in range(len(alignments[0].seq)):
        amino_acids = dict(list())
        for alignment in alignments:
            amino_acid = alignment.seq[i]
            if amino_acid not in amino_acids:
                amino_acids[amino_acid] = [alignment.id.split('|')[1]]
            else:
                amino_acids[amino_acid].append(alignment.id.split('|')[1])
    
        if '-' in amino_acids:
            if len(amino_acids) == 2:
                different_possibilities[i] = -1
            if len(amino_acids) > 2:
                different_possibilities[i] = len(amino_acids)-1
        else:
            if len(amino_acids) > 1:
                different_possibilities[i] = len(amino_acids)
            else:
                different_possibilities[i] = 1

    i = 0
    exon_start_index = -1
    max_exon_length = 0
    exon_found = False
    while i < len(different_possibilities):
        if different_possibilities[i] == 2:
            if exon_found:
                raise ValueError("There are multiple exons in the sequence, currently the tool just supports 1 different exon.")
            current_position = i
            exon_start_index = i
            while i > 0:
                if different_possibilities[i] == -1:
                    exon_start_index = i
                    i -= 1
                else:
                    break
            i = current_position
            while i < len(different_possibilities):
                if different_possibilities[i] == -1 or different_possibilities[i] == 2:
                    i += 1
                else:
                    break
            exon_end_index = i
            min_alignment_offset = -1
            for alignment in alignments:
                alignemnt_offset = alignment.seq[exon_start_index:exon_end_index].count('-')
                if min_alignment_offset == -1 or alignemnt_offset < min_alignment_offset:
                    min_alignment_offset = alignemnt_offset
                    
            max_exon_length = exon_end_index - exon_start_index - min_alignment_offset
            # filters minor differences in aa sequence which should not count as a different exon
            if max_exon_length < MIN_EXON_LENGTH:
                max_exon_length = 0
                exon_start_index = -1
            else:
                exon_found = True
                exon_1 = None
                exon_2 = None
                exon_1_isoforms = []
                exon_2_isoforms = []
                exon_none_isoforms = []
                exon_1_length = -1
                exon_2_length = -1
                for alignment in alignments:
                    exon = alignment.seq[exon_start_index:exon_end_index].replace('-', '')
                    isoform = alignment.id.split('|')[1]
                    if exon != '' and len(exon) > MIN_EXON_LENGTH:
                        if exon_1 is None:
                            exon_1 = exon
                            exon_1_isoforms.append(isoform)
                            exon_1_length = len(exon)
                        elif levenshtein_distance(exon_1, exon, MIN_EXON_LENGTH):
                            exon_1_isoforms.append(isoform)
                        elif exon_2 is None:
                            exon_2 = exon
                            exon_2_isoforms.append(isoform)
                            exon_2_length = len(exon)
                        elif levenshtein_distance(exon_1, exon, MIN_EXON_LENGTH):
                            exon_2_isoforms.append(isoform)
                        else:
                            raise ValueError("There are more than 2 different exons in the sequence.")
                    else:
                        exon_none_isoforms.append(isoform)
        i += 1

    if exon_found:
        # exon_start_index starts with 0 (so the first amino acid in the exon is +1)
        return True, exon_start_index, exon_end_index, max_exon_length, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length, exon_none_isoforms, max_sequence_length
    
    else:
        return False, -1, -1, -1, [], -1, [], -1, [alignment.id.split('|')[1] for alignment in alignments], max_sequence_length
