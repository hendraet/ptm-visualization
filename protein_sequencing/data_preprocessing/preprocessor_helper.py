"""Helper functions for the preprocessor module"""
from typing import Tuple
import csv

def process_tau_file(fasta_file, aligned_fasta_file):
    """Extracts the sequences from the fasta file and the aligned fasta file
    and returns them as a list of tuples sorted by the sequence length."""
    headers = []
    aligned_sequences = {}
    with open(aligned_fasta_file, 'r', encoding="utf-8") as file:
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

    with open(fasta_file, 'r', encoding="utf-8") as file:
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


def extract_index(string):
    """Extracts the index from a string."""
    return int(string.split('@')[1].split('_')[0])

def extract_mod_location(mod_string):
    """Extracts the location of the modification from the modification string."""
    return mod_string.split('(')[1].split(')')[0]+mod_string.split('@')[1].split('_')[0]

def extract_cleavage_location(cleavage_string):
    """Extracts the location of the cleavage from the cleavage string."""
    return cleavage_string.split('@')[1].split('_')[0]

def extract_cleavages_ranges(all_cleavages):
    """Extracts the cleavages ranges from the cleavages list."""
    cleavages_with_ranges = []
    i = 0
    while i < len(all_cleavages):
        first_cleavage = all_cleavages[i]
        isoform = first_cleavage.split('_')[1]
        cleavage = all_cleavages[i]
        while i < len(all_cleavages)-1 and int(extract_cleavage_location(cleavage))+1 == int(extract_cleavage_location(all_cleavages[i+1])) and isoform == all_cleavages[i+1].split('_')[1]:
            i += 1
            cleavage = all_cleavages[i]
        if int(extract_cleavage_location(first_cleavage)) != int(extract_cleavage_location(cleavage)):
            cleavage_range = extract_cleavage_location(first_cleavage) + '-' + extract_cleavage_location(cleavage)
            cleavages_with_ranges.append(f'{cleavage_range}_{isoform}')
        else:
            location = extract_cleavage_location(first_cleavage)
            cleavages_with_ranges.append(f'{location}_{isoform}')
        i += 1
    return cleavages_with_ranges

def parse_ranges(ranges_list):
    """Parses the ranges list and returns a list of ranges."""
    ranges = []
    for part in ranges_list:
        part = part.split('_')[0]
        if '-' in part:
            start, end = map(int, part.split('-'))
            ranges.append(list(range(start, end + 1)))
        else:
            ranges.append([int(part)])
    return ranges

def cleavage_score(ranges, cleavages):
    """Calculates the cleavage score for the given ranges and cleavages."""
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

def get_accession(accession: str, peptide: str, sorted_isoform_headers) -> Tuple[str, str, int, str]:
    """Get the isoform, sequence, offset and aligned sequence for the given accession and peptide."""
    offset = 0
    sequence = None
    aligned_sequence = ""
    for header in sorted_isoform_headers:
        if peptide in header[1]:
            isoform = header[0]
            sequence = header[1]
            offset = sequence.index(peptide)
            aligned_sequence = header[2]
            break
    if sequence is not None:
        return isoform, sequence, offset, aligned_sequence
    raise ValueError(f"Peptide {peptide} with accession {accession} not found in fasta file")

def count_missing_amino_acids(peptide: str, aligned_sequence: str, peptide_offset: int, exon_start_index: int, exon_end_index: int) -> int:
    """Count the missing amino acids in the aligned sequence."""
    missing = 0
    stop_count = False
    for i in range(peptide_offset):
        # -1 beacuse of 1 based index for exon_start_index and exon_end_index
        if exon_start_index-1 <= i < exon_end_index:
            continue
        if aligned_sequence[i] == '-':
            missing += 1

    j = 0
    for c in aligned_sequence[peptide_offset-1:]:
        stop_count = exon_start_index-1 <= i < exon_end_index
        if c == '-':
            if not stop_count:
                missing += 1
        elif peptide[j] == c:
            j+=1
        if j == len(peptide):
            break
    return missing

def write_results(all_mods, mods_for_exp, cleavages_with_ranges, cleavages_for_exp, output_folder, groups_df):
    """Write modification and cleavage strings to csv files."""
    with open(f"{output_folder}/result_max_quant_mods.csv", 'w', newline='', encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Group'] + all_mods)
        writer.writerow(['', ''] + [mod.split('(')[0] for mod in all_mods])
        writer.writerow(['', ''] + [extract_mod_location(mod) for mod in all_mods])
        writer.writerow(['', ''] + [mod.split('_')[1] for mod in all_mods])
        for key, value in mods_for_exp.items():
            row = [1 if mod in value else 0 for mod in all_mods]
            if key not in groups_df['file_name'].values:
                raise ValueError(f"File {key} not found in groups file")
            group = groups_df.loc[groups_df['file_name'] == key]['group_name'].values[0]
            writer.writerow([key, group] + row)

    with open(f"{output_folder}/result_max_quant_cleavages.csv", 'w', newline='', encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Group'] + cleavages_with_ranges)
        writer.writerow(['', ''] + ['Non-Tryptic' for _ in cleavages_with_ranges])
        writer.writerow(['', ''] + [cleavage.split('_')[0] for cleavage in cleavages_with_ranges])
        writer.writerow(['', ''] + [cleavage.split('_')[1] for cleavage in cleavages_with_ranges])
        ranges = parse_ranges(cleavages_with_ranges)
        for key, value in cleavages_for_exp.items():
            indexes = [extract_index(cleavage) for cleavage in value]
            row = cleavage_score(ranges, indexes)
            if key not in groups_df['file_name'].values:
                raise ValueError(f"File {key} not found in groups file")
            group = groups_df.loc[groups_df['file_name'] == key]['group_name'].values[0]
            writer.writerow([key, group] + row)

def calculate_exon_offset(offset: int, isoform: str, exon_found: bool, exon_end_index: int, exon_1_isoforms: list, exon_2_isoforms: list, exon_1_length: int, exon_2_length: int, exon_length: int) -> int:
    """Calculate the exon offset. Starting index is 1."""
    if exon_found:
        if isoform in exon_1_isoforms:
            if exon_1_length < exon_length and offset > exon_end_index:
                return offset + exon_length - exon_1_length
            return offset
        elif isoform in exon_2_isoforms:
            if exon_2_length < exon_length and offset > exon_end_index:
                return offset + exon_length - exon_2_length
            return offset
        else:
            if offset < exon_end_index - exon_length + 1:
                return offset
            return offset + exon_length
    else:
        return offset

def count_missing_aa_in_exon(aligned_sequence: str, exon_start_index: int, exon_end_index: int, offset: int) -> int:
    """Count the missing amino acids in the exon."""
    missing = 0
    for i in range(offset):
        if i >= exon_start_index-1 and i < exon_end_index:
            if aligned_sequence[i] == '-':
                missing += 1
    return missing

def check_N_term_cleavage(peptide: str, accession: str, sorted_isoform_headers, exon_found: bool, exon_start_index: int, exon_end_index: int, exon_1_isoforms: list, exon_2_isoforms: list, exon_1_length: int, exon_2_length: int, exon_length: int) -> str:
    """Check if the N-term cleavage is possible for the given peptide and accession."""
    isoform, sequence, offset, aligned_sequence = get_accession(accession, peptide, sorted_isoform_headers)
    amino_acid_first = peptide[0]
    amino_acid_before = ""
    missing_aa = 0
    if len(sequence) != len(aligned_sequence):
        missing_aa = count_missing_amino_acids(peptide, aligned_sequence, offset, exon_start_index, exon_end_index)

    if offset > 0:
        amino_acid_before = sequence[offset - 1]
    if amino_acid_before != "K" and amino_acid_before != "R":
        offset = calculate_exon_offset(offset+missing_aa, isoform, exon_found, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
        iso = get_isoform_for_offset(isoform, offset, exon_start_index, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length)
        return f"{amino_acid_first}@{offset+1}_{iso}"

    return ""

def check_C_term_cleavage(peptide: str, accession: str, sorted_isoform_headers,exon_found: bool, exon_start_index: int, exon_end_index: int, exon_1_isoforms: list, exon_2_isoforms: list, exon_1_length: int, exon_2_length: int, exon_length: int) -> str:
    """Check if the C-term cleavage is possible for the given peptide and accession."""
    isoform, sequence, offset, aligned_sequence = get_accession(accession, peptide, sorted_isoform_headers)
    missing_aa = 0
    if len(sequence) != len(aligned_sequence):
        missing_aa = count_missing_amino_acids(peptide, aligned_sequence, offset, exon_start_index, exon_end_index)
    amino_acid_last = peptide[-1]
    if amino_acid_last not in ["K", "R"]:
        offset = calculate_exon_offset(offset+len(peptide)+missing_aa, isoform, exon_found, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
        iso = get_isoform_for_offset(isoform, offset, exon_start_index, exon_1_isoforms, exon_1_length, exon_2_isoforms, exon_2_length)
        return f"{amino_acid_last}@{offset}_{iso}"

    return ""

def get_isoform_for_offset(isoform: str, offset: int, exon_start_index: int, exon_1_isoforms: list, exon_1_length: int, exon_2_isoforms: list, exon_2_length: int) -> str:
    """Get the isoform for the given offset."""
    iso = 'general'
    if isoform in exon_1_isoforms:
        if offset >= exon_start_index and offset <= exon_start_index + exon_1_length:
            iso = 'exon1'
    elif isoform in exon_2_isoforms:
        if offset >= exon_start_index and offset <= exon_start_index + exon_2_length:
            iso = 'exon2'
    return iso

def sort_by_index_and_exons(entries):
    """Sort the entries by index and exons."""
    before = []
    after = []
    exon1 = []
    exon2 = []
    exon = False
    for entry in entries:
        _, type_part = entry.split("_")
        if type_part == "general" and not exon:
            before.append(entry)
        elif type_part == "exon1":
            exon1.append(entry)
            exon = True
        elif type_part == "exon2":
            exon2.append(entry)
            exon = True
        else:
            after.append(entry)
    result = before + exon1 + exon2 + after
    return result
