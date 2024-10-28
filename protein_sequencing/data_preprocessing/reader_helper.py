from typing import Tuple
from protein_sequencing import uniprot_align

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


def extract_index(string):
    return int(string.split('@')[1].split('_')[0])

def extract_mod_location(mod_string):
    return mod_string.split('(')[1].split(')')[0]+mod_string.split('@')[1].split('_')[0]

def extract_cleavage_location(cleavage_string):
    return cleavage_string.split('@')[1].split('_')[0]

def extract_cleavages_ranges(all_cleavages):
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
    offset = 0
    sequence = None
    for header in sorted_isoform_headers:
        if peptide in header[1]:
            isoform = header[0]
            sequence = header[1]
            offset = sequence.index(peptide)
            break
    if sequence is not None:
        return isoform, sequence, offset, header[2]
    else:
        raise ValueError(f"Peptide {peptide} with accession {accession} not found in fasta file")
    

def count_missing_amino_acids(peptide: str, aligned_sequence: str, peptide_offset: int, exon_start_index: int, exon_end_index: int) -> int:
    missing = 0
    stop_count = False
    for i in range(peptide_offset):
        # -1 beacuse of 1 based index for exon_start_index and exon_end_index
        if i >= exon_start_index-1 and i < exon_end_index:
            continue
        if aligned_sequence[i] == '-':
            missing += 1

    j = 0
    for i in range(peptide_offset, len(aligned_sequence)):
        if i >= exon_start_index-1 and i < exon_end_index:
            stop_count = True
        else:
            stop_count = False
        if aligned_sequence[i] == '-':
            if not stop_count:
                missing += 1
        elif peptide[j] == aligned_sequence[i]:
            j+=1
        if j == len(peptide):
            break
    return missing

# calculates with 1 based index
def calculate_exon_offset(offset: int, isoform: str, exon_found: bool, exon_end_index: int, exon_1_isoforms: list, exon_2_isoforms: list, exon_1_length: int, exon_2_length: int, exon_length: int) -> int:
    if exon_found and offset > exon_end_index:
        if isoform in exon_1_isoforms:
            return offset - exon_1_length + exon_length
        elif isoform in exon_2_isoforms:
            return offset - exon_2_length + exon_length
        else:
            return offset + exon_length
    else:
        return offset
    
def check_N_term_cleavage(peptide: str, accession: str, sorted_isoform_headers, exon_found: bool, exon_start_index: int, exon_end_index: int, exon_1_isoforms: list, exon_2_isoforms: list, exon_1_length: int, exon_2_length: int, exon_length: int) -> str:
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
        return f"{amino_acid_first}@{offset+1}_{isoform}"

    return ""

def check_C_term_cleavage(peptide: str, accession: str, sorted_isoform_headers,exon_found: bool, exon_start_index: int, exon_end_index: int, exon_1_isoforms: list, exon_2_isoforms: list, exon_1_length: int, exon_2_length: int, exon_length: int) -> str:
    isoform, sequence, offset, aligned_sequence = get_accession(accession, peptide, sorted_isoform_headers)
    missing_aa = 0
    if len(sequence) != len(aligned_sequence):
        missing_aa = count_missing_amino_acids(peptide, aligned_sequence, offset, exon_start_index, exon_end_index)
    amino_acid_last = peptide[-1]
    if amino_acid_last not in ["K", "R"]:
        offset = calculate_exon_offset(offset+len(peptide)+missing_aa, isoform, exon_found, exon_end_index, exon_1_isoforms, exon_2_isoforms, exon_1_length, exon_2_length, exon_length)
        return f"{amino_acid_last}@{offset}_{isoform}"

    return ""