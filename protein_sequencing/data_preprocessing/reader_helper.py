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
        cleavage = all_cleavages[i]
        while i < len(all_cleavages)-1 and int(extract_cleavage_location(cleavage))+1 == int(extract_cleavage_location(all_cleavages[i+1])):
            i += 1
            cleavage = all_cleavages[i]
        if int(extract_cleavage_location(first_cleavage)) != int(extract_cleavage_location(cleavage)):
            cleavage_range = extract_cleavage_location(first_cleavage) + '-' + extract_cleavage_location(cleavage)
            cleavages_with_ranges.append(cleavage_range)
        else:
            cleavages_with_ranges.append(extract_cleavage_location(first_cleavage))
        i += 1
    return cleavages_with_ranges

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