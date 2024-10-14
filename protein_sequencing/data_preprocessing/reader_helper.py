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