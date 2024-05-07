from Bio import AlignIO, SeqIO, Seq, Align

import os


def align_protein(input_file: str | os.PathLike):
    # code from: https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length

    records = list(SeqIO.parse(input_file, 'fasta'))

    # get max sequence length
    max_length = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) < max_length:
            sequence = str(record.seq).ljust(max_length, '.')
            record.seq = Seq.Seq(sequence)

    # check that all sequences have the same length         
    assert all(len(record.seq) == max_length for record in records)

    # write to temporary file and do alignment
    output_file = f'{os.path.splitext(input_file)[0]}_padded.fasta'
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')

    alignment = AlignIO.parse(output_file, "fasta")
    os.remove(output_file)

    # alignment ready
    return alignment