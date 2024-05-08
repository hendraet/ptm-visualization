import plotly.graph_objects as go
import os
import uniprot_align
import parameters
import numpy as np

def create_plot(input_file: str | os.PathLike, output_path: str | os.PathLike) -> str:
    output_file = f"{output_path}_figure.png"

    alignments = list(uniprot_align.get_alignment(input_file))
    max_length = 0
    for alignment in alignments:
        if len(alignment.seq) > max_length:
            max_length = len(alignment.seq)
        # print(alignment.seq)

    assert all(len(alignment.seq) == max_length for alignment in alignments)

    different_possibilities = [-1]*max_length
    for i in range(len(alignments[0].seq)):
        proteins = set()
        for alignment in alignments:
            protein = alignment.seq[i]
            proteins.add(protein)
        
        if '-' in proteins:
            if len(proteins) == 2:
                different_possibilities[i] = -1
            if len(proteins) > 2:
                different_possibilities[i] = len(proteins)-1
        else:
            different_possibilities[i] = len(proteins)

    print(max_length)
    print(different_possibilities)

    width = max_length
    height = 50
    rectangle = np.zeros((height, width))
    for i, value in enumerate(different_possibilities):        
        rectangle[:, i] = value

    fig = go.Figure(data=go.Heatmap(z=rectangle))

    fig.show()
    fig.write_image(f'{output_path}/fig1.png')

    clean_up()

    return output_file

def clean_up():
    directory = 'data/tmp'

    files = os.listdir(directory)

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

create_plot('data/uniprot_data/tau_isoforms2N4R.fasta', 'output')