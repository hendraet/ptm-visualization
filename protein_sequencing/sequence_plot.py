import plotly.graph_objects as go
import os
import uniprot_align
import parameters

def create_plot(input_file: str | os.PathLike, output_path: str | os.PathLike) -> str:
    output_file = f"{output_path}_figure.png"

    alignments = uniprot_align.align_protein(input_file)

    layout = go.Layout(
        font=dict(size=parameters.FONT_SIZE),
        width = parameters.FIGURE_WIDTH,
        height = parameters.FIGURE_HEIGHT
    )

    fig = go.Figure(data = [], layout=layout)
    

    fig.show()
    fig.write_image(f'{output_path}/fig1.png')

    return output_file

create_plot('data/uniprot_data/tau_isoforms2N4R.fasta', 'output')