import plotly.graph_objects as go
import os
import uniprot_align
import parameters
import numpy as np

def create_plot(input_file: str | os.PathLike, output_path: str | os.PathLike) -> str:
    output_file = f"{output_path}_figure.png"

    alignments = list(uniprot_align.get_alignment(input_file))
    max_sequence_length = 0
    for alignment in alignments:
        if len(alignment.seq) > max_sequence_length:
            max_sequence_length = len(alignment.seq)
        # print(alignment.seq)

    assert all(len(alignment.seq) == max_sequence_length for alignment in alignments)

    different_possibilities = [-1]*max_sequence_length
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

    # print(max_sequence_length)
    # print(different_possibilities)

    # TODO: needed later for different starts/ endings
    count, i = 0, 0
    while i < len(different_possibilities):
        if different_possibilities[i] == 2:
            count += 1
            while i < len(different_possibilities) and different_possibilities[i+1] == 2:
                i += 1
        i += 1

    # basis for all pixel calculations
    max_sequence_length_pixels = parameters.FIGURE_WIDTH * (1 - parameters.LEFT_MARGIN - parameters.RIGHT_MARGIN)
    pixels_per_protein = int(max_sequence_length_pixels // max_sequence_length)

    sequence_length_pixels = max_sequence_length * pixels_per_protein

    # calculate region boundaries in pixels
    region_boundaries = []
    region_end_pixel = 0
    for region_name, region_end, region_color in parameters.REGIONS:
        region_start_pixel = region_end_pixel
        region_end_pixel = region_end * pixels_per_protein + 1
        region_boundaries.append((region_name, region_start_pixel, region_end_pixel, region_color))

    width = sequence_length_pixels
    height = parameters.SEQUENCE_PLOT_HEIGHT

    fig = create_sequence_plot(width, height, region_boundaries)
    
    fig.show()
    fig.write_image(f'{output_path}/fig1.png')

    clean_up()

    return output_file

def create_sequence_plot(width: int, height: int, region_boundaries: list[tuple[str, int, int, str]]) -> go.Figure:
    fig = go.Figure()

    for region_name, region_start_pixel, region_end_pixel, region_color in region_boundaries:
        x = (region_start_pixel + region_end_pixel) / 2
        fig.add_shape(
            type="rect",
            x0=region_start_pixel,
            y0=0,
            x1=region_end_pixel,
            y1=height,
            line=dict(color="darkgrey", width=1),
            fillcolor=region_color
        )

        # Add region name in the middle of the region
        fig.add_annotation(
            x=x,
            y=height / 2,
            text=region_name,
            showarrow=False,
            font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="black")
        )

        fig.update_layout(
            title="Rectangle Divided into Regions",
            width = parameters.FIGURE_WIDTH,
            height = parameters.FIGURE_HEIGHT,
            xaxis=dict(range=[0, parameters.FIGURE_WIDTH], autorange=False),
            yaxis=dict(range=[0, parameters.FIGURE_HEIGHT], autorange=False),
            plot_bgcolor="white"
        )

    return fig

def clean_up():
    directory = 'data/tmp'

    files = os.listdir(directory)

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

create_plot('data/uniprot_data/tau_isoforms2N4R.fasta', 'output')