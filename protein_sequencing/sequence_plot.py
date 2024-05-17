from collections import defaultdict
import time
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

    # TODO: needed later for different starts/ endings
    count, i = 0, 0
    while i < len(different_possibilities):
        if different_possibilities[i] == 2:
            count += 1
            while i < len(different_possibilities) and different_possibilities[i+1] == 2:
                i += 1
        i += 1

    # For debugging purposes
    # different_possibilities_plot(max_sequence_length, 50, different_possibilities)

    # basis for all pixel calculations
    margins = parameters.LEFT_MARGIN + parameters.RIGHT_MARGIN if parameters.FIGURE_ORIENTATION == 0 else parameters.TOP_MARGIN + parameters.BOTTOM_MARGIN
    max_sequence_length_pixels = parameters.FIGURE_WIDTH * (1 - margins)
    pixels_per_protein = int(max_sequence_length_pixels // max_sequence_length)

    sequence_length_pixels = max_sequence_length * pixels_per_protein

    # calculate region boundaries in pixels
    region_boundaries = []
    region_end_pixel = 0
    for region_name, region_end, region_color in parameters.REGIONS:
        region_start_pixel = region_end_pixel
        region_end_pixel = region_end * pixels_per_protein + 1
        region_boundaries.append((region_name, region_start_pixel, region_end_pixel, parameters.SEQUENCE_REGION_COLORS[region_color]))

    # TODO: create this file
    input_file = 'data/chris/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_reduced.csv'

    fig = create_sequence_plot(pixels_per_protein, parameters.SEQUENCE_PLOT_HEIGHT, region_boundaries, input_file)
    
    fig.show()
    fig.write_image(f'{output_path}/fig1.png')

    clean_up()

    return output_file

def different_possibilities_plot(width: int, height: int, different_possibilities: list[int]):
    rectangle = np.zeros((height, width))
    for i, value in enumerate(different_possibilities):        
        rectangle[:, i] = value
    fig = go.Figure(data=go.Heatmap(z=rectangle))
    fig.show()

def create_sequence_plot(pixels_per_protein: int, sequence_height: int, region_boundaries: list[tuple[str, int, int, str]], file_path: str) -> go.Figure:
    fig = go.Figure()
    
    width = parameters.FIGURE_WIDTH
    height = parameters.FIGURE_HEIGHT

    if parameters.FIGURE_ORIENTATION == 1:
        width, height = height, width

    left_margin = parameters.LEFT_MARGIN * width
    top_margin = parameters.TOP_MARGIN * height

    # General Layout
    fig.update_layout(
        title="Plot",
        width = width,
        height = height,
        xaxis=dict(range=[0, width], autorange=False),
        yaxis=dict(range=[0, height], autorange=False),
        plot_bgcolor="white"
    )

    for region_name, region_start_pixel, region_end_pixel, region_color in region_boundaries:
        
        x0 = region_start_pixel
        x1 = region_end_pixel
        y0 = 0
        y1 = sequence_height
        
        if parameters.FIGURE_ORIENTATION == 0:
            y0 = height/2 - sequence_height/2
            y1 += y0
            x1 += left_margin
            x0 += left_margin
        else:
            y0 = width/2 - sequence_height/2
            y1 += y0
            x0, x1, y0, y1 = y0, y1, height-x0, height-x1
            y0 -= top_margin
            y1 -= top_margin


        # Region rects
        fig.add_shape(
            type="rect",
            x0=x0,
            y0=y0,
            x1=x1,
            y1=y1,
            line=dict(color="darkgrey", width=2),
            fillcolor=region_color
        )

        # Labels
        x = (x0 + x1) / 2
        y = (y0 + y1) / 2
        fig.add_annotation(
            x=x,
            y=y,
            text=region_name,
            showarrow=False,
            font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="black"),
            textangle= 90 if parameters.FIGURE_ORIENTATION == 1 else 0
        )

    # Protein Annotations
    with open(file_path, 'r') as f:

        rows = f.readlines()[1:3]
        modification_types = rows[0].strip().split(',')
        labels = rows[1].strip().split(',')

        groups_by_position = defaultdict(list)
        for i, (label) in enumerate(labels):
            if label == '':
                continue
            position = int(label[1:])
            letter = label[0]
            if letter in parameters.EXCLUDED_MODIFICATIONS:
                continue
            groups_by_position[position].append((label, modification_types[i]))
        
        for protein_position, modification_types in groups_by_position.items():
            line_plotted_A, line_plotted_B = False, False
            vertical_text_offset = 0
            for i, modification_type in enumerate(modification_types):
                if modification_type[1] not in parameters.MODIFICATIONS:
                    continue
                
                if parameters.FIGURE_ORIENTATION == 0:
                    x_start = (protein_position * pixels_per_protein) + left_margin
                    x_end = x_start
                    y_start = y1
                    y_end = y_start + parameters.SEQUENCE_MIN_LINE_LENGTH
                else:
                    x_start = x1
                    x_end = x_start + parameters.SEQUENCE_MIN_LINE_LENGTH
                    y_start = height - (protein_position * pixels_per_protein) - top_margin
                    y_end = y_start
                    vertical_text_offset = parameters.FONT_SIZE/2 * len(modification_type[0])

                if parameters.MODIFICATIONS[modification_type[1]][2] == 'A' and parameters.FIGURE_ORIENTATION == 1:
                    x_start = x0
                    x_end = x_start - parameters.SEQUENCE_MIN_LINE_LENGTH
                if parameters.MODIFICATIONS[modification_type[1]][2] == 'B' and parameters.FIGURE_ORIENTATION == 0:
                    y_start = y0
                    y_end = y_start - parameters.SEQUENCE_MIN_LINE_LENGTH

                if not line_plotted_A:
                    fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='black', width=1), showlegend=False, hoverinfo='none'))  
                    line_plotted_A = True

                if not line_plotted_B:
                    fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='black', width=1), showlegend=False, hoverinfo='none'))  
                    line_plotted_B = True

                if parameters.FIGURE_ORIENTATION == 0:
                    fig.add_annotation(
                        x=x_end,
                        y=y_end+5 if parameters.MODIFICATIONS[modification_type[1]][2] == 'A' else y_end-5,
                        text=modification_type[0],
                        showarrow=False,
                        font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color=parameters.MODIFICATIONS[modification_type[1]][1]),
                    )
                else:
                    fig.add_annotation(
                        x=x_end-vertical_text_offset-5 if parameters.MODIFICATIONS[modification_type[1]][2] == 'A' else x_end+vertical_text_offset+5,
                        y=y_end,
                        text=modification_type[0],
                        showarrow=False,
                        font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color=parameters.MODIFICATIONS[modification_type[1]][1]),
                    )

    return fig

def clean_up():
    directory = 'data/tmp'

    files = os.listdir(directory)

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

create_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)