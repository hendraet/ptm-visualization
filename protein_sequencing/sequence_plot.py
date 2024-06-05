from collections import defaultdict
import time
import plotly.graph_objects as go
import os
from protein_sequencing import utils
import uniprot_align
import parameters
import numpy as np

def create_plot(input_file: str | os.PathLike) -> go.Figure:

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
    utils.PIXELS_PER_PROTEIN = int(max_sequence_length_pixels // max_sequence_length)

    sequence_length_pixels = max_sequence_length * utils.PIXELS_PER_PROTEIN

    # calculate region boundaries in pixels
    region_boundaries = []
    region_end_pixel = 0
    region_start = 1
    for region_name, region_end, region_color in parameters.REGIONS:
        region_start_pixel = region_end_pixel
        region_end_pixel = region_end * utils.PIXELS_PER_PROTEIN + 1
        region_boundaries.append((region_name, region_start_pixel, region_end_pixel, parameters.SEQUENCE_REGION_COLORS[region_color], region_start, region_end))
        region_start = region_end + 1

    fig = create_sequence_plot(parameters.SEQUENCE_PLOT_HEIGHT, region_boundaries)
    
    return fig

def create_sequence_plot(sequence_height: int, region_boundaries: list[tuple[str, int, int, str, int, int]]) -> go.Figure:
    fig = go.Figure()
    
    width = utils.get_width()
    height = utils.get_height()
    left_margin = utils.get_left_margin()
    top_margin = utils.get_top_margin()

    # General Layout
    fig.update_layout(
        title="",
        width = width,
        height = height,
        xaxis=dict(range=[0, width], autorange=False),
        yaxis=dict(range=[0, height], autorange=False),
        plot_bgcolor="white",
        font_family=parameters.FONT,
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)

    for i, modification in enumerate(parameters.MODIFICATIONS.values()):
        fig.add_trace(go.Scatter(x=[0], y=[height - i*utils.get_label_height()], mode='text', text=modification[0], textposition="bottom right", showlegend=False, hoverinfo='none', textfont=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color=modification[1])))

    fig = plot_regions(fig, region_boundaries, sequence_height, width, height, left_margin, top_margin)

    return fig

def plot_regions(fig, region_boundaries, sequence_height, width, height, left_margin, top_margin):
    for i, (region_name, region_start_pixel, region_end_pixel, region_color, region_start, region_end) in enumerate(region_boundaries):
        
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

        if i == 0:
            if parameters.FIGURE_ORIENTATION == 0:
                x = x0 - utils.get_label_length(str(region_start))
                y = y
            else:
                x = x
                y = y0 + utils.get_label_height()
            fig.add_annotation(
                x=x,
                y=y,
                text='1',
                showarrow=False,
                font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="gray"),
                textangle= 0
            )
    if parameters.FIGURE_ORIENTATION == 0:
        x = x1 + utils.get_label_length(str(region_end))
        y = y
    else:
        x = x
        y = y1 - utils.get_label_height()
    fig.add_annotation(
        x=x,
        y=y,
        text= region_end,
        showarrow=False,
        font=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color="gray"),
        textangle= 0
    )

    utils.SEQUENCE_BOUNDARIES = (x0, x1, y0, y1)
    return fig
