"""Utility functions for protein sequencing tool."""

import os
import importlib
from collections import defaultdict
import plotly.graph_objects as go
import numpy as np

CONFIG = importlib.import_module('configs.default_config', 'configs')

# x0, x1, y0, y1
SEQUENCE_BOUNDARIES = {'x0': 0, 'x1': 0, 'y0': 0, 'y1': 0}
PIXELS_PER_AA = 0
SEQUENCE_OFFSET = 0
EXON_1_OFFSET = {'index_start': -1,
                'index_end': -1,
                'pixel_start': -1,
                'pixel_end': -1}
EXON_2_OFFSET = {'index_start': -1,
                'index_end': -1,
                'pixel_start': -1,
                'pixel_end': -1}

ISOFORM_IDS = []

def get_width():
    """Return width of the plot, based on user settings in default_config.py."""
    if CONFIG.FIGURE_ORIENTATION == 0:
        return CONFIG.FIGURE_WIDTH
    return CONFIG.FIGURE_HEIGHT

def get_height():
    """Return height of the plot, based on user settings in default_config.py."""
    if CONFIG.FIGURE_ORIENTATION == 0:
        return CONFIG.FIGURE_HEIGHT
    return CONFIG.FIGURE_WIDTH

def get_left_margin():
    """Return the left margin for the sequence plot.
    Calculated based on the longest label in the legend."""
    longest_text = CONFIG.MODIFICATION_LEGEND_TITLE
    for mod in CONFIG.MODIFICATIONS:
        if len(CONFIG.MODIFICATIONS[mod][0]) > len(longest_text):
            longest_text = CONFIG.MODIFICATIONS[mod][0]
    return int((get_label_length(longest_text) / get_width() * 1.05) * get_width())

def get_top_margin():
    """Return the top margin for the sequence plot.
    Calculated based on the number of modifications in the legend."""
    legend_height = (len(CONFIG.MODIFICATIONS)+3) * get_label_height()
    return int((legend_height / get_height() * 1.05) *get_height())

def get_label_length(label):
    """Approximate the length of a label in pixels based on font size and label length."""
    return int(CONFIG.FONT_SIZE/1.5 * len(label))

def get_label_height():
    """Approximate the height of a label in pixels based on font size."""
    return CONFIG.FONT_SIZE+CONFIG.FONT_SIZE//5

def separate_by_group(groups_by_position_and_isoform):
    """Separate the modification sights into two groups based on the user defined groups."""
    group_a = defaultdict(list)
    group_b = defaultdict(list)

    for key, value in groups_by_position_and_isoform.items():
        for modification_sight in value:
            if modification_sight[2] == 'A':  # group A
                group_a[key].append(modification_sight)
            else:  # group B
                group_b[key].append(modification_sight)

    return group_a, group_b

def different_possibilities_plot(width: int, height: int, different_possibilities: list[int]):
    """Debug option. Plot the different possibilities of the sequence in a heatmap."""
    rectangle = np.zeros((height, width))
    for i, value in enumerate(different_possibilities):
        rectangle[:, i] = value
    fig = go.Figure(data=go.Heatmap(z=rectangle))
    fig.show()

def clean_up():
    """Remove all files from the temporary directory."""
    directory = 'data/tmp'

    files = os.listdir(directory)

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

def show_plot(fig, output_path):
    """Show the plot and save it as a .png and .svg file."""
    # TODO: hardcoded paths -.-
    output_svg = f"{output_path}/figure1.svg"
    output_png = f"{output_path}/figure1.png"
    fig.show()
    fig.write_image(output_png)
    fig.write_image(output_svg)

    clean_up()

    return output_png

def get_position_with_offset(position, isoform):
    """Return the position in the rendering index based on sequence position and isoform."""
    if isoform == 'exon2':
        position += EXON_1_OFFSET['index_end'] - EXON_1_OFFSET['index_start'] + 1
    elif position > max(EXON_1_OFFSET['index_end'], EXON_2_OFFSET['index_end']):
        if isoform != 'general':
            raise ValueError(f"Position {position} is out of range for isoform {isoform}")
        exon_1_length = EXON_1_OFFSET['index_end'] - EXON_1_OFFSET['index_start'] + 1
        exon_2_length = EXON_2_OFFSET['index_end'] - EXON_2_OFFSET['index_start'] + 1
        position += max(exon_1_length, exon_2_length)

    return position

def offset_line_for_exon(line_position, aa_position, oritentation):
    """Offset the line position based on the exon boundaries."""
    if aa_position >= EXON_1_OFFSET['index_start'] and EXON_1_OFFSET['index_start'] != -1:
        if oritentation == 0:
            line_position += CONFIG.EXONS_GAP
        else:
            line_position -= CONFIG.EXONS_GAP
    if aa_position > EXON_1_OFFSET['index_end'] and EXON_1_OFFSET['index_start'] != -1:
        if oritentation == 0:
            line_position += CONFIG.EXONS_GAP
        else:
            line_position -= CONFIG.EXONS_GAP

    return line_position
