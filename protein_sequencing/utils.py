from collections import defaultdict
import os
import plotly.graph_objects as go

import numpy as np
import importlib

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
    if CONFIG.FIGURE_ORIENTATION == 0:
        return CONFIG.FIGURE_WIDTH
    return CONFIG.FIGURE_HEIGHT

def get_height():
    if CONFIG.FIGURE_ORIENTATION == 0:
        return CONFIG.FIGURE_HEIGHT
    return CONFIG.FIGURE_WIDTH

def get_left_margin():
    longest_text = CONFIG.MODIFICATION_LEGEND_TITLE
    for mod in CONFIG.MODIFICATIONS:
        if len(CONFIG.MODIFICATIONS[mod][0]) > len(longest_text):
            longest_text = CONFIG.MODIFICATIONS[mod][0]
    return int((get_label_length(longest_text) / CONFIG.FIGURE_WIDTH * 1.05) * get_width())

def get_top_margin(): 
    return int(CONFIG.TOP_MARGIN * get_height())

def get_right_margin():
    return int(CONFIG.RIGHT_MARGIN * get_width())

def get_bottom_margin():
    return int(CONFIG.BOTTOM_MARGIN * get_height())

def get_label_length(label):
    return int(CONFIG.FONT_SIZE/1.5 * len(label))

def get_label_height():
    return CONFIG.FONT_SIZE+CONFIG.FONT_SIZE//5

def separate_by_group(groups_by_position_and_isoform):
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
    rectangle = np.zeros((height, width))
    for i, value in enumerate(different_possibilities):        
        rectangle[:, i] = value
    fig = go.Figure(data=go.Heatmap(z=rectangle))
    fig.show()

def clean_up():
    directory = 'data/tmp'

    files = os.listdir(directory)

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

def show_plot(fig, output_path):
    output_svg = f"{output_path}/figure1.svg"
    output_png = f"{output_path}/figure1.png"
    fig.show()
    fig.write_image(output_png)
    fig.write_image(output_svg)

    clean_up()

    return output_png

def get_position_with_offset(position, isoform):
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