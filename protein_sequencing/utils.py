from collections import defaultdict
import os
import plotly.graph_objects as go

import numpy as np
from protein_sequencing import parameters

SEQUENCE_BOUNDARIES = (None, None, None, None)
PIXELS_PER_PROTEIN = 0
SEQUENCE_OFFSET = 0

def get_width():
    if parameters.FIGURE_ORIENTATION == 0:
        return parameters.FIGURE_WIDTH
    return parameters.FIGURE_HEIGHT

def get_height():
    if parameters.FIGURE_ORIENTATION == 0:
        return parameters.FIGURE_HEIGHT
    return parameters.FIGURE_WIDTH

def get_left_margin():
    return parameters.LEFT_MARGIN * get_width()

def get_top_margin(): 
    return parameters.TOP_MARGIN * get_height()

def get_right_margin():
    return parameters.RIGHT_MARGIN * get_width()

def get_bottom_margin():
    return parameters.BOTTOM_MARGIN * get_height()

def get_label_length(label):
    if parameters.FIGURE_ORIENTATION == 1:
        return parameters.FONT_SIZE/1.5 * len(label) + 4
    return parameters.FONT_SIZE/1.5 * len(label)

def get_label_height():
    return parameters.FONT_SIZE+parameters.FONT_SIZE/4

def separate_by_group(groups_by_position):
    group_a = defaultdict(list)
    group_b = defaultdict(list)

    for protein_position in groups_by_position.keys():
        modification_sights = groups_by_position[protein_position]
        for modification_sight in modification_sights:
            if modification_sight[2] == 'A':  # group A
                group_a[protein_position].append(modification_sight)
            else:  # group B
                group_b[protein_position].append(modification_sight)
    
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
    output_file = f"{output_path}/figure1.png"
    fig.show()
    fig.write_image(output_file)

    clean_up()

    return output_file