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
        plot_bgcolor="white",
        font_family=parameters.FONT,
    )

    fig, y0, y1 = plot_regions(fig, region_boundaries, sequence_height, width, height, left_margin, top_margin)
    fig = plot_labels(fig, pixels_per_protein, file_path, left_margin, y0, y1)

    return fig

def plot_regions(fig, region_boundaries, sequence_height, width, height, left_margin, top_margin):
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
    return fig, y0, y1

def plot_labels(fig, pixels_per_protein, file_path, left_margin, y0, y1):
    # Protein Annotations
    with open(file_path, 'r') as f:

        rows = f.readlines()[1:3]
        modification_types = rows[0].strip().split(',')
        labels = rows[1].strip().split(',')

        modifications_by_position = defaultdict(list)
        for i, (label) in enumerate(labels):
            if label == '':
                continue
            position = int(label[1:])
            letter = label[0]
            if letter in parameters.EXCLUDED_MODIFICATIONS:
                continue
            if modification_types[i] not in parameters.MODIFICATIONS:
                continue
            modifications_by_position[position].append((label, modification_types[i], parameters.MODIFICATIONS[modification_types[i]][2]))
        for position, mods in modifications_by_position.items():
            modifications_by_position[position] = list(set(mods))

        label_offsets_with_orientation = get_label_offsets_with_orientation(modifications_by_position, pixels_per_protein)
        for protein_position in label_offsets_with_orientation.keys():
            line_plotted_A, line_plotted_B = False, False
            for height_offset, group, label, modification_type, orientation in label_offsets_with_orientation[protein_position]:
                if modification_type not in parameters.MODIFICATIONS:
                    continue

                if parameters.FIGURE_ORIENTATION == 0:
                    x_position_line = (protein_position * pixels_per_protein) + left_margin

                    y_length = parameters.SEQUENCE_MIN_LINE_LENGTH + height_offset * get_label_height()
                    y_beginning_line = y0 if group == 'B' else y1
                    y_end_line = y_beginning_line - y_length if group == 'B' else y_beginning_line + y_length
                    
                    if not line_plotted_A and group == 'A':
                        plot_line(fig, x_position_line, x_position_line, y_beginning_line, y_end_line)
                        line_plotted_A = True
                    if not line_plotted_B and group == 'B':
                        plot_line(fig, x_position_line, x_position_line, y_beginning_line, y_end_line)
                        line_plotted_B = True

                    position_label = 'top center'
                    if group == 'B':
                            position_label = 'bottom '+orientation
                    if group == 'A':
                        position_label = 'top '+orientation
                    
                    plot_label(fig, x_position_line, y_end_line, label, modification_type, position_label)
                else:
                    # TODO: implement vertical orientation
                    pass
    return fig

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

def get_distance_groups(group, pixels_per_protein):
    distance_groups = []
    last_sight = {'position': None, 'mod': None}
    distance_groups = [defaultdict(list)]
    new_group = True
    for protein_position in group.keys():
        for modification in group[protein_position]:
            current_sight = {'position': protein_position, 'mod': modification}
            if last_sight['position'] is not None:
                if check_distance(last_sight, current_sight, pixels_per_protein) == -1:
                    if new_group:
                        distance_groups[-1][last_sight['position']].append(last_sight['mod'])
                        new_group = False
                    distance_groups[-1][protein_position].append(modification)
                else:
                    if new_group:
                        distance_groups[-1][last_sight['position']].append(last_sight['mod'])
                    new_group = True
                    distance_groups.append(defaultdict(list))
            last_sight = current_sight
    return distance_groups

# TODO refactor
def get_offsets_with_orientations(distance_group, label_offsets_with_orientation, group_label):
    n = len(distance_group)
    mid = n // 2

    for i, position in enumerate(distance_group.keys()):
        mod = distance_group[position]

        if i < mid:
            # Entries from start to mid
            offset = i
            orientation = 'left'
        elif i > mid:
            # Entries from mid to end
            offset = n - 1 - i
            orientation = 'right'
        else:
            # Center entry (handle both even and uneven cases)
            offset = mid-1 if n % 2 == 0 else mid
            orientation = 'center' if n % 2 != 0 else ('right' if i == mid else 'left')

        label_offsets_with_orientation[position].append((offset, group_label, mod[0], mod[1], orientation))
    
    return label_offsets_with_orientation

def find_nearest_positions(label_offsets_with_orientation, distance_group):
    first_position = min(distance_group.keys())
    last_position = max(distance_group.keys())

    nearest_smaller = None
    nearest_larger = None
    min_distance_smaller = float('inf')
    min_distance_larger = float('inf')

    for position in label_offsets_with_orientation.keys():
        if position < first_position:
            distance = abs(position - first_position)
            if distance < min_distance_smaller:
                min_distance_smaller = distance
                nearest_smaller = position
        elif position > last_position:
            distance = abs(position - last_position)
            if distance < min_distance_larger:
                min_distance_larger = distance
                nearest_larger = position

    return nearest_smaller, nearest_larger

def get_label_offsets_with_orientation(groups_by_position, pixels_per_protein):
    group_a, group_b = separate_by_group(groups_by_position)
    label_offsets_with_orientation = defaultdict(list)

    for group in [group_a, group_b]:
        group_label = 'A' if group == group_a else 'B'
        group_distance_groups = get_distance_groups(group, pixels_per_protein)
        for distance_group in sorted(group_distance_groups, key=len, reverse=True):
            nearest_left, nearest_right = find_nearest_positions(label_offsets_with_orientation, distance_group)
            get_offsets_with_orientations(distance_group, label_offsets_with_orientation, group_label)
            # means label_offsets_with_orientation is empty
            # if not nearest_left and not nearest_right:
            #     get_offsets_with_orientations(distance_group, label_offsets_with_orientation, group_label)
            # else:
            #     if nearest_left:
            #         if check_distance({'position': nearest_left, 'mod': label_offsets_with_orientation[nearest_left][0]}, distance_group[0], pixels_per_protein) == 1:
            #             get_offsets_with_orientations(distance_group, label_offsets_with_orientation, group_label)
            #     if nearest_right:
            #         pass
       

    return label_offsets_with_orientation

# distance between positions, -1 if label must be positioned left or right, 0 if label must be positioned in the center, 1 if label can be positioned anywhere
def check_distance(first_modification, second_modification, pixels_per_protein):
    distance_between_modifications = abs(first_modification['position'] - second_modification['position']) * pixels_per_protein
    label_length = get_label_length(first_modification['mod'][0]) + get_label_length(second_modification['mod'][0])
    if distance_between_modifications < label_length/2:
        return -1
    if distance_between_modifications < label_length:
        return 0
    if distance_between_modifications > label_length*2:
        return 2
    return 1

def get_label_length(label):
    return parameters.FONT_SIZE/1.5 * len(label)

def get_label_height():
    return parameters.FONT_SIZE+parameters.FONT_SIZE/4

def plot_line(fig, x_start, x_end, y_start, y_end):
    fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='black', width=1), showlegend=False, hoverinfo='none'))

def plot_label(fig, x, y, text, modification_type, position_label):
    # Label bounding box for debugging purposes
    x0 = x
    x1 = x-get_label_length(text)
    y1 = y+get_label_height()
    if 'bottom' in position_label:
        y1 = y - get_label_height()
    if 'right' in position_label:
        x1 = x + get_label_length(text)
    if 'center' in position_label:
        x1 = x + get_label_length(text)/2
        x0 = x - get_label_length(text)/2
    fig.add_shape(
            type="rect",
            x0=x0,
            y0=y,
            x1=x1,
            y1=y1,
            line=dict(color="red", width=1),
        )
    fig.add_trace(go.Scatter(x=[x], y=[y], mode='text', text=text, textposition=position_label, showlegend=False, hoverinfo='none', textfont=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color=parameters.MODIFICATIONS[modification_type][1])))

def clean_up():
    directory = 'data/tmp'

    files = os.listdir(directory)

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

create_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)