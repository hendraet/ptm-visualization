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

    fig, x0, x1, y0, y1 = plot_regions(fig, region_boundaries, sequence_height, width, height, left_margin, top_margin)
    fig = plot_labels(fig, pixels_per_protein, file_path, left_margin, top_margin, x0, x1, y0, y1)

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
    return fig, x0, x1, y0, y1

def plot_labels(fig, pixels_per_protein, file_path, left_margin, top_margin, x0, x1, y0, y1):
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
                if parameters.EXCLUDED_MODIFICATIONS[letter] is None:
                    continue
                if modification_types[i] in parameters.EXCLUDED_MODIFICATIONS[letter]:
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
                    y_position_line = parameters.FIGURE_WIDTH - (protein_position * pixels_per_protein) - top_margin

                    x_length = parameters.SEQUENCE_MIN_LINE_LENGTH + height_offset * get_label_length(label)
                    x_beginning_line = x0 if group == 'B' else x1
                    x_end_line = x_beginning_line - x_length if group == 'B' else x_beginning_line + x_length

                    if not line_plotted_A and group == 'A':
                        plot_line(fig, x_beginning_line, x_end_line, y_position_line, y_position_line)
                        line_plotted_A = True
                    if not line_plotted_B and group == 'B':
                        plot_line(fig, x_beginning_line, x_end_line, y_position_line, y_position_line)
                        line_plotted_B = True

                    position_label = 'middle'
                    if orientation == 'left':
                        position_label = 'top'
                    if orientation == 'right':
                        position_label = 'bottom'
                    if group == 'B':
                        position_label = position_label + ' left'
                    if group == 'A':
                        position_label = position_label + ' right'

                    plot_label(fig, x_end_line, y_position_line, label, modification_type, position_label)
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
    result = []
    last_sight = {'position': None, 'mod': None}
    distance_group = defaultdict(list)
    size = 0
    new_group = True
    for protein_position in group.keys():
        for modification in group[protein_position]:
            current_sight = {'position': protein_position, 'mod': modification}
            if last_sight['position'] is not None:
                if check_distance(last_sight, current_sight, pixels_per_protein) <= 0:
                    if new_group:
                        distance_group[last_sight['position']].append(last_sight['mod'])
                        size += 1
                        new_group = False
                    distance_group[protein_position].append(modification)
                    size += 1
                else:
                    if new_group:
                        distance_group[last_sight['position']].append(last_sight['mod'])
                        size += 1
                    result.append((size, distance_group))
                    new_group = True
                    size = 0
                    distance_group = defaultdict(list)
            last_sight = current_sight
    if not new_group:
        result.append((size, distance_group))
    else:
        result.append((1, {last_sight['position']: [last_sight['mod']]}))
    return result

# TODO refactor
def get_offsets_with_orientations(distance_group, label_offsets_with_orientation, group_label, nearest_left, nearest_right):
    nearest_left_offset = 0
    nearest_right_offset = 0
    if nearest_left[0] is not None:
        nearest_left_offset = nearest_left[1]+1
    if nearest_right[0] is not None:
        nearest_right_offset = nearest_right[1]+1
    n = distance_group[0]
    mid_count = (n - nearest_left[1] + nearest_right[1]) // 2
    sum = 0
    for position in distance_group[1].keys():
        mid_position = position
        if sum >= mid_count:
            break
        sum += len(distance_group[1][position])
    left_offset, right_offset = -1, -1

    additional_offset = 0
    for i, position in enumerate(sorted((k for k in distance_group[1].keys() if k < mid_position))):
        for mod in distance_group[1][position]:
            offset = i+additional_offset+nearest_left_offset
            orientation = 'left'
            additional_offset += 1
            label_offsets_with_orientation[position].append((offset, group_label, mod[0], mod[1], orientation))
            left_offset = max(left_offset, offset)
        additional_offset -= 1
    
    additional_offset = 0
    for i, position in enumerate(sorted((k for k in distance_group[1].keys() if k > mid_position), reverse=True)):
        for mod in distance_group[1][position]:
            offset = i + additional_offset + nearest_right_offset
            orientation = 'right'
            additional_offset += 1
            label_offsets_with_orientation[position].append((offset, group_label, mod[0], mod[1], orientation))
            right_offset = max(right_offset, offset)
        additional_offset -= 1
        
    additional_offset = 0
    if nearest_left_offset > 0:
        left_offset = max(left_offset, nearest_left_offset-1)
    if nearest_right_offset > 0:
        right_offset = max(right_offset, nearest_right_offset-1)
    for mod in distance_group[1][mid_position]:
        if right_offset < left_offset:
            offset = right_offset + additional_offset + 1
            orientation = 'right'
        elif right_offset > left_offset:
            offset = left_offset + additional_offset + 1
            orientation = 'left'
        else:
            orientation = 'center'
            offsets = []
            if right_offset != -1:
                offsets.append(right_offset)
            if left_offset != -1:
                offsets.append(left_offset)
            if len(offsets) == 0:
                offset = 0 + additional_offset
            else:
                offset = min(offsets) + additional_offset + 1
        additional_offset += 1
        label_offsets_with_orientation[mid_position].append((offset, group_label, mod[0], mod[1], orientation))        
    
    return label_offsets_with_orientation

def find_nearest_positions(label_offsets_with_orientation, distance_group, pixels_per_protein):
    first_position = min(distance_group[1].keys())
    last_position = max(distance_group[1].keys())

    nearest_smaller = None
    nearest_larger = None
    smaller_offset = 0
    larger_offset = 0

    for position in sorted((k for k in label_offsets_with_orientation.keys() if k < first_position), reverse=True):
        if position < first_position:
            distance = check_distance({'position': position, 'mod': (label_offsets_with_orientation[position][-1][2], label_offsets_with_orientation[position][-1][3])},
                                      {'position': first_position, 'mod': distance_group[1][first_position][0]},
                                      pixels_per_protein)
            if distance < 2:
                nearest_smaller = position
                smaller_offset = label_offsets_with_orientation[position][-1][0]
            else:
                break 

    for position in sorted(k for k in label_offsets_with_orientation.keys() if k > last_position):
        if position > last_position:
            distance = check_distance({'position': position, 'mod': (label_offsets_with_orientation[position][-1][2], label_offsets_with_orientation[position][-1][3])},
                                      {'position': last_position, 'mod': distance_group[1][last_position][0]},
                                      pixels_per_protein)
            if distance < 2:
                nearest_larger = position
                larger_offset = label_offsets_with_orientation[position][-1][0]
            else:
                break

    return (nearest_smaller, smaller_offset), (nearest_larger, larger_offset)

def get_label_offsets_with_orientation(groups_by_position, pixels_per_protein):
    group_a, group_b = separate_by_group(groups_by_position)
    label_offsets_with_orientation_a = defaultdict(list)
    label_offsets_with_orientation_b = defaultdict(list)

    for group, label_offsets_with_orientation in [(group_a, label_offsets_with_orientation_a), (group_b, label_offsets_with_orientation_b)]:
        group_label = 'A' if group == group_a else 'B'
        distance_groups = get_distance_groups(group, pixels_per_protein)
        for distance_group in sorted(distance_groups, key=lambda x: x[0], reverse=True):
            nearest_left, nearest_right = find_nearest_positions(label_offsets_with_orientation, distance_group, pixels_per_protein)
            get_offsets_with_orientations(distance_group, label_offsets_with_orientation, group_label, nearest_left, nearest_right)

    return {**label_offsets_with_orientation_a, **label_offsets_with_orientation_b}

# distance between positions, -1 if label must be positioned left or right, 0 if label must be positioned in the center, 1 if label can be positioned anywhere, 2 if there is enogh space for both labels to be positioned left and right
def check_distance(first_modification, second_modification, pixels_per_protein):
    label_length = get_label_length(first_modification['mod'][0]) if parameters.FIGURE_ORIENTATION == 0 else get_label_height()
    distance_between_modifications = abs(first_modification['position'] - second_modification['position']) * pixels_per_protein
    if distance_between_modifications < label_length/2:
        return -1
    if distance_between_modifications < label_length:
        return 0
    second_label_length = get_label_length(second_modification['mod'][0]) if parameters.FIGURE_ORIENTATION == 0 else get_label_height()
    if distance_between_modifications > label_length + second_label_length:
        return 2
    return 1

def get_label_length(label):
    return parameters.FONT_SIZE/1.5 * len(label)+2

def get_label_height():
    return parameters.FONT_SIZE+parameters.FONT_SIZE/4+2

def plot_line(fig, x_start, x_end, y_start, y_end):
    fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='black', width=1), showlegend=False, hoverinfo='none'))

def plot_label(fig, x, y, text, modification_type, position_label):
    # Label bounding box for debugging purposes
    # x0 = x
    # y0 = y
    # x1 = x-get_label_length(text)
    # y1 = y+get_label_height()
    # if 'bottom' in position_label:
    #     y1 = y - get_label_height()
    # if 'right' in position_label:
    #     x1 = x + get_label_length(text)
    # if 'center' in position_label:
    #     x1 = x + get_label_length(text)/2
    #     x0 = x - get_label_length(text)/2
    # if 'middle' in position_label:
    #     y0 = y - get_label_height()/2
    #     y1 = y + get_label_height()/2

    # fig.add_shape(
    #         type="rect",
    #         x0=x0,
    #         y0=y0,
    #         x1=x1,
    #         y1=y1,
    #         line=dict(color="red", width=1),
    #     )
    fig.add_trace(go.Scatter(x=[x], y=[y], mode='text', text=text, textposition=position_label, showlegend=False, hoverinfo='none', textfont=dict(size=parameters.SEQUENCE_PLOT_FONT_SIZE, color=parameters.MODIFICATIONS[modification_type][1])))

def clean_up():
    directory = 'data/tmp'

    files = os.listdir(directory)

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)

create_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)