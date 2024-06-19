from collections import defaultdict
import os

import plotly.graph_objects as go
from protein_sequencing import parameters, utils, sequence_plot as sequence



def plot_labels(fig, file_path):
    x0 = utils.SEQUENCE_BOUNDARIES['x0']
    x1 = utils.SEQUENCE_BOUNDARIES['x1']
    y0 = utils.SEQUENCE_BOUNDARIES['y0']
    y1 = utils.SEQUENCE_BOUNDARIES['y1']
    modifications_by_position = utils.get_modifications_per_position(file_path)

    label_offsets_with_orientation = get_label_offsets_with_orientation(modifications_by_position, utils.PIXELS_PER_PROTEIN)
    for protein_position in label_offsets_with_orientation.keys():
        line_plotted_A, line_plotted_B = False, False
        for height_offset, group, label, modification_type, orientation in label_offsets_with_orientation[protein_position]:
            if parameters.FIGURE_ORIENTATION == 0:
                x_position_line = (protein_position * utils.PIXELS_PER_PROTEIN) + utils.SEQUENCE_OFFSET

                y_length = parameters.SEQUENCE_MIN_LINE_LENGTH + height_offset * utils.get_label_height()
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
                y_position_line = parameters.FIGURE_WIDTH - (protein_position * utils.PIXELS_PER_PROTEIN) - utils.SEQUENCE_OFFSET

                x_length = parameters.SEQUENCE_MIN_LINE_LENGTH + height_offset * utils.get_label_length(label)
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
    group_a, group_b = utils.separate_by_group(groups_by_position)
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
    label_length = utils.get_label_length(first_modification['mod'][0]) if parameters.FIGURE_ORIENTATION == 0 else utils.get_label_height()
    distance_between_modifications = abs(first_modification['position'] - second_modification['position']) * pixels_per_protein
    if distance_between_modifications < label_length/2:
        return -1
    if distance_between_modifications < label_length:
        return 0
    second_label_length = utils.get_label_length(second_modification['mod'][0]) if parameters.FIGURE_ORIENTATION == 0 else utils.get_label_height()
    if distance_between_modifications > label_length + second_label_length:
        return 2
    return 1

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
    fig.add_trace(go.Scatter(x=[x], y=[y], mode='text',
                             text=text,
                             textposition=position_label,
                             showlegend=False,
                             hoverinfo='none',
                             textfont=dict(
                                 family=parameters.FONT,
                                 size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                                 color=parameters.MODIFICATIONS[modification_type][1])))

def create_overview_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
    fig = sequence.create_plot(input_file)

    # TODO: create this file
    input_file = 'data/chris/overview_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_reduced.csv'
    fig = plot_labels(fig, input_file)

    utils.show_plot(fig, output_path)

create_overview_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)
    