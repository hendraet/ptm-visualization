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

        groups_by_position = defaultdict(list)
        for i, (label) in enumerate(labels):
            if label == '':
                continue
            position = int(label[1:])
            letter = label[0]
            if letter in parameters.EXCLUDED_MODIFICATIONS:
                continue
            groups_by_position[position].append((label, modification_types[i]))
        for position, mods in groups_by_position.items():
            groups_by_position[position] = list(set(mods))

        label_offsets_with_orientation = get_label_offsets_with_orientation(groups_by_position, pixels_per_protein)
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

def separate_by_group(resulting_offsets_with_orientation):
    group_a = defaultdict(list)
    group_b = defaultdict(list)

    for protein_position, offsets in resulting_offsets_with_orientation.items():
        for offset in offsets:
            if offset[1] == 'A':  # group A
                group_a[protein_position].append(offset)
            else:  # group B
                group_b[protein_position].append(offset)
    
    return group_a, group_b

def check_and_adjust_overlaps(resulting_offsets_with_orientation, pixels_per_protein):
    def process_group(group):
        sorted_positions = sorted(group.keys())
        for i in range(len(sorted_positions)):
            current_position = sorted_positions[i]
            current_labels = group[current_position]

            # Only need to get the length of one label, since they are all the same at this position
            current_label_length = get_label_length(current_labels[0][2])
            additional_offset = 0
            
            # Compare with subsequent labels within potential overlap range
            for k in range(i + 1, len(sorted_positions)):
                compare_position = sorted_positions[k]
                
                # Check if the positions are within the overlap range
                if (compare_position - current_position) * pixels_per_protein < current_label_length + get_label_length(group[compare_position][0][2]):
                    if current_labels[0][4] == 'right' and group[compare_position][0][4] == 'left':
                        compare_labels = group[compare_position]
                        
                        # Adjust the height offset of all labels at the compare position uniformly
                        for l in range(len(compare_labels)):
                            compare_height_offset, compare_group, compare_label, compare_modification_type, compare_orientation = compare_labels[l]
                            new_height_offset = current_labels[-1][0]+additional_offset
                            group[compare_position][l] = (
                                new_height_offset, compare_group, compare_label, compare_modification_type, compare_orientation
                            )
                else:
                    # No need to check further if positions are beyond the overlap range
                    additional_offset = 0
                    break

        return group

    group_a, group_b = separate_by_group(resulting_offsets_with_orientation)

    adjusted_group_a = process_group(group_a)
    adjusted_group_b = process_group(group_b)

    # Combine adjusted groups back into one dictionary
    adjusted_offsets_with_orientation = {**adjusted_group_a, **adjusted_group_b}

    return adjusted_offsets_with_orientation

def get_label_offsets_with_orientation(groups_by_position, pixels_per_protein):
    # new algo idea:
    # calculate height offset for each label from smallest position to highest if labels placed left to the line
    # calculate height offset for each label from highest position to smallest if labels placed right to the line
    # retrieve minimal height offset for each label based on previous two steps, check if labels overlap due to left and right placement
    # check if label can be placed in the center between two lines, if yes do so
    # check if due to previous step surrounding labels can be dropped a layer
    # last iteration go through all labels, if left and right are smaller position label in the center

    # step 1 and 2
    height_offsets = calculate_height_offset(groups_by_position, pixels_per_protein, False)
    height_offsets_reversed = calculate_height_offset(groups_by_position, pixels_per_protein, True)

    # step 3.1
    resulting_offsets_with_orientation = defaultdict(list)
    for protein_position in height_offsets.keys():
        for i, (height_offset, group, label, modification_type) in enumerate(height_offsets[protein_position]):
            height_offset_reversed = height_offsets_reversed[protein_position][i][0]
            if height_offset_reversed < height_offset:
                resulting_offsets_with_orientation[protein_position].append((height_offset_reversed, group, label, modification_type, 'right'))
            else:
                resulting_offsets_with_orientation[protein_position].append((height_offset, group, label, modification_type, 'left'))

    # step 3.2
    adjusted_offsets = check_and_adjust_overlaps(resulting_offsets_with_orientation, pixels_per_protein)

    return adjusted_offsets

def calculate_height_offset(groups_by_position, pixels_per_protein, reversed):
    previous_labels = {'A': {}, 'B': {}}
    height_offsets = defaultdict(list)
    for protein_position in sorted(groups_by_position.keys(), reverse=reversed):
        modification_types = groups_by_position[protein_position]
        
        for modification_type in modification_types:
            if modification_type[1] not in parameters.MODIFICATIONS:
                continue
            # TODO check if this is correct for all font sizes

            horizontal_text_offset = get_label_length(modification_type[0])
            vertical_text_offset = get_label_height()
            height_offset = 0

            # Check for overlaps
            group = parameters.MODIFICATIONS[modification_type[1]][2]
            previous_label = previous_labels[group]
            if 'protein_position' in previous_label:
                if parameters.FIGURE_ORIENTATION == 0:
                    if abs(protein_position-previous_label['protein_position'])*pixels_per_protein < horizontal_text_offset:
                        height_offset = previous_label['height_offset'] + 1
                else:
                    # TODO implement vertical orientation
                    pass
            
            height_offsets[protein_position].append((height_offset, group, modification_type[0], modification_type[1]))
            previous_labels[group] = {'protein_position': protein_position, 'height_offset': height_offset}

    return height_offsets

def get_label_length(label):
    return parameters.FONT_SIZE/1.5 * len(label)

def get_label_height():
    return parameters.FONT_SIZE+parameters.FONT_SIZE/4

def plot_line(fig, x_start, x_end, y_start, y_end):
    fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='black', width=1), showlegend=False, hoverinfo='none'))

def plot_label(fig, x, y, text, modification_type, position_label):
    # Label bounding box for debugging purposes
    x1 = x-get_label_length(text)
    y1 = y+get_label_height()
    if 'bottom' in position_label:
        y1 = y - get_label_height()
    if 'right' in position_label:
        x1 = x + get_label_length(text)  
    fig.add_shape(
            type="rect",
            x0=x,
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