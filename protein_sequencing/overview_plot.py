"""Module to generate overview plot for protein sequences."""
import importlib
from collections import defaultdict
import plotly.graph_objects as go
from protein_sequencing import utils, sequence_plot as sequence

class OverviewPlot:
    def __init__(self, config, plot_config, input_file, output_path):
        self.CONFIG = importlib.import_module(config, 'configs')
        self.PLOT_CONFIG = importlib.import_module(plot_config, 'configs')
        self.input_file = input_file
        self.output_path = output_path
        self.create_overview_plot()

    def get_present_modifications(self, mod_file):
        """Get present modifications"""
        with open(mod_file, 'r', encoding="utf-8") as f:
            rows = f.readlines()[1:4]
            modification_types = rows[0].strip().split(',')
            present_modifications = set()
            for i, (label) in enumerate(rows[1].strip().split(',')):
                if label == '':
                    continue
                aa = label[0]
                if modification_types[i] not in self.PLOT_CONFIG.MODIFICATIONS_GROUP:
                    continue
                if self.CONFIG.INCLUDED_MODIFICATIONS.get(modification_types[i]):
                    if aa not in self.CONFIG.INCLUDED_MODIFICATIONS[modification_types[i]]:
                        continue
                    if aa == 'R' and modification_types[i] == 'Deamidated':
                        modification_types[i] = 'Citrullination'
                present_modifications.add(modification_types[i])
        return present_modifications

    def get_modifications_per_position(self, mod_file):
        """Get modifications per amino acid position"""
        with open(mod_file, 'r', encoding='utf-8') as f:
            rows = f.readlines()[1:4]
            modification_types = rows[0].strip().split(',')
            labels = rows[1].strip().split(',')
            isoforms = rows[2].strip().split(',')
            modifications_by_position = defaultdict(list)
            for i, (label) in enumerate(labels):
                if label == '':
                    continue
                aa = label[0]
                if modification_types[i] not in self.PLOT_CONFIG.MODIFICATIONS_GROUP:
                    continue
                if self.CONFIG.INCLUDED_MODIFICATIONS.get(modification_types[i]):
                    if aa not in self.CONFIG.INCLUDED_MODIFICATIONS[modification_types[i]]:
                        continue
                    if aa == 'R' and modification_types[i] == 'Deamidated':
                        modification_types[i] = 'Citrullination'
                isoform = isoforms[i]
                position = utils.get_position_with_offset(int(label[1:]), isoform)
                modifications_by_position[position].append((label, modification_types[i], self.PLOT_CONFIG.MODIFICATIONS_GROUP[modification_types[i]], isoform))
            for position, mods in modifications_by_position.items():
                modifications_by_position[position] = list(set(mods))
        return modifications_by_position

    def plot_labels(self, fig, modifications_by_position):
        """Main plotting function. Plots labels for modifications at correspinding positions."""
        x0 = utils.SEQUENCE_BOUNDARIES['x0']
        x1 = utils.SEQUENCE_BOUNDARIES['x1']
        y0 = utils.SEQUENCE_BOUNDARIES['y0']
        y1 = utils.SEQUENCE_BOUNDARIES['y1']

        label_offsets_with_orientation = self.get_label_offsets_with_orientation(modifications_by_position)
        for aa_position in label_offsets_with_orientation.keys():
            line_plotted_a, line_plotted_b = False, False
            for height_offset, group, label, modification_type, orientation in label_offsets_with_orientation[aa_position]:
                if self.CONFIG.FIGURE_ORIENTATION == 0:
                    x_position_line = (aa_position * utils.PIXELS_PER_AA) + utils.SEQUENCE_OFFSET
                    x_position_line = utils.offset_line_for_exon(x_position_line, int(label[1:]), self.CONFIG.FIGURE_ORIENTATION)
                    y_length = self.PLOT_CONFIG.SEQUENCE_MIN_LINE_LENGTH + height_offset * utils.get_label_height()
                    y_beginning_line = y0 if group == 'B' else y1
                    y_end_line = y_beginning_line - y_length if group == 'B' else y_beginning_line + y_length

                    if not line_plotted_a and group == 'A':
                        self.plot_line(fig, x_position_line, x_position_line, y_beginning_line, y_end_line)
                        line_plotted_a = True
                    if not line_plotted_b and group == 'B':
                        self.plot_line(fig, x_position_line, x_position_line, y_beginning_line, y_end_line)
                        line_plotted_b = True

                    position_label = 'top center'
                    if group == 'B':
                        position_label = 'bottom '+orientation
                    if group == 'A':
                        position_label = 'top '+orientation

                    self.plot_label(fig, x_position_line, y_end_line, label, modification_type, position_label)
                else:
                    y_position_line = self.CONFIG.FIGURE_WIDTH - (aa_position * utils.PIXELS_PER_AA) - utils.SEQUENCE_OFFSET
                    y_position_line = utils.offset_line_for_exon(y_position_line, int(label[1:]), self.CONFIG.FIGURE_ORIENTATION)

                    x_length = self.PLOT_CONFIG.SEQUENCE_MIN_LINE_LENGTH + height_offset * utils.get_label_length(label)
                    x_beginning_line = x0 if group == 'B' else x1
                    x_end_line = x_beginning_line - x_length if group == 'B' else x_beginning_line + x_length

                    if not line_plotted_a and group == 'A':
                        self.plot_line(fig, x_beginning_line, x_end_line, y_position_line, y_position_line)
                        line_plotted_a = True
                    if not line_plotted_b and group == 'B':
                        self.plot_line(fig, x_beginning_line, x_end_line, y_position_line, y_position_line)
                        line_plotted_b = True

                    position_label = 'middle'
                    if orientation == 'left':
                        position_label = 'top'
                    if orientation == 'right':
                        position_label = 'bottom'
                    if group == 'B':
                        position_label = position_label + ' left'
                    if group == 'A':
                        position_label = position_label + ' right'

                    self.plot_label(fig, x_end_line, y_position_line, label, modification_type, position_label)
        return fig

    def get_distance_groups(self, group):
        """Get distance groups for modifications."""
        result = []
        last_sight = {'position': None, 'mod': None}
        distance_group = defaultdict(list)
        size = 0
        new_group = True
        for aa_pos in group.keys():
            for modification in group[aa_pos]:
                current_sight = {'position': aa_pos, 'mod': modification}
                if last_sight['position'] is not None:
                    if self.check_distance(last_sight, current_sight) <= 0:
                        if new_group:
                            distance_group[last_sight['position']].append(last_sight['mod'])
                            size += 1
                            new_group = False
                        distance_group[aa_pos].append(modification)
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

    def get_offsets_with_orientations(self, distance_group, label_offsets_with_orientation, group_label, nearest_left, nearest_right):
        """Get offsets with orientations for modifications."""
        nearest_left_offset = 0
        nearest_right_offset = 0
        if nearest_left[0] is not None:
            nearest_left_offset = nearest_left[1]+1
        if nearest_right[0] is not None:
            nearest_right_offset = nearest_right[1]+1
        n = distance_group[0]
        mid_count = (n - nearest_left[1] + nearest_right[1]) // 2
        mod_count = 0
        for position in distance_group[1].keys():
            mid_position = position
            if mod_count >= mid_count:
                break
            mod_count += len(distance_group[1][position])
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

    def find_nearest_positions(self, label_offsets_with_orientation, distance_group):
        """Find nearest positions for modifications."""
        first_position = min(distance_group[1].keys())
        last_position = max(distance_group[1].keys())

        nearest_smaller = None
        nearest_larger = None
        smaller_offset = 0
        larger_offset = 0

        for position in sorted((k for k in label_offsets_with_orientation.keys() if k < first_position), reverse=True):
            if position < first_position:
                distance = self.check_distance({'position': position, 'mod': (label_offsets_with_orientation[position][-1][2], label_offsets_with_orientation[position][-1][3])},
                                               {'position': first_position, 'mod': distance_group[1][first_position][0]})
                if distance < 2:
                    nearest_smaller = position
                    smaller_offset = label_offsets_with_orientation[position][-1][0]
                else:
                    break

        for position in sorted(k for k in label_offsets_with_orientation.keys() if k > last_position):
            if position > last_position:
                distance = self.check_distance({'position': position, 'mod': (label_offsets_with_orientation[position][-1][2], label_offsets_with_orientation[position][-1][3])},
                                               {'position': last_position, 'mod': distance_group[1][last_position][0]})
                if distance < 2:
                    nearest_larger = position
                    larger_offset = label_offsets_with_orientation[position][-1][0]
                else:
                    break

        return (nearest_smaller, smaller_offset), (nearest_larger, larger_offset)

    def get_label_offsets_with_orientation(self, groups_by_position_and_isoform):
        """Get label offsets with orientation for modifications."""
        group_a, group_b = utils.separate_by_group(groups_by_position_and_isoform)
        label_offsets_with_orientation_a = defaultdict(list)
        label_offsets_with_orientation_b = defaultdict(list)

        for group, label_offsets_with_orientation in [(group_a, label_offsets_with_orientation_a), (group_b, label_offsets_with_orientation_b)]:
            group_label = 'A' if group == group_a else 'B'
            if len(group) == 0:
                continue
            distance_groups = self.get_distance_groups(group)
            for distance_group in sorted(distance_groups, key=lambda x: x[0], reverse=True):
                nearest_left, nearest_right = self.find_nearest_positions(label_offsets_with_orientation, distance_group)
                self.get_offsets_with_orientations(distance_group, label_offsets_with_orientation, group_label, nearest_left, nearest_right)

        return {**label_offsets_with_orientation_a, **label_offsets_with_orientation_b}

    def check_distance(self, first_modification, second_modification):
        """Check distance between modifications.
        Returns:
            -1 if label must be positioned left or right
            0 if label must be positioned in the center
            1 if label can be positioned anywhere
            2 if there is enough space for both labels to be positioned left and right"""
        first_position = int(first_modification['position'])
        second_position = int(second_modification['position'])
        label_length = utils.get_label_length(first_modification['mod'][0]) if self.CONFIG.FIGURE_ORIENTATION == 0 else utils.get_label_height()
        distance_between_modifications = abs(first_position - second_position) * utils.PIXELS_PER_AA
        if distance_between_modifications < label_length/2:
            return -1
        if distance_between_modifications < label_length:
            return 0
        second_label_length = utils.get_label_length(second_modification['mod'][0]) if self.CONFIG.FIGURE_ORIENTATION == 0 else utils.get_label_height()
        if distance_between_modifications > label_length + second_label_length:
            return 2
        return 1

    def plot_line(self, fig, x_start, x_end, y_start, y_end):
        """Plot single line for modifications."""
        fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='black', width=1), showlegend=False, hoverinfo='none'))

    def plot_label(self, fig, x, y, text, modification_type, position_label):
        """Plots single label for modification."""
        #Label bounding box for highlitghted PTMs
        if f'{modification_type}({text[0]})@{text[1:]}' in self.CONFIG.PTMS_TO_HIGHLIGHT:
            x0 = x+1
            y0 = y-1
            x1 = x-utils.get_label_length(text)+2
            y1 = y+utils.get_label_height()-1
            if 'bottom' in position_label:
                y1 = y - utils.get_label_height() + 1
                y0 = y + 1
            if 'right' in position_label:
                x1 = x + utils.get_label_length(text)-2
                x0 = x - 1
            if 'center' in position_label:
                x1 = x + utils.get_label_length(text)//2
                x0 = x - utils.get_label_length(text)//2
            if 'middle' in position_label:
                y0 = y - utils.get_label_height()//2
                y1 = y + utils.get_label_height()//2

            fig.add_shape(
                    type="rect",
                    x0=x0,
                    y0=y0,
                    x1=x1,
                    y1=y1,
                    layer='below',
                    fillcolor=self.CONFIG.PTM_HIGHLIGHT_LABEL_COLOR,
                    line=dict(width=0),
                )
        fig.add_trace(go.Scatter(x=[x], y=[y], mode='text',
                                text=text,
                                textposition=position_label,
                                showlegend=False,
                                hoverinfo='none',
                                textfont=dict(
                                    family=self.CONFIG.FONT,
                                    size=self.CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                                    color=self.CONFIG.MODIFICATIONS[modification_type][1])))

    def create_overview_plot(self):
        """Create overview plot for protein sequences."""
        present_modifications = self.get_present_modifications(self.PLOT_CONFIG.INPUT_FILE)
        groups_present = {self.PLOT_CONFIG.MODIFICATIONS_GROUP[mod] for mod in present_modifications if mod in self.PLOT_CONFIG.MODIFICATIONS_GROUP}
        if not 'A' in groups_present:
            fig = sequence.create_plot(self.input_file, present_modifications, 'A', 'A')
        elif not 'B' in groups_present:
            fig = sequence.create_plot(self.input_file, present_modifications, 'B', 'B')
        else:
            fig = sequence.create_plot(self.input_file, present_modifications, None, 'A')

        modifications_by_position = self.get_modifications_per_position(self.PLOT_CONFIG.INPUT_FILE)
        fig = self.plot_labels(fig, modifications_by_position)

        utils.show_plot(fig, self.output_path)
        