"""Module to create bar plots for protein sequences"""
from collections import defaultdict
import plotly.graph_objects as go
import pandas as pd
from protein_sequencing import utils,sequence_plot

# TODO: linie mittig schwerpunk in plot

class BarPlot:
    def __init__(self, config, plot_config, input_file, output_path):
        self.CONFIG = config
        self.PLOT_CONFIG = plot_config
        self.input_file = input_file
        self.output_path = output_path
        # calculate exons and central sequence boundaries
        sequence_plot.create_plot(self.input_file, None, 'A')
        self.create_bar_plot()

    def get_bar_positions(self, modification_sights_all_a: dict[int, list[tuple[int, str, str]]], modification_sights_all_b: dict[int, list[tuple[int, str, str]]]):
        """Get positions for bar plot."""
        positions_a = []
        positions_b = []
        for aa_position in modification_sights_all_a.keys():
            for _ in modification_sights_all_a[aa_position]:
                positions_a.append(aa_position)
        for aa_position in modification_sights_all_b.keys():
            for _ in modification_sights_all_b[aa_position]:
                positions_b.append(aa_position)

        return positions_a, positions_b

    def get_bar_plot_width(self, group_size_a: int, group_size_b: int) -> int:
        """Get width of bar plot."""
        legend_height = utils.get_label_height() * (max(value.count('<br>') + 1 for value in self.PLOT_CONFIG.BAR_GROUPS.values()))
        if self.CONFIG.FIGURE_ORIENTATION == 0:
            legend_height += utils.get_label_length('100%')
            bar_width = (utils.get_width() - legend_height) // max(group_size_a, group_size_b)
            bar_plot_width = bar_width * max(group_size_a, group_size_b)
        else:
            legend_height += utils.get_label_height()
            bar_width = (utils.get_height() - legend_height) // max(group_size_a, group_size_b)
            bar_plot_width = bar_width * max(group_size_a, group_size_b)
        return bar_plot_width

    def add_bar_plot(self, fig: go.Figure, above: str, modification_sights_all: dict[int, list[tuple[int, str, str, str]]], modification_sights_relevant: dict[int, list[tuple[int, str, str, str]]], df, group_positions: list, bar_plot_width: int, label_plot_height: int) -> go.Figure:
        """Add bar plot to sequence plot."""
        group_direction = 1 if above == 'A' else -1
        bar_width = bar_plot_width // len(group_positions)
        assert bar_width >= self.CONFIG.FONT_SIZE, f"Too many bars to plot! Bar width: {bar_width} < FONT_SIZE: {self.CONFIG.FONT_SIZE}."
        bar_plot_margin = self.CONFIG.FONT_SIZE

        height_offset = 0
        modifications_visited = 0
        positions_visited = 0
        bar_percentages = {group: [] for group in self.PLOT_CONFIG.BAR_GROUPS.keys()}

        for aa_position in sorted(modification_sights_all.keys(), reverse=True):
            for modification_sight in modification_sights_all[aa_position]:
                if modification_sight not in modification_sights_relevant[aa_position]:
                    for group in self.PLOT_CONFIG.BAR_GROUPS.keys():
                        bar_percentages[group].append(0)
                    modifications_visited += 1
                    continue
                label, modification_type, _, isoform = modification_sight
                if self.CONFIG.FIGURE_ORIENTATION == 0:
                    # x position for protein sequence
                    position = utils.get_position_with_offset(aa_position, isoform)
                    x_0_line = position * utils.PIXELS_PER_AA + utils.SEQUENCE_OFFSET
                    x_0_line = utils.offset_line_for_exon(x_0_line, aa_position, self.CONFIG.FIGURE_ORIENTATION)
                    # x position for bar plot
                    x_1_line = utils.get_width() - (modifications_visited * bar_width + bar_width//2)
                    y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if above == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
                    y_1_line = y_0_line + group_direction * height_offset
                    y_3_line = y_0_line + label_plot_height * group_direction - (utils.get_label_length(label) + 10) * group_direction
                    y_2_line = y_3_line - 10 * group_direction
                    y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction
                    y_bar = y_3_line + (utils.get_label_length(label) + 5) * group_direction

                    # plot line with label
                    self.plot_line_with_label_horizontal(fig,
                                        x_0_line, x_1_line,
                                        y_0_line, y_1_line, y_2_line, y_3_line,
                                        y_label,
                                        self.CONFIG.MODIFICATIONS[modification_type][1],
                                        label, modification_type)

                    space_above_sequence = utils.get_height()-y_0_line if above == 'A' else y_0_line
                    space_per_group = (space_above_sequence - label_plot_height)  // (len(df["Group"].unique())-1)
                    max_bar_height = space_per_group - 2*bar_plot_margin
                    # plot bars
                    for i, group in enumerate(self.PLOT_CONFIG.BAR_GROUPS.keys()):
                        modification_column = [col for col in df.columns if col not in ['ID', 'Group'] and df[col].iloc[1] == label]
                        bar_df = df[['ID', 'Group'] + modification_column]
                        bar_df = bar_df[bar_df['Group'] == group]
                        values = bar_df.iloc[:, 2].astype(int)
                        percentage = values.mean()
                        height = max_bar_height * percentage
                        bar_percentages[group].append(percentage)
                        x0 = x_1_line - bar_width // 2 * self.PLOT_CONFIG.BAR_WIDTH
                        x1 = x_1_line + bar_width // 2 * self.PLOT_CONFIG.BAR_WIDTH
                        y0 = y_bar + (i * space_per_group + 5) * group_direction
                        y1 = y0 + height * group_direction
                        if self.PLOT_CONFIG.INVERT_AXIS_GROUP_B and above == 'B':
                            y0 = y_bar + (i * space_per_group + 5 + max_bar_height) * group_direction
                            y1 = y0 - height * group_direction

                        fig.add_shape(
                            type="rect",
                            x0=x0,
                            y0=y0,
                            x1=x1,
                            y1=y1,
                            line={'color': "black", 'width': 1},
                            fillcolor=self.CONFIG.MODIFICATIONS[modification_type][1]
                        )
                else:
                    position = utils.get_position_with_offset(aa_position, isoform)
                    y_0_line = utils.get_height() - (position * utils.PIXELS_PER_AA + utils.SEQUENCE_OFFSET)
                    y_0_line = utils.offset_line_for_exon(y_0_line, aa_position, self.CONFIG.FIGURE_ORIENTATION)
                    y_1_line = modifications_visited * bar_width + bar_width//2
                    x_0_line = utils.SEQUENCE_BOUNDARIES['x1'] if above == 'A' else utils.SEQUENCE_BOUNDARIES['x0']
                    x_1_line = x_0_line + group_direction * height_offset
                    x_3_line = x_0_line + label_plot_height * group_direction - (utils.get_label_length(label) + 10) * group_direction
                    x_2_line = x_3_line - 10 * group_direction
                    x_label = x_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction
                    x_bar = x_3_line + (utils.get_label_length(label) + 5) * group_direction

                    # plot line with label
                    self.plot_line_with_label_vertical(fig,
                                                x_0_line, x_1_line, x_2_line, x_3_line, x_label,
                                                y_0_line, y_1_line,
                                                self.CONFIG.MODIFICATIONS[modification_type][1],
                                                label, modification_type)

                    space_above_sequence = utils.get_width()-x_0_line if above == 'A' else x_0_line
                    space_per_group = (space_above_sequence - label_plot_height)  // (len(df["Group"].unique())-1)
                    max_bar_height = space_per_group - 2*bar_plot_margin
                    # plot bars
                    for i, group in enumerate(self.PLOT_CONFIG.BAR_GROUPS.keys()):
                        modification_column = [col for col in df.columns if col not in ['ID', 'Group'] and df[col].iloc[1] == label]
                        bar_df = df[['ID', 'Group'] + modification_column]
                        bar_df = bar_df[bar_df['Group'] == group]
                        values = bar_df.iloc[:, 2].astype(int)
                        percentage = values.mean()
                        height = max_bar_height * percentage
                        bar_percentages[group].append(percentage)
                        y0 = y_1_line - bar_width // 2 * self.PLOT_CONFIG.BAR_WIDTH
                        y1 = y_1_line + bar_width // 2 * self.PLOT_CONFIG.BAR_WIDTH
                        x0 = x_bar + (i * space_per_group + 5) * group_direction
                        x1 = x0 + height * group_direction
                        if self.PLOT_CONFIG.INVERT_AXIS_GROUP_B and above == 'B':
                            x0 = x_bar + (i * space_per_group + 5 + max_bar_height) * group_direction
                            x1 = x0 - height * group_direction

                        fig.add_shape(
                            type="rect",
                            x0=x0,
                            y0=y0,
                            x1=x1,
                            y1=y1,
                            line={'color': "black", 'width': 1},
                            fillcolor=self.CONFIG.MODIFICATIONS[modification_type][1]
                        )

                modifications_visited += 1
            positions_visited += 1
            if positions_visited >= len(modification_sights_relevant)//2:
                height_offset -= 2
            else:
                height_offset += 2

        max_lines = max(value.count('<br>')+1 for value in self.PLOT_CONFIG.BAR_GROUPS.values())
        if self.CONFIG.FIGURE_ORIENTATION == 0:
            for i, group in enumerate(self.PLOT_CONFIG.BAR_GROUPS.keys()):
                y_group = y_bar + (i * space_per_group + 5) * group_direction
                x_line_start = utils.get_width()-bar_plot_width
                fig.add_annotation(x=x_line_start-utils.get_label_length('100%')-max_lines*utils.get_label_height()+3,
                                    y = y_group + space_per_group//2 * group_direction,
                                    width=space_per_group,
                                    text=self.PLOT_CONFIG.BAR_GROUPS[group],
                                    textangle=-90,
                                    font={'size': self.CONFIG.SEQUENCE_PLOT_FONT_SIZE, 'family': self.CONFIG.FONT, 'color':"black"},
                                    showarrow=False)
                for j in range(5):
                    if j == 0:
                        y_trace = y_group
                    elif j == 4:
                        y_trace = y_group + max_bar_height * group_direction
                    else:
                        y_trace = y_group + round(max_bar_height/4, 1) * j * group_direction
                    fig.add_trace(go.Scatter(x=[x_line_start, x_line_start+bar_plot_width],
                                            y=[y_trace, y_trace],
                                            mode='lines',
                                            line={'color':'lightgray', 'width':1},
                                            showlegend=False,
                                            hoverinfo='none'))
                if j % 2 == 0:
                    text = str(j*25) + '%'
                    if self.PLOT_CONFIG.INVERT_AXIS_GROUP_B and above == 'B':
                        text = str(100-j*25) + '%'
                    fig.add_annotation(x=x_line_start-utils.get_label_length(text)//2-5,
                                        y=y_trace,
                                        text=text,
                                        font={'size': self.CONFIG.SEQUENCE_PLOT_FONT_SIZE, 'family': self.CONFIG.FONT, 'color':"black"},
                                        showarrow=False)
        else:
            for i, group in enumerate(self.PLOT_CONFIG.BAR_GROUPS.keys()):
                x_group = x_bar + (i * space_per_group + 5) * group_direction
                y_line_start = bar_plot_width
                fig.add_annotation(x=x_group + space_per_group//2 * group_direction,
                                y=y_line_start+utils.get_label_height()+max_lines*utils.get_label_height()-5,
                                text=self.PLOT_CONFIG.BAR_GROUPS[group],
                                font={'size': self.CONFIG.SEQUENCE_PLOT_FONT_SIZE, 'family': self.CONFIG.FONT, 'color':"black"},
                                showarrow=False)
                for j in range(5):
                    if j == 0:
                        x_trace = x_group
                    elif j == 4:
                        x_trace = x_group + max_bar_height * group_direction
                    else:
                        x_trace = x_group + round(max_bar_height/4, 1) * j * group_direction
                    fig.add_trace(go.Scatter(x=[x_trace, x_trace],
                                        y=[y_line_start, 0],
                                        mode='lines',
                                        line={'color':'lightgray', 'width':1},
                                        showlegend=False,
                                        hoverinfo='none'))
                    if j % 2 == 0:
                        text = str(j*25) + '%'
                        if self.PLOT_CONFIG.INVERT_AXIS_GROUP_B and above == 'B':
                            text = str(100-j*25) + '%'
                        if j == 0:
                            x_pos = x_trace - utils.get_label_length(text)//2 if above == 'B' else x_trace + utils.get_label_length(text)//2
                        elif j == 4:
                            x_pos = x_trace + utils.get_label_length(text)//2 if above == 'B' else x_trace - utils.get_label_length(text)//2
                        else:
                            x_pos = x_trace
                        fig.add_annotation(x=x_pos,
                                                y=y_line_start+5,
                                                text=text,
                                                font={'size': self.CONFIG.SEQUENCE_PLOT_FONT_SIZE, 'family': self.CONFIG.FONT, 'color':"black"},
                                                showarrow=False)

        return fig

    def plot_line_with_label_horizontal(self, fig, x_0, x_1, y_0, y_1, y_2, y_3, y_label, color, label, modification_type):
        """Plot single line with label in horizontal orientation."""
        fig.add_trace(go.Scatter(x=[x_0, x_0, x_1, x_1],
                                    y=[y_0, y_1, y_2, y_3],
                                    mode='lines',
                                    line={'color':color, 'width':1}, showlegend=False, hoverinfo='none'))
        fig.add_annotation(x=x_1, y=y_label,
                            text=label,
                            showarrow=False,
                            textangle=-90,
                            font={'family': self.CONFIG.FONT, 'size': self.CONFIG.SEQUENCE_PLOT_FONT_SIZE, 'color':color})
        if f'{modification_type}({label[0]})@{label[1:]}' in self.CONFIG.PTMS_TO_HIGHLIGHT:
            fig.add_shape(type='rect',
                                x0 = x_1-utils.get_label_height()//2-1,
                                x1 = x_1+utils.get_label_height()//2+1,
                                y0 = y_label-utils.get_label_length(label)//2-3,
                                y1 = y_label+utils.get_label_length(label)//2+3,
                                line={'width': 0},
                                fillcolor=self.CONFIG.PTM_HIGHLIGHT_LABEL_COLOR,
                                showlegend=False)
        return fig

    def plot_line_with_label_vertical(self, fig, x_0, x_1, x_2, x_3, x_label, y_0, y_1, color, label, modification_type):
        """Plot single line with label in vertical orientation."""
        fig.add_trace(go.Scatter(x=[x_0, x_1, x_2, x_3],
                                y=[y_0, y_0, y_1, y_1],
                                mode='lines',
                                line={'color': color, 'width': 1}, showlegend=False, hoverinfo='none'))
        fig.add_annotation(x=x_label, y=y_1,
                        text=label,
                        showarrow=False,
                        font={'family': self.CONFIG.FONT, 'size': self.CONFIG.SEQUENCE_PLOT_FONT_SIZE, 'color': color})
        if f'{modification_type}({label[0]})@{label[1:]}' in self.CONFIG.PTMS_TO_HIGHLIGHT:
            fig.add_shape(type='rect',
                                x0 = x_label-utils.get_label_length(label)//2-3,
                                x1 = x_label+utils.get_label_length(label)//2+3,
                                y0 = y_1-utils.get_label_height()//2-1,
                                y1 = y_1+utils.get_label_height()//2+1,
                                line={'width': 0},
                                fillcolor=self.CONFIG.PTM_HIGHLIGHT_LABEL_COLOR,
                                showlegend=False)
        return fig

    def filter_relevant_modification_sights(self, helper_file: str):
        """Filter relevant modification sights from input file based on user defined filters."""
        # read csv into df
        df = pd.read_csv(helper_file)

        # only keep first two columns and columns that are in MODIFICATIONS
        columns_to_keep = list(self.CONFIG.INCLUDED_MODIFICATIONS.keys())
        df = df[[col for col in df.columns if df[col][0] in columns_to_keep or col in ['ID', 'Group']]]

        # create dict for all modification sights
        all_modification_sights = defaultdict(list)
        for column in df.columns:
            if column in ['ID', 'Group']:
                continue
            column_modification = df[column][0]
            column_label = df[column][1][0]
            if column_modification not in self.PLOT_CONFIG.MODIFICATIONS_GROUP:
                continue
            if self.CONFIG.INCLUDED_MODIFICATIONS.get(column_modification):
                if column_label not in self.CONFIG.INCLUDED_MODIFICATIONS[column_modification]:
                    continue
                if column_label == 'R' and column_modification == 'Deamidated':
                    column_modification = 'Citrullination'
            isoform = df[column][2]
            all_modification_sights[int(df[column][1][1:])].append((df[column][1], column_modification, self.PLOT_CONFIG.MODIFICATIONS_GROUP[column_modification], isoform))

        # remove rows where groups is not in self.PLOT_CONFIG.BAR_GROUPS
        header_rows = df.iloc[:4, :]
        df = df[df['Group'].isin(self.PLOT_CONFIG.BAR_GROUPS.keys())]
        df = pd.concat([header_rows, df], ignore_index=True)

        # filter out columns that have no modifications observed for relevant groups
        data_to_filter = df.iloc[4:, 2:]
        header_columns = df.iloc[:, :2]
        for col in data_to_filter.columns:
            data_to_filter[col] = data_to_filter[col].astype(int)
        df = df[[col for col in data_to_filter if data_to_filter[col].sum(axis=0) > 0]]
        df = pd.concat([header_columns, df], axis=1)

        # create new dict for modification sights
        relevant_modification_sights = defaultdict(list)
        for column in df.columns:
            if column in ['ID', 'Group']:
                continue
            column_modification = df[column][0]
            column_label = df[column][1][0]
            if column_modification not in self.PLOT_CONFIG.MODIFICATIONS_GROUP:
                continue
            if self.CONFIG.INCLUDED_MODIFICATIONS.get(column_modification):
                if column_label not in self.CONFIG.INCLUDED_MODIFICATIONS[column_modification]:
                    continue
                if column_label == 'R' and column_modification == 'Deamidated':
                    column_modification = 'Citrullination'
            isoform = df[column][2]
            relevant_modification_sights[int(df[column][1][1:])].append((df[column][1], column_modification, self.PLOT_CONFIG.MODIFICATIONS_GROUP[column_modification], isoform))
        return all_modification_sights, relevant_modification_sights, df

    def get_relevant_mod_types(self, relevant_positions: dict[int, list[tuple[int, str, str, str]]]) -> set[str]:
        """Get relevant modification types."""
        present_mod_types = set()
        for aa_position in relevant_positions.keys():
            for modification_sight in relevant_positions[aa_position]:
                present_mod_types.add(modification_sight[1])
        return present_mod_types

    def create_bar_plot(self):
        """Main function to create bar plot."""
        all_positions, relevant_positions, df = self.filter_relevant_modification_sights(self.PLOT_CONFIG.BAR_INPUT_FILE)
        above_all, below_all = utils.separate_by_group(all_positions)
        above_relevant, below_relevant = utils.separate_by_group(relevant_positions)
        present_mod_types = self.get_relevant_mod_types(relevant_positions)


        if len(above_relevant) == 0:
            fig = sequence_plot.create_plot(self.input_file, present_mod_types, 'A', 'A')
        elif len(below_relevant) == 0:
            fig = sequence_plot.create_plot(self.input_file, present_mod_types, 'B', 'B')
        else:
            legend = 'A' if self.CONFIG.FIGURE_ORIENTATION == 0 else 'B'
            fig = sequence_plot.create_plot(self.input_file, present_mod_types, None, legend)

        positions_a, positions_b = self.get_bar_positions(above_all, below_all)
        group_size_a = len(positions_a)
        group_size_b = len(positions_b)
        bar_plot_width = self.get_bar_plot_width(group_size_a, group_size_b)
        highest_position = positions_a[-1] if group_size_a > group_size_b else positions_b[-1]
        label_plot_height = max(group_size_a, group_size_b) + utils.get_label_length(f'X{highest_position}') + 30

        for (group_label, group_all, group_relevant, group_positions) in [('A', above_all, above_relevant, positions_a), ('B', below_all, below_relevant, positions_b)]:
            if len(group_positions) == 0:
                continue
            fig = self.add_bar_plot(fig, group_label, group_all, group_relevant, df, group_positions, bar_plot_width, label_plot_height)

        utils.show_plot(fig, self.output_path)
