"""Module for plotting cleavages and PTMs on the sequence plot."""

from pathlib import Path

import logging
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from protein_sequencing.plotter import Plotter


class DetailsPlotter(Plotter):
    """Class to plot cleavages and PTMs on the sequence plot."""

    def __init__(self, config, plot_config, input_file, output_path):
        super().__init__(config, plot_config)

        self.input_file = input_file
        self.output_path = output_path
        if not Path(self.output_path).exists():
            Path(self.output_path).mkdir(parents=True, exist_ok=True)

    def get_present_regions(self, positions, isoforms):
        ranges = []
        for position_range in positions:
            if "-" in str(position_range):
                start, end = map(int, position_range.split("-"))
                ranges.append((start, end))
            else:
                start = end = int(position_range)
                ranges.append((start, end))

        region_ranges = []
        region_start = 1
        for _, region_end, _, _ in self.REGIONS:
            region_ranges.append((region_start, region_end))
            region_start = region_end + 1

        regions_present = [False] * len(region_ranges)
        region_index = 0
        for i, position_range in enumerate(ranges):
            if isoforms[i] == "exon1":
                index = next(
                    (
                        index
                        for index, region in enumerate(self.REGIONS)
                        if region[1] == self.EXON_1_OFFSET["index_end"]
                    ),
                    None,
                )
                if index:
                    regions_present[index] = True
            elif isoforms[i] == "exon2":
                index = next(
                    (
                        index
                        for index, region in enumerate(self.REGIONS)
                        if region[1] == self.EXON_2_OFFSET["index_end"]
                    ),
                    None,
                )
                if index:
                    regions_present[index] = True
            # TODO: add check here that this is not overflowing
            while position_range[0] > region_ranges[region_index][1]:
                region_index += 1
            regions_present[region_index] = True
        return regions_present

    def get_present_regions_cleavage(self, cleavage_df: pd.DataFrame):
        cleavages = cleavage_df.iloc[1:2, 2:].values[0].tolist()
        isoforms = cleavage_df.iloc[2:3, 2:].values[0].tolist()
        return self.get_present_regions(cleavages, isoforms)

    def get_present_regions_ptm(self, ptm_df: pd.DataFrame):
        ptms = ptm_df.iloc[1:2, 2:].values[0].tolist()
        ptms = [ptm[1:] for ptm in ptms]
        isoforms = ptm_df.iloc[2:3, 2:].values[0].tolist()
        return self.get_present_regions(ptms, isoforms)

    def plot_line_with_label_horizontal(
        self,
        fig: go.Figure,
        x_0: int,
        x_1: int,
        y_0: int,
        y_1: int,
        y_2: int,
        y_3: int,
        y_label: int,
        label: str,
        ptm: bool,
        ptm_color: str | None = None,
        ptm_modification: str | None = None,
    ):
        line_color = "black"
        if ptm:
            line_color = ptm_color
        fig.add_trace(
            go.Scatter(
                x=[x_0, x_0, x_1, x_1],
                y=[y_0, y_1, y_2, y_3],
                mode="lines",
                line={"color": line_color, "width": 1},
                showlegend=False,
                hoverinfo="none",
            )
        )
        if ptm:
            color = ptm_color
            if f"{ptm_modification}({label[0]})@{label[1:]}" in self.PTMS_TO_HIGHLIGHT:
                fig.add_shape(
                    type="rect",
                    x0=x_1 - self.get_label_height() // 2 - 1,
                    x1=x_1 + self.get_label_height() // 2 + 1,
                    y0=y_label - self.get_label_length(label) // 2 - 3,
                    y1=y_label + self.get_label_length(label) // 2 + 3,
                    line={"width": 0},
                    fillcolor=self.PTM_HIGHLIGHT_LABEL_COLOR,
                    showlegend=False,
                )
        else:
            color = self.plot_config.CLEAVAGE_LABEL_COLOR
            if label in self.plot_config.CLEAVAGES_TO_HIGHLIGHT:
                color = self.plot_config.CLEAVAGE_HIGHLIGHT_COLOR
        fig.add_annotation(
            x=x_1,
            y=y_label,
            text=label,
            showarrow=False,
            textangle=-90,
            font=dict(
                family=self.FONT,
                size=self.sequence_plot_font_size,
                color=color,
            ),
        )
        return fig

    def plot_line_with_label_vertical(
        self,
        fig: go.Figure,
        x_0: int,
        x_1: int,
        x_2: int,
        x_3: int,
        y_0: int,
        y_1: int,
        x_label: int,
        label: str,
        ptm: bool,
        ptm_color: str | None = None,
        ptm_modification: str | None = None,
    ):
        line_color = "black"
        if ptm:
            line_color = ptm_color
        fig.add_trace(
            go.Scatter(
                x=[x_0, x_1, x_2, x_3],
                y=[y_0, y_0, y_1, y_1],
                mode="lines",
                line={"color": line_color, "width": 1},
                showlegend=False,
                hoverinfo="none",
            )
        )
        if ptm:
            color = ptm_color
            if f"{ptm_modification}({label[0]})@{label[1:]}" in self.PTMS_TO_HIGHLIGHT:
                fig.add_shape(
                    type="rect",
                    x0=x_label - self.get_label_length(label) // 2 - 3,
                    x1=x_label + self.get_label_length(label) // 2 + 3,
                    y0=y_1 - self.get_label_height() // 2 - 1,
                    y1=y_1 + self.get_label_height() // 2 + 1,
                    line={"width": 0},
                    fillcolor=self.PTM_HIGHLIGHT_LABEL_COLOR,
                    showlegend=False,
                )
        else:
            color = self.plot_config.CLEAVAGE_LABEL_COLOR
            if label in self.plot_config.CLEAVAGES_TO_HIGHLIGHT:
                color = self.plot_config.CLEAVAGE_HIGHLIGHT_COLOR
        fig.add_annotation(
            x=x_label,
            y=y_1,
            text=label,
            showarrow=False,
            font={
                "family": self.FONT,
                "size": self.sequence_plot_font_size,
                "color": color,
            },
        )
        return fig

    def plot_range_with_label_horizontal(
        self,
        fig: go.Figure,
        x_0_start: int,
        x_0_end: int,
        x_1: int,
        y_0: int,
        y_1: int,
        y_2: int,
        y_3: int,
        y_label: int,
        label: str,
    ):
        fig.add_trace(
            go.Scatter(
                x=[x_0_start, x_0_start, x_1, x_1, x_1, x_0_end, x_0_end],
                y=[y_0, y_1, y_2, y_3, y_2, y_1, y_0],
                mode="lines",
                fill="toself",
                line={"color": "black", "width": 1},
                showlegend=False,
                hoverinfo="none",
            )
        )

        color = self.plot_config.CLEAVAGE_LABEL_COLOR
        if label in self.plot_config.CLEAVAGES_TO_HIGHLIGHT:
            color = self.plot_config.CLEAVAGE_HIGHLIGHT_COLOR
        fig.add_annotation(
            x=x_1,
            y=y_label,
            text=label,
            showarrow=False,
            textangle=-90,
            font={
                "family": self.FONT,
                "size": self.sequence_plot_font_size,
                "color": color,
            },
        )
        return fig

    def plot_range_with_label_vertical(
        self,
        fig: go.Figure,
        x_0: int,
        x_1: int,
        x_2: int,
        x_3: int,
        y_0_start: int,
        y_0_end: int,
        y_1: int,
        x_label: int,
        label: str,
    ):
        fig.add_trace(
            go.Scatter(
                x=[x_0, x_1, x_2, x_3, x_2, x_1, x_0],
                y=[y_0_start, y_0_start, y_1, y_1, y_1, y_0_end, y_0_end],
                mode="lines",
                fill="toself",
                line={"color": "black", "width": 1},
                showlegend=False,
                hoverinfo="none",
            )
        )
        color = self.plot_config.CLEAVAGE_LABEL_COLOR
        if label in self.plot_config.CLEAVAGES_TO_HIGHLIGHT:
            color = self.plot_config.CLEAVAGE_HIGHLIGHT_COLOR
        fig.add_annotation(
            x=x_label,
            y=y_1,
            text=label,
            showarrow=False,
            font={
                "family": self.FONT,
                "size": self.sequence_plot_font_size,
                "color": color,
            },
        )
        return fig

    def plot_groups_horizontal(
        self,
        fig: go.Figure,
        df: pd.DataFrame,
        x_0_groups: int,
        y_0_groups: int,
        dx: int,
        dy: int,
        x_label: int,
        y_label: int,
        last_region: int,
        group_dircetion: int,
        ptm: bool,
    ):
        x_margin = 0
        if dx % 2 != 0:
            x_margin = 1
        color_low = self.plot_config.CLEAVAGE_SCALE_COLOR_LOW
        color_mid = self.plot_config.CLEAVAGE_SCALE_COLOR_MID
        color_high = self.plot_config.CLEAVAGE_SCALE_COLOR_HIGH
        if ptm:
            color_low = self.plot_config.PTM_SCALE_COLOR_LOW
            color_mid = self.plot_config.PTM_SCALE_COLOR_MID
            color_high = self.plot_config.PTM_SCALE_COLOR_HIGH
        fig.add_shape(
            type="rect",
            x0=x_0_groups - dx // 2 - x_margin,
            y0=y_0_groups,
            x1=x_0_groups + dx * len(df.iloc[0:1, :].columns) - dx // 2,
            y1=y_0_groups + dy * len(df.index) + 1,
            fillcolor="grey",
            line={"color": "grey", "width": 1},
            showlegend=False,
            layer="below",
        )
        df.columns = np.arange(len(df.columns))  # Otherwise Heatmap will complain
        fig.add_trace(
            go.Heatmap(
                z=df,
                x0=x_0_groups,
                y0=y_0_groups + dy // 2,
                dx=dx,
                dy=dy,
                showscale=False,
                hoverinfo="none",
                xgap=1,
                ygap=1,
                zmin=0,
                zmax=1,
                zmid=0.5,
                colorscale=[[0, color_low], [0.5, color_mid], [1, color_high]],
            )
        )
        yanchor = "bottom"
        xanchor = "left"
        if group_dircetion == -1:
            yanchor = "top"
            xanchor = "right"
        fig.add_annotation(
            x=x_label - self.get_label_height() * group_dircetion,
            y=y_label
            - int((8 / self.offset_region_label_from_angle()) * group_dircetion),
            text=self.REGIONS[last_region][3],
            showarrow=False,
            textangle=-self.plot_config.REGION_LABEL_ANGLE_GROUPS,
            xanchor=xanchor,
            yanchor=yanchor,
            font={
                "family": self.FONT,
                "size": self.sequence_plot_font_size,
                "color": "black",
            },
        )
        return fig

    def plot_groups_vertical(
        self,
        fig: go.Figure,
        df: pd.DataFrame,
        x_0_groups: int,
        y_0_groups: int,
        dx: int,
        dy: int,
        x_label: int,
        y_label: int,
        last_region: int,
        group_dircetion: int,
        ptm: bool,
    ):
        y_margin = 0
        if dy % 2 != 0:
            y_margin = 1
        color_low = self.plot_config.CLEAVAGE_SCALE_COLOR_LOW
        color_mid = self.plot_config.CLEAVAGE_SCALE_COLOR_MID
        color_high = self.plot_config.CLEAVAGE_SCALE_COLOR_HIGH
        if ptm:
            color_low = self.plot_config.PTM_SCALE_COLOR_LOW
            color_mid = self.plot_config.PTM_SCALE_COLOR_MID
            color_high = self.plot_config.PTM_SCALE_COLOR_HIGH
        fig.add_shape(
            type="rect",
            x0=x_0_groups,
            y0=y_0_groups + dy // 2 + y_margin,
            x1=x_0_groups + dx * len(df.index) + 1,
            y1=y_0_groups - dy * len(df.iloc[0:1, :].columns) + dy // 2,
            fillcolor="grey",
            line={"color": "grey", "width": 1},
            showlegend=False,
            layer="below",
        )

        fig.add_trace(
            go.Heatmap(
                z=df.T,
                x0=x_0_groups + dx // 2,
                y0=y_0_groups,
                dx=dx,
                dy=-dy,
                showscale=False,
                hoverinfo="none",
                xgap=1,
                ygap=1,
                zmin=0,
                zmax=1,
                zmid=0.5,
                colorscale=[[0, color_low], [0.5, color_mid], [1, color_high]],
            )
        )
        xanchor = "left"
        if group_dircetion == -1:
            xanchor = "right"
        fig.add_annotation(
            x=x_label,
            y=y_label,
            text=self.REGIONS[last_region][3],
            showarrow=False,
            textangle=-self.plot_config.REGION_LABEL_ANGLE_GROUPS + 90,
            xanchor=xanchor,
            font={
                "family": self.FONT,
                "size": self.sequence_plot_font_size,
                "color": "black",
            },
        )
        return fig

    def plot_group_labels_horizontal(
        self, fig: go.Figure, mean_values: pd.DataFrame, y_0_groups: int, dy: int
    ):
        for i, group in enumerate(mean_values.index):
            y_0_rect = y_0_groups + i * dy
            x_1_rect = self.calculate_group_space()
            fig.add_shape(
                type="rect",
                x0=0,
                x1=self.calculate_group_space(),
                y0=y_0_rect,
                y1=y_0_rect + dy,
                fillcolor=self.plot_config.GROUPS[group][1],
                line={"width": 0},
                showlegend=False,
                layer="below",
            )
            color = self.get_label_color(group)

            fig.add_annotation(
                x=x_1_rect // 2,
                y=y_0_rect + dy // 2,
                text=group,
                showarrow=False,
                align="center",
                font={
                    "family": self.FONT,
                    "size": self.sequence_plot_font_size,
                    "color": color,
                },
            )

    def get_label_color(self, group: str):
        # based on https://stackoverflow.com/questions/3942878/
        red, green, blue = tuple(
            int(self.plot_config.GROUPS[group][1][i : i + 2], 16) for i in (1, 3, 5)
        )
        return (
            "#000000" if red * 0.299 + green * 0.587 + blue * 0.114 > 130 else "#ffffff"
        )

    def plot_group_labels_vertical(
        self, fig: go.Figure, mean_values: pd.DataFrame, x_0_groups: int, dx: int
    ):
        for i, group in enumerate(mean_values.index):
            x_0_rect = x_0_groups + i * dx
            y_0_rect = self.get_height()
            y_rect = self.calculate_group_space()
            fig.add_shape(
                type="rect",
                x0=x_0_rect,
                x1=x_0_rect + dx,
                y0=y_0_rect,
                y1=y_0_rect - y_rect,
                fillcolor=self.plot_config.GROUPS[group][1],
                line=dict(width=0),
                showlegend=False,
                layer="below",
            )

            color = self.get_label_color(group)

            fig.add_annotation(
                x=x_0_rect + dx // 2,
                y=y_0_rect - y_rect // 2,
                text=group,
                showarrow=False,
                align="center",
                textangle=90,
                font=dict(
                    family=self.FONT, size=self.sequence_plot_font_size, color=color
                ),
            )

    def preprocess_groups(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[tuple[str, str]]]:
        df.columns = df.iloc[0]
        labels = df.iloc[1:2, 2:].values.flatten().tolist()
        regions = df.iloc[2:3, 2:].values.flatten().tolist()
        label_region_type_map = list(zip(labels, regions))
        df = df.iloc[3:]

        reverse_group_mapping = {}
        for k, v in self.plot_config.GROUPS.items():
            values = v[0]
            for value in values:
                reverse_group_mapping[value] = k
        df.iloc[:, 1] = df.iloc[:, 1].map(reverse_group_mapping)
        mean_values = df.iloc[:, 2:].astype(float).groupby(df.iloc[:, 1]).mean()
        mean_values = mean_values.reindex([*self.plot_config.GROUPS])

        return mean_values, label_region_type_map

    def offset_region_label_from_angle(self):
        longest_label = ""
        for _, _, _, region_label_short in self.REGIONS:
            if self.get_label_length(region_label_short) > self.get_label_length(
                longest_label
            ):
                longest_label = region_label_short

        length = self.get_label_length(longest_label)
        height = self.get_label_height()

        angle_radians = math.radians(-self.plot_config.REGION_LABEL_ANGLE_GROUPS)
        dy = abs((length / 2) * math.sin(angle_radians)) + abs(
            (height / 2) * math.cos(angle_radians)
        )

        return int(dy) + 10

    def _calculate_cleavage_line_coordinates_horizontal(
        self,
        above: str,
        label_plot_height: int,
        longest_label: str,
        group_direction: int,
    ):
        y_0_line = (
            self.SEQUENCE_BOUNDARIES["y1"]
            if above == "A"
            else self.SEQUENCE_BOUNDARIES["y0"]
        )
        y_1_line = y_0_line + 10 * group_direction
        y_2_line = (
            y_0_line
            + (label_plot_height - self.get_label_length(longest_label) - 10)
            * group_direction
        )
        return y_0_line, y_1_line, y_2_line

    def _calculate_cleavage_line_coordinates_vertical(
        self,
        above: str,
        label_plot_height: int,
        longest_label: str,
        group_direction: int,
    ):
        x_0_line = (
            self.SEQUENCE_BOUNDARIES["x1"]
            if above == "A"
            else self.SEQUENCE_BOUNDARIES["x0"]
        )
        x_1_line = x_0_line + 10 * group_direction
        x_2_line = (
            x_0_line
            + (label_plot_height - self.get_label_length(longest_label) - 10)
            * group_direction
        )
        return x_0_line, x_1_line, x_2_line

    def _setup_cleavage_group_spacing_horizontal(
        self,
        pixels_per_cleavage: int,
        y_0_line: int,
        label_plot_height: int,
        above: str,
        group_direction: int,
        mean_values: pd.DataFrame,
    ):
        y_0_groups = y_0_line + (label_plot_height + 10) * group_direction
        vertical_space_left = (
            self.get_height() - y_0_groups if above == "A" else y_0_groups
        )
        # offset for border around heatmap
        vertical_space_left -= 2
        # offset for label for region
        dy_label = self.offset_region_label_from_angle()
        vertical_space_left -= dy_label * 2
        dx = pixels_per_cleavage
        dy = vertical_space_left // len(mean_values.index) * group_direction
        return dx, dy, y_0_groups, vertical_space_left

    def _setup_cleavage_group_spacing_vertical(
        self,
        pixels_per_cleavage: int,
        x_0_line: int,
        label_plot_height: int,
        above: str,
        group_direction: int,
        mean_values: pd.DataFrame,
    ):
        x_0_groups = x_0_line + (label_plot_height + 10) * group_direction
        horizontal_space_left = (
            self.get_width() - x_0_groups if above == "A" else x_0_groups
        )
        # offset for border around heatmap
        horizontal_space_left -= 2
        # offset for label for region
        dx_label = self.offset_region_label_from_angle()
        horizontal_space_left -= dx_label * 2
        dy = pixels_per_cleavage
        dx = horizontal_space_left // len(mean_values.index) * group_direction
        return dx, dy, x_0_groups, horizontal_space_left

    def _plot_region_divider_horizontal(
        self,
        fig: go.Figure,
        item_idx: int,
        first_item_in_region: int,
        i: int,
        pixels_per_item: int,
        dx: int,
        dy: int,
        mean_values: pd.DataFrame,
        y_0_groups: int,
        last_region: int,
        group_direction: int,
        is_ptm: bool,
    ):
        start_idx = item_idx - (i - first_item_in_region)
        x_0_groups = start_idx * pixels_per_item + self.get_horizontal_offset(dx)
        x_divider = item_idx * pixels_per_item + self.get_horizontal_offset(dx)
        x_label = x_0_groups + (x_divider - x_0_groups) // 2 - dx // 2
        y_label = (
            y_0_groups
            + len(mean_values.index) * dy
            + (5 + self.get_label_height() // 2) * group_direction
        )

        self.plot_groups_horizontal(
            fig,
            mean_values.iloc[:, first_item_in_region:i],
            x_0_groups,
            y_0_groups,
            dx,
            dy,
            x_label,
            y_label,
            last_region,
            group_direction,
            is_ptm,
        )

        fig.add_trace(
            go.Scatter(
                x=[x_divider, x_divider],
                y=[y_0_groups, y_0_groups + len(mean_values.index) * dy],
                mode="lines",
                line=dict(color="black", width=3),
                showlegend=False,
                hoverinfo="none",
            )
        )

    def _plot_region_divider_vertical(
        self,
        fig: go.Figure,
        item_idx: int,
        first_item_in_region: int,
        i: int,
        pixels_per_item: int,
        dx: int,
        dy: int,
        mean_values: pd.DataFrame,
        x_0_groups: int,
        last_region: int,
        group_direction: int,
        is_ptm: bool,
    ):
        start_idx = item_idx - (i - first_item_in_region)
        y_0_groups = (
            self.get_height()
            - start_idx * pixels_per_item
            - self.get_vertical_offset(dy)
        )
        y_divider = (
            self.get_height()
            - item_idx * pixels_per_item
            - self.get_vertical_offset(dy)
        )
        y_label = y_0_groups - (y_0_groups - y_divider) // 2 + dy // 2
        x_label = (
            x_0_groups
            + len(mean_values.index) * dx
            + (5 + self.get_label_height() // 2) * group_direction
        )

        self.plot_groups_vertical(
            fig,
            mean_values.iloc[:, first_item_in_region:i],
            x_0_groups,
            y_0_groups,
            dx,
            dy,
            x_label,
            y_label,
            last_region,
            group_direction,
            is_ptm,
        )

        fig.add_trace(
            go.Scatter(
                x=[x_0_groups, x_0_groups + len(mean_values.index) * dx],
                y=[y_divider, y_divider],
                mode="lines",
                line=dict(color="black", width=3),
                showlegend=False,
                hoverinfo="none",
            )
        )

    def _plot_single_cleavage_horizontal(
        self,
        fig: go.Figure,
        start: int,
        end: int,
        isoform: str,
        cleavage_idx: int,
        pixels_per_cleavage: int,
        dx: int,
        y_0_line: int,
        y_1_line: int,
        y_2_line: int,
        label_plot_height: int,
        group_direction: int,
    ):
        if start == end:
            label = str(start)
            position = self.get_position_with_offset(start, isoform)
            x_0_line = position * self.PIXELS_PER_AA + self.SEQUENCE_OFFSET
            x_0_line = self.offset_line_for_exon(
                x_0_line, start, self.FIGURE_ORIENTATION, isoform
            )
            x_1_line = cleavage_idx * pixels_per_cleavage + self.get_horizontal_offset(
                dx
            )
            y_3_line = (
                y_0_line
                + (label_plot_height - self.get_label_length(label)) * group_direction
            )
            y_label = (
                y_3_line + (self.get_label_length(label) // 2 + 5) * group_direction
            )

            self.plot_line_with_label_horizontal(
                fig,
                x_0_line,
                x_1_line,
                y_0_line,
                y_1_line,
                y_2_line,
                y_3_line,
                y_label,
                label,
                False,
                None,
                None,
            )
        else:
            label = f"{start}-{end}"
            start_position = self.get_position_with_offset(start, isoform)
            end_position = self.get_position_with_offset(end, isoform)
            x_0_start_line = start_position * self.PIXELS_PER_AA + self.SEQUENCE_OFFSET
            x_0_end_line = end_position * self.PIXELS_PER_AA + self.SEQUENCE_OFFSET
            x_0_start_line = self.offset_line_for_exon(
                x_0_start_line, start, self.FIGURE_ORIENTATION, isoform
            )
            x_0_end_line = self.offset_line_for_exon(
                x_0_end_line, end, self.FIGURE_ORIENTATION, isoform
            )
            x_1_line = cleavage_idx * pixels_per_cleavage + self.get_horizontal_offset(
                dx
            )
            y_3_line = (
                y_0_line
                + (label_plot_height - self.get_label_length(label)) * group_direction
            )
            y_label = (
                y_3_line + (self.get_label_length(label) // 2 + 5) * group_direction
            )

            self.plot_range_with_label_horizontal(
                fig,
                x_0_start_line,
                x_0_end_line,
                x_1_line,
                y_0_line,
                y_1_line,
                y_2_line,
                y_3_line,
                y_label,
                label,
            )

    def _plot_single_cleavage_vertical(
        self,
        fig: go.Figure,
        start: int,
        end: int,
        isoform: str,
        cleavage_idx: int,
        pixels_per_cleavage: int,
        dy: int,
        x_0_line: int,
        x_1_line: int,
        x_2_line: int,
        label_plot_height: int,
        group_direction: int,
    ):
        if start == end:
            label = str(start)
            position = self.get_position_with_offset(start, isoform)
            y_0_line = (
                self.get_height() - position * self.PIXELS_PER_AA - self.SEQUENCE_OFFSET
            )
            y_0_line = self.offset_line_for_exon(
                y_0_line, start, self.FIGURE_ORIENTATION, isoform
            )
            y_1_line = (
                self.get_height()
                - cleavage_idx * pixels_per_cleavage
                - self.get_vertical_offset(dy)
            )
            x_3_line = (
                x_0_line
                + (label_plot_height - self.get_label_length(label)) * group_direction
            )
            x_label = (
                x_3_line + (self.get_label_length(label) // 2 + 5) * group_direction
            )

            self.plot_line_with_label_vertical(
                fig,
                x_0_line,
                x_1_line,
                x_2_line,
                x_3_line,
                y_0_line,
                y_1_line,
                x_label,
                label,
                False,
                None,
                None,
            )
        else:
            label = f"{start}-{end}"
            start_position = self.get_position_with_offset(start, isoform)
            end_position = self.get_position_with_offset(end, isoform)
            y_0_start_line = (
                self.get_height()
                - start_position * self.PIXELS_PER_AA
                - self.SEQUENCE_OFFSET
            )
            y_0_end_line = (
                self.get_height()
                - end_position * self.PIXELS_PER_AA
                - self.SEQUENCE_OFFSET
            )
            y_0_start_line = self.offset_line_for_exon(
                y_0_start_line, start, self.FIGURE_ORIENTATION, isoform
            )
            y_0_end_line = self.offset_line_for_exon(
                y_0_end_line, end, self.FIGURE_ORIENTATION, isoform
            )
            y_1_line = (
                self.get_height()
                - cleavage_idx * pixels_per_cleavage
                - self.get_vertical_offset(dy)
            )
            x_3_line = (
                x_0_line
                + (label_plot_height - self.get_label_length(label)) * group_direction
            )
            x_label = (
                x_3_line + (self.get_label_length(label) // 2 + 5) * group_direction
            )

            self.plot_range_with_label_vertical(
                fig,
                x_0_line,
                x_1_line,
                x_2_line,
                x_3_line,
                y_0_start_line,
                y_0_end_line,
                y_1_line,
                x_label,
                label,
            )

    def _plot_last_region_horizontal(
        self,
        fig: go.Figure,
        item_idx: int,
        last_i: int,
        first_item_in_region: int,
        pixels_per_item: int,
        dx: int,
        dy: int,
        mean_values: pd.DataFrame,
        y_0_groups: int,
        last_region: int,
        group_direction: int,
        vertical_space_left: int,
        is_ptm: bool,
    ):
        start_idx = item_idx - (last_i - first_item_in_region) - 1
        x_0_groups = start_idx * pixels_per_item + self.get_horizontal_offset(dx)
        region_length = len(mean_values.iloc[0:1, first_item_in_region:].columns)
        x_label = x_0_groups + (region_length * pixels_per_item) // 2 - dx // 2
        y_label = (
            y_0_groups
            + len(mean_values.index) * dy
            + (5 + self.get_label_height() // 2) * group_direction
        )

        self.plot_groups_horizontal(
            fig,
            mean_values.iloc[:, first_item_in_region:],
            x_0_groups,
            y_0_groups,
            dx,
            dy,
            x_label,
            y_label,
            last_region,
            group_direction,
            is_ptm,
        )

        self.create_custom_colorscale(
            fig,
            vertical_space_left,
            group_direction,
            x_0_groups,
            y_0_groups,
            region_length,
            pixels_per_item,
            is_ptm,
        )

    def _plot_last_region_vertical(
        self,
        fig: go.Figure,
        item_idx: int,
        last_i: int,
        first_item_in_region: int,
        pixels_per_item: int,
        dx: int,
        dy: int,
        mean_values: pd.DataFrame,
        x_0_groups: int,
        last_region: int,
        group_direction: int,
        horizontal_space_left: int,
        is_ptm: bool,
    ):
        start_idx = item_idx - (last_i - first_item_in_region) - 1
        y_0_groups = (
            self.get_height()
            - start_idx * pixels_per_item
            - self.get_vertical_offset(dy)
        )
        region_length = len(mean_values.iloc[0:1, first_item_in_region:].columns)
        y_label = y_0_groups - (region_length * pixels_per_item) // 2 + dy // 2
        x_label = (
            x_0_groups
            + len(mean_values.index) * dx
            + (5 + self.get_label_height() // 2) * group_direction
        )

        self.plot_groups_vertical(
            fig,
            mean_values.iloc[:, first_item_in_region:],
            x_0_groups,
            y_0_groups,
            dx,
            dy,
            x_label,
            y_label,
            last_region,
            group_direction,
            is_ptm,
        )

        self.create_custom_colorscale(
            fig,
            horizontal_space_left,
            group_direction,
            x_0_groups,
            y_0_groups,
            region_length,
            pixels_per_item,
            is_ptm,
        )

    def plot_cleavages(
        self,
        fig: go.Figure,
        cleavage_df: pd.DataFrame,
        pixels_per_cleavage: int,
        label_plot_height: int,
        above: str,
    ):
        # Validate and preprocess data
        mean_values, cleavage_region_type_map = self.preprocess_groups(cleavage_df)
        if len(cleavage_region_type_map) == 0 or len(mean_values.columns) == 0:
            logging.warning("No groups found in cleavage data, skipping cleavage plot.")
            return

        isoforms = cleavage_df.iloc[2:3, 2:].values.flatten().tolist()
        group_direction = 1 if above == "A" else -1

        # Inverse index for group B
        if above == "B":
            mean_values = mean_values.iloc[::-1]

        # Find longest label
        longest_label = ""
        for cleavage in list(zip(*cleavage_region_type_map))[0][::-1]:
            if self.get_label_length(str(cleavage)) > self.get_label_length(
                longest_label
            ):
                longest_label = str(cleavage)

        # Calculate line coordinates based on orientation
        if self.FIGURE_ORIENTATION == 0:
            y_0_line, y_1_line, y_2_line = (
                self._calculate_cleavage_line_coordinates_horizontal(
                    above, label_plot_height, longest_label, group_direction
                )
            )
            dx, dy, y_0_groups, vertical_space_left = (
                self._setup_cleavage_group_spacing_horizontal(
                    pixels_per_cleavage,
                    y_0_line,
                    label_plot_height,
                    above,
                    group_direction,
                    mean_values,
                )
            )
            self.plot_group_labels_horizontal(fig, mean_values, y_0_groups, dy)
        else:
            x_0_line, x_1_line, x_2_line = (
                self._calculate_cleavage_line_coordinates_vertical(
                    above, label_plot_height, longest_label, group_direction
                )
            )
            dx, dy, x_0_groups, horizontal_space_left = (
                self._setup_cleavage_group_spacing_vertical(
                    pixels_per_cleavage,
                    x_0_line,
                    label_plot_height,
                    above,
                    group_direction,
                    mean_values,
                )
            )
            self.plot_group_labels_vertical(fig, mean_values, x_0_groups, dx)

        # Initialize tracking variables
        first_cleavage_in_region = 0
        cleavage_idx = 0
        last_end = self.REGIONS[0][1]
        last_region_idx = 0
        previous_index = 0
        previous_region_type = None
        last_i = 0

        # Process each cleavage
        for i, (cleavage, region_type) in enumerate(cleavage_region_type_map):
            if "-" in str(cleavage):
                start, end = map(int, cleavage.split("-"))
            else:
                start = end = int(cleavage)

            # Check if we're entering a new region
            if (
                start > last_end
                or start < previous_index
                or (previous_region_type is not None and region_type != previous_region_type)
            ):
                # Plot groups for the previous region and add divider
                if self.FIGURE_ORIENTATION == 0:
                    self._plot_region_divider_horizontal(
                        fig,
                        cleavage_idx,
                        first_cleavage_in_region,
                        i,
                        pixels_per_cleavage,
                        dx,
                        dy,
                        mean_values,
                        y_0_groups,
                        last_region_idx,
                        group_direction,
                        False,
                    )
                else:
                    self._plot_region_divider_vertical(
                        fig,
                        cleavage_idx,
                        first_cleavage_in_region,
                        i,
                        pixels_per_cleavage,
                        dx,
                        dy,
                        mean_values,
                        x_0_groups,
                        last_region_idx,
                        group_direction,
                        False,
                    )

                # Update region tracking
                if start < previous_index:
                    last_region_idx += 1
                    last_end = self.REGIONS[last_region_idx][1]
                else:
                    while start > last_end:
                        last_region_idx += 1
                        last_end = self.REGIONS[last_region_idx][1]
                cleavage_idx += 1
                first_cleavage_in_region = i

            # Plot the current cleavage line and label
            if self.FIGURE_ORIENTATION == 0:
                self._plot_single_cleavage_horizontal(
                    fig,
                    start,
                    end,
                    isoforms[i],
                    cleavage_idx,
                    pixels_per_cleavage,
                    dx,
                    y_0_line,
                    y_1_line,
                    y_2_line,
                    label_plot_height,
                    group_direction,
                )
            else:
                self._plot_single_cleavage_vertical(
                    fig,
                    start,
                    end,
                    isoforms[i],
                    cleavage_idx,
                    pixels_per_cleavage,
                    dy,
                    x_0_line,
                    x_1_line,
                    x_2_line,
                    label_plot_height,
                    group_direction,
                )

            cleavage_idx += 1
            previous_index = start
            previous_region_type = region_type
            last_i = i

        # Finalize region tracking
        last_region_idx = len(self.REGIONS) - 1

        # Plot groups for the final region
        if self.FIGURE_ORIENTATION == 0:
            self._plot_last_region_horizontal(
                fig,
                cleavage_idx,
                last_i,
                first_cleavage_in_region,
                pixels_per_cleavage,
                dx,
                dy,
                mean_values,
                y_0_groups,
                last_region_idx,
                group_direction,
                vertical_space_left,
                False,
            )
        else:
            self._plot_last_region_vertical(
                fig,
                cleavage_idx,
                last_i,
                first_cleavage_in_region,
                pixels_per_cleavage,
                dx,
                dy,
                mean_values,
                x_0_groups,
                last_region_idx,
                group_direction,
                horizontal_space_left,
                False,
            )

    def get_horizontal_offset(self, dx):
        return self.calculate_group_space() + dx // 2

    def get_vertical_offset(self, dy):
        return self.calculate_group_space() + dy // 2

    def _calculate_ptm_line_coordinates_horizontal(
        self,
        above: str,
        label_plot_height: int,
        label_length: int,
        group_direction: int,
        second_row: bool,
    ):
        y_0_line = (
            self.SEQUENCE_BOUNDARIES["y1"]
            if above == "A"
            else self.SEQUENCE_BOUNDARIES["y0"]
        )
        y_1_line = y_0_line + 10 * group_direction
        y_2_line = (
            y_0_line
            + (
                label_plot_height
                - label_length
                - 10
                - self.plot_config.PTM_RECT_LENGTH
                - 10
            )
            * group_direction
        )
        if second_row:
            y_2_line = (
                y_0_line
                + (
                    label_plot_height
                    - 2 * (label_length + 10)
                    - self.plot_config.PTM_RECT_LENGTH
                    - 5
                )
                * group_direction
            )
        return y_0_line, y_1_line, y_2_line

    def _calculate_ptm_line_coordinates_vertical(
        self,
        above: str,
        label_plot_height: int,
        label_length: int,
        group_direction: int,
        second_row: bool,
    ):
        x_0_line = (
            self.SEQUENCE_BOUNDARIES["x1"]
            if above == "A"
            else self.SEQUENCE_BOUNDARIES["x0"]
        )
        x_1_line = x_0_line + 10 * group_direction
        x_2_line = (
            x_0_line
            + (
                label_plot_height
                - label_length
                - 10
                - self.plot_config.PTM_RECT_LENGTH
                - 10
            )
            * group_direction
        )
        if second_row:
            x_2_line = (
                x_0_line
                + (
                    label_plot_height
                    - 2 * (label_length + 10)
                    - self.plot_config.PTM_RECT_LENGTH
                    - 5
                )
                * group_direction
            )
        return x_0_line, x_1_line, x_2_line

    def _setup_ptm_group_spacing_horizontal(
        self,
        pixels_per_ptm: int,
        y_0_line: int,
        label_plot_height: int,
        above: str,
        group_direction: int,
        mean_values: pd.DataFrame,
    ):
        dx = pixels_per_ptm
        y_0_groups = y_0_line + (label_plot_height + 10) * group_direction
        vertical_space_left = (
            self.get_height() - y_0_groups if above == "A" else y_0_groups
        )
        # offset for region label
        dy_label = self.offset_region_label_from_angle()
        vertical_space_left -= dy_label * 2
        dy = vertical_space_left // len(mean_values.index) * group_direction
        return dx, dy, y_0_groups, vertical_space_left

    def _setup_ptm_group_spacing_vertical(
        self,
        pixels_per_ptm: int,
        x_0_line: int,
        label_plot_height: int,
        above: str,
        group_direction: int,
        mean_values: pd.DataFrame,
    ):
        dy = pixels_per_ptm
        x_0_groups = x_0_line + (label_plot_height + 10) * group_direction
        horizontal_space_left = (
            self.get_width() - x_0_groups if above == "A" else x_0_groups
        )
        # offset for region label
        dx_label = self.offset_region_label_from_angle()
        horizontal_space_left -= dx_label * 2
        dx = horizontal_space_left // len(mean_values.index) * group_direction
        return dx, dy, x_0_groups, horizontal_space_left

    def _update_region_tracking(
        self, ptm_position: int, previous_ptm: int, last_end: int, last_region: int
    ):
        if ptm_position < previous_ptm:
            last_region += 1
            last_end = self.REGIONS[last_region][1]
        else:
            while ptm_position > last_end:
                last_region += 1
                last_end = self.REGIONS[last_region][1]
        return last_end, last_region

    def _plot_single_ptm_horizontal(
        self,
        fig: go.Figure,
        ptm: str,
        ptm_position: int,
        isoform: str,
        ptm_idx: int,
        pixels_per_ptm: int,
        dx: int,
        i: int,
        second_row: bool,
        y_0_line: int,
        y_1_line: int,
        y_2_line: int,
        label_length: int,
        label_plot_height: int,
        group_direction: int,
        ptm_df: pd.DataFrame,
    ):
        position = self.get_position_with_offset(ptm_position, isoform)
        x_0_line = position * self.PIXELS_PER_AA + self.SEQUENCE_OFFSET
        x_0_line = self.offset_line_for_exon(
            x_0_line, ptm_position, self.FIGURE_ORIENTATION, isoform
        )
        x_1_line = ptm_idx * pixels_per_ptm + self.get_horizontal_offset(dx)
        y_3_line = y_2_line + 10 * group_direction
        if second_row and i % 2 == 1:
            y_3_line = y_2_line + (label_length + 10 + 5) * group_direction
        y_label = y_3_line + (self.get_label_length(ptm) + 10) // 2 * group_direction
        text_color = self.MODIFICATIONS[str(ptm_df.iloc[0, i + 2])][1]

        self.plot_line_with_label_horizontal(
            fig,
            x_0_line,
            x_1_line,
            y_0_line,
            y_1_line,
            y_2_line,
            y_3_line,
            y_label,
            ptm,
            True,
            text_color,
            str(ptm_df.iloc[0, i + 2]),
        )

        x_0_rect = x_1_line - dx // 2
        fig.add_shape(
            type="rect",
            x0=x_0_rect,
            x1=x_0_rect + dx,
            y0=y_0_line
            + (label_plot_height - self.plot_config.PTM_RECT_LENGTH) * group_direction,
            y1=y_0_line + label_plot_height * group_direction,
            fillcolor=text_color,
            line=dict(width=1, color="grey"),
            showlegend=False,
        )

    def _plot_single_ptm_vertical(
        self,
        fig: go.Figure,
        ptm: str,
        ptm_position: int,
        isoform: str,
        ptm_idx: int,
        pixels_per_ptm: int,
        dy: int,
        i: int,
        second_row: bool,
        x_0_line: int,
        x_1_line: int,
        x_2_line: int,
        label_length: int,
        label_plot_height: int,
        group_direction: int,
        ptm_df: pd.DataFrame,
    ):
        position = self.get_position_with_offset(ptm_position, isoform)
        y_0_line = (
            self.get_height() - position * self.PIXELS_PER_AA - self.SEQUENCE_OFFSET
        )
        y_0_line = self.offset_line_for_exon(
            y_0_line, ptm_position, self.FIGURE_ORIENTATION, isoform
        )
        y_1_line = (
            self.get_height() - ptm_idx * pixels_per_ptm - self.get_vertical_offset(dy)
        )
        x_3_line = x_2_line + 10 * group_direction
        if second_row and i % 2 == 1:
            x_3_line = x_2_line + (label_length + 10 + 5) * group_direction
        x_label = x_3_line + (self.get_label_length(ptm) + 10) // 2 * group_direction
        text_color = self.MODIFICATIONS[str(ptm_df.iloc[0, i + 2])][1]

        self.plot_line_with_label_vertical(
            fig,
            x_0_line,
            x_1_line,
            x_2_line,
            x_3_line,
            y_0_line,
            y_1_line,
            x_label,
            ptm,
            True,
            text_color,
            str(ptm_df.iloc[0, i + 2]),
        )

        y_0_rect = y_1_line - dy // 2
        fig.add_shape(
            type="rect",
            x0=x_0_line
            + (label_plot_height - self.plot_config.PTM_RECT_LENGTH) * group_direction,
            x1=x_0_line + label_plot_height * group_direction,
            y0=y_0_rect,
            y1=y_0_rect + dy,
            fillcolor=text_color,
            line=dict(width=1, color="grey"),
            showlegend=False,
        )

    def plot_ptms(
        self,
        fig: go.Figure,
        ptm_df: pd.DataFrame,
        pixels_per_ptm: int,
        label_plot_height: int,
        above: str,
        second_row: bool,
    ):
        # TODO: figure out why last region is mislabeled (and why this is not the case for cleavages)
        # TODO: fix problem that gap is nor properly accounted for - also for cleavages
        # Validate and preprocess data
        group_direction = 1 if above == "A" else -1
        mean_values, ptm_region_type_map = self.preprocess_groups(ptm_df)
        if len(mean_values) == 0:
            logging.warning("No groups found in PTM data, skipping PTM plot.")
            return
        if len(ptm_region_type_map) == 0:
            logging.warning("No PTMs to plot - will be omitted.")
            return

        isoforms = ptm_df.iloc[2:3, 2:].values.flatten().tolist()
        label_length = self.get_label_length(ptm_region_type_map[-1][0])

        # Inverse index for group B
        if above == "B":
            mean_values = mean_values.iloc[::-1]

        # Calculate line coordinates based on orientation
        if self.FIGURE_ORIENTATION == 0:
            y_0_line, y_1_line, y_2_line = (
                self._calculate_ptm_line_coordinates_horizontal(
                    above, label_plot_height, label_length, group_direction, second_row
                )
            )
            dx, dy, y_0_groups, vertical_space_left = (
                self._setup_ptm_group_spacing_horizontal(
                    pixels_per_ptm,
                    y_0_line,
                    label_plot_height,
                    above,
                    group_direction,
                    mean_values,
                )
            )
            self.plot_group_labels_horizontal(fig, mean_values, y_0_groups, dy)
        else:
            x_0_line, x_1_line, x_2_line = (
                self._calculate_ptm_line_coordinates_vertical(
                    above, label_plot_height, label_length, group_direction, second_row
                )
            )
            dx, dy, x_0_groups, horizontal_space_left = (
                self._setup_ptm_group_spacing_vertical(
                    pixels_per_ptm,
                    x_0_line,
                    label_plot_height,
                    above,
                    group_direction,
                    mean_values,
                )
            )
            self.plot_group_labels_vertical(fig, mean_values, x_0_groups, dx)

        # Initialize tracking variables
        last_end = self.REGIONS[0][1]
        first_ptm_in_region = 0
        ptm_idx = 0
        last_region_idx = 0
        previous_region_type_label = None
        previous_ptm_idx = 0
        last_i = 0

        for i, (ptm, region_type_label) in enumerate(ptm_region_type_map):
            ptm_position = int(ptm[1:])
            # Check if we're entering a new region
            if (
                ptm_position > last_end
                or ptm_position < previous_ptm_idx
                or (
                    previous_region_type_label is not None
                    and region_type_label != previous_region_type_label
                )
            ):
                # Plot groups for the previous region and add divider
                if self.FIGURE_ORIENTATION == 0:
                    self._plot_region_divider_horizontal(
                        fig,
                        ptm_idx,
                        first_ptm_in_region,
                        i,
                        pixels_per_ptm,
                        dx,
                        dy,
                        mean_values,
                        y_0_groups,
                        last_region_idx,
                        group_direction,
                        is_ptm=True,
                    )
                else:
                    self._plot_region_divider_vertical(
                        fig,
                        ptm_idx,
                        first_ptm_in_region,
                        i,
                        pixels_per_ptm,
                        dx,
                        dy,
                        mean_values,
                        x_0_groups,
                        last_region_idx,
                        group_direction,
                        is_ptm=True,
                    )

                # Update region tracking
                last_end, last_region_idx = self._update_region_tracking(
                    ptm_position, previous_ptm_idx, last_end, last_region_idx
                )
                ptm_idx += 1
                first_ptm_in_region = i

            # Plot the current PTM line, label, and marker box
            if self.FIGURE_ORIENTATION == 0:
                self._plot_single_ptm_horizontal(
                    fig,
                    ptm,
                    ptm_position,
                    isoforms[i],
                    ptm_idx,
                    pixels_per_ptm,
                    dx,
                    i,
                    second_row,
                    y_0_line,
                    y_1_line,
                    y_2_line,
                    label_length,
                    label_plot_height,
                    group_direction,
                    ptm_df,
                )
            else:
                self._plot_single_ptm_vertical(
                    fig,
                    ptm,
                    ptm_position,
                    isoforms[i],
                    ptm_idx,
                    pixels_per_ptm,
                    dy,
                    i,
                    second_row,
                    x_0_line,
                    x_1_line,
                    x_2_line,
                    label_length,
                    label_plot_height,
                    group_direction,
                    ptm_df,
                )

            ptm_idx += 1
            previous_ptm_idx = ptm_position
            previous_region_type_label = region_type_label
            last_i = i

        # Finalize region tracking
        last_region_idx = len(self.REGIONS) - 1

        # Plot groups for the final region
        if self.FIGURE_ORIENTATION == 0:
            self._plot_last_region_horizontal(
                fig,
                ptm_idx,
                last_i,
                first_ptm_in_region,
                pixels_per_ptm,
                dx,
                dy,
                mean_values,
                y_0_groups,
                last_region_idx,
                group_direction,
                vertical_space_left,
                is_ptm=True,
            )
        else:
            self._plot_last_region_vertical(
                fig,
                ptm_idx,
                last_i,
                first_ptm_in_region,
                pixels_per_ptm,
                dx,
                dy,
                mean_values,
                x_0_groups,
                last_region_idx,
                group_direction,
                horizontal_space_left,
                is_ptm=True,
            )

    def create_custom_colorscale(
        self,
        fig: go.Figure,
        vertical_space_left: int,
        group_direction: int,
        x_0_groups: int,
        y_0_groups: int,
        region_length: int,
        pixels_per_step: int,
        ptm: bool,
    ):
        if ptm:
            colorscale = [
                [0.0, self.plot_config.PTM_SCALE_COLOR_LOW],
                [0.5, self.plot_config.PTM_SCALE_COLOR_MID],
                [1.0, self.plot_config.PTM_SCALE_COLOR_HIGH],
            ]
            label = self.plot_config.PTM_LEGEND_TITLE
        else:
            colorscale = [
                [0.0, self.plot_config.CLEAVAGE_SCALE_COLOR_LOW],
                [0.5, self.plot_config.CLEAVAGE_SCALE_COLOR_MID],
                [1.0, self.plot_config.CLEAVAGE_SCALE_COLOR_HIGH],
            ]
            label = self.plot_config.CLEAVAGE_LEGEND_TITLE
        # Create a heatmap
        z = np.linspace(0, 1, 100).reshape(100, 1)

        if self.FIGURE_ORIENTATION == 0:
            dx = 15
            dy = 1
            scale_height = dy * 100 + 10 + self.get_label_height() * label.count("<br>")
            y_offset = (vertical_space_left - scale_height) // 2 * group_direction
            x_bar = x_0_groups + region_length * pixels_per_step + 10
            y_bar = y_0_groups + y_offset
            if group_direction == -1:
                y_bar -= scale_height
        else:
            z = z.T
            dx = 1
            dy = 15
            y_scale = dy + 5
            x_offset = vertical_space_left // 2 * group_direction
            y_bar = y_0_groups - region_length * pixels_per_step - 5
            x_bar = x_0_groups + x_offset - dx * 50
        fig.add_trace(
            go.Heatmap(
                x0=x_bar,
                y0=y_bar,
                z=z,
                dx=dx,
                dy=dy,
                colorscale=colorscale,
                showscale=False,
                hoverinfo="none",
            )
        )
        for i in range(3):
            percentage_label = f"{i * 50}%"
            if self.FIGURE_ORIENTATION == 0:
                x_scale = x_bar + 15 + self.get_label_length(percentage_label) // 2
                y_scale = y_bar + i * 100 * dy / 2
            else:
                x_scale = x_bar + i * 100 * dx / 2
                y_scale = y_bar - self.get_label_height()
            fig.add_annotation(
                x=x_scale,
                y=y_scale,
                text=percentage_label,
                showarrow=False,
                font=dict(
                    family=self.FONT,
                    size=self.sequence_plot_font_size,
                    color="black",
                ),
            )
        longest_label = ""
        for string in label.split("<br>"):
            if self.get_label_length(string) > self.get_label_length(longest_label):
                longest_label = string
        if self.FIGURE_ORIENTATION == 0:
            x_legend_title = x_bar + self.get_label_length(longest_label) // 2 - 15
            y_legend_title = y_bar + scale_height
        else:
            x_legend_title = x_bar + dx * 50
            y_legend_title = y_scale - self.get_label_height() * (
                label.count("<br>") + 1
            )
        fig.add_annotation(
            x=x_legend_title,
            y=y_legend_title,
            text=label,
            showarrow=False,
            font=dict(
                family=self.FONT,
                size=self.sequence_plot_font_size,
                color="black",
            ),
        )

    def filter_relevant_modification_sites(self, ptm_file: str, threshold: int):
        df = pd.read_csv(ptm_file, dtype={"ID": str, "Group": str})
        columns_to_keep = []
        for col in df.columns:
            if self.INCLUDED_MODIFICATIONS.get(df[col].iloc[0]):
                if df[col].iloc[1][:1] not in self.INCLUDED_MODIFICATIONS.get(
                    df[col].iloc[0]
                ):
                    continue
                if df[col].iloc[0] not in self.MODIFICATIONS:
                    continue
                if df[col].iloc[1][:1] == "R" and df[col].iloc[0] == "Deamidation":
                    df[col].iloc[0] = "Citrullination"
                columns_to_keep.append(col)
        df_filtered = df[columns_to_keep]
        df_values = df_filtered.iloc[3:].astype(int)
        sums = df_values.sum()
        filtered_columns = sums[sums >= threshold].index
        filtered_df = df[filtered_columns]

        result_df = pd.concat([df.iloc[:, :2], filtered_df], axis=1)

        return result_df

    def calculate_group_space(self):
        longest_label = ""
        for key in self.plot_config.GROUPS.keys():
            if self.get_label_length(key) > self.get_label_length(longest_label):
                longest_label = key
        return self.get_label_length(longest_label) + 10

    def calculate_legend_space(self, ptm: bool):
        if self.FIGURE_ORIENTATION == 0:
            longest_label = ""
            if ptm:
                for string in self.plot_config.PTM_LEGEND_TITLE.split("<br>"):
                    if self.get_label_length(string) > self.get_label_length(
                        longest_label
                    ):
                        longest_label = string
            else:
                for string in self.plot_config.CLEAVAGE_LEGEND_TITLE.split("<br>"):
                    if self.get_label_length(string) > self.get_label_length(
                        longest_label
                    ):
                        longest_label = string
            if self.get_label_length("100%") + 10 > self.get_label_length(
                longest_label
            ):
                return self.get_label_length("100%") + 10
            return self.get_label_length(longest_label)
        else:
            if ptm:
                title_height = self.get_label_height() * (
                    self.plot_config.PTM_LEGEND_TITLE.count("<br") + 1
                )
            else:
                title_height = self.get_label_height() * (
                    self.plot_config.CLEAVAGE_LEGEND_TITLE.count("<br") + 1
                )
            return self.get_label_height() + title_height + 10

    def get_present_mod_types(self):
        for above in self.plot_config.INPUT_FILES.values():
            if above[0] == "PTM":
                ptm_df = self.filter_relevant_modification_sites(
                    above[1], self.plot_config.MODIFICATION_THRESHOLD
                )
                return set(ptm_df.iloc[0:1, 2:].values.flatten().tolist())
        return set()

    def create_details_plot(self):
        legend = None
        messages = []

        present_mod_types = self.get_present_mod_types()
        mod_file = [
            f[1] for f in self.plot_config.INPUT_FILES.values() if f[0] == "PTM"
        ][0]
        detected_modifications = self.get_modifications_from_file(mod_file)
        if detected_modifications > present_mod_types:
            messages.append(
                {
                    "level": logging.WARNING,
                    "msg": "More modifications were detected than are present in the settings. Only the modifications "
                    "present in the modification settings are shown in the plot. You can see the additional "
                    'modifications in the "Tables" section.',
                }
            )

        if not "A" in self.plot_config.INPUT_FILES.keys():
            if self.plot_config.INPUT_FILES["B"][0] == "PTM":
                legend = "B"
            groups_missing = "A"
        elif not "B" in self.plot_config.INPUT_FILES.keys():
            if self.plot_config.INPUT_FILES["A"][0] == "PTM":
                legend = "A"
            groups_missing = "B"
        else:
            if self.plot_config.INPUT_FILES["A"][0] == "PTM":
                legend = "A"
            if self.plot_config.INPUT_FILES["B"][0] == "PTM":
                legend = "B"
            groups_missing = None

        fig = self._create_plot(
            input_file=self.input_file,
            present_modifications=present_mod_types,
            groups_missing=groups_missing,
            legend_positioning=legend,
            out_dir=self.output_path,
        )

        cleavage_df = None
        ptm_df = None
        relevant_groups = set()
        for above in self.plot_config.INPUT_FILES.keys():
            match self.plot_config.INPUT_FILES[above][0]:
                case "Cleavage":
                    cleavage_df = pd.read_csv(
                        self.plot_config.INPUT_FILES[above][1],
                        dtype={"ID": str, "Group": str},
                    )
                    cleavage_above = above
                    relevant_groups.update(cleavage_df["Group"].unique().tolist())
                case "PTM":
                    ptm_df = self.filter_relevant_modification_sites(
                        self.plot_config.INPUT_FILES[above][1],
                        self.plot_config.MODIFICATION_THRESHOLD,
                    )
                    ptm_above = above
                    relevant_groups.update(ptm_df["Group"].unique().tolist())
        assert (
            len(relevant_groups) > 0
        ), "No relevant groups found in the provided data."
        # Since we are using only the relevant parts of the metadata, we have to filter the groups here. The original
        # tool was not developed with variable metadata in mind, so this has to be corrected for
        for group in self.plot_config.GROUPS.copy():
            if group not in relevant_groups:
                self.plot_config.GROUPS.pop(group)

        if self.FIGURE_ORIENTATION == 0:
            plot_space = self.get_width() - self.SEQUENCE_BOUNDARIES["x0"]
        else:
            # first we calculate the missing space above the sequence and then subtract it from the total height
            plot_space = self.get_height() - (
                self.get_height() - self.SEQUENCE_BOUNDARIES["y0"]
            )

        label_plot_height = 150

        if cleavage_df is not None:
            present_regions = self.get_present_regions_cleavage(cleavage_df)
            number_of_cleavages = len(cleavage_df.columns)
            number_of_dividers = present_regions.count(True) - 1
            cleavage_space = (
                plot_space
                - self.calculate_legend_space(False)
                - self.calculate_group_space()
            )
            pixels_per_cleavage = cleavage_space // (
                number_of_cleavages + number_of_dividers
            )
            assert pixels_per_cleavage >= self.FONT_SIZE

            self.plot_cleavages(
                fig, cleavage_df, pixels_per_cleavage, label_plot_height, cleavage_above
            )

        if ptm_df is not None:
            present_regions = self.get_present_regions_ptm(ptm_df)
            number_of_ptms = len(ptm_df.columns)
            number_of_dividers = present_regions.count(True) - 1
            second_row = False
            ptm_space = (
                plot_space
                - self.calculate_legend_space(True)
                - self.calculate_group_space()
            )
            pixels_per_ptm = ptm_space // (number_of_ptms + number_of_dividers)
            if (
                number_of_ptms + number_of_dividers
            ) * self.get_label_height() > 2 * ptm_space:
                raise ValueError("Too many PTMs to fit in plot")
            if (
                number_of_ptms + 2 * number_of_dividers
            ) * self.get_label_height() > ptm_space:
                second_row = True

            self.plot_ptms(
                fig, ptm_df, pixels_per_ptm, label_plot_height, ptm_above, second_row
            )

        self.finalize_plotting(
            fig,
            self.output_path,
            save_plot=self.plot_config.SAVE_PLOT,
            show_plot=self.plot_config.SHOW_PLOT,
        )
        return fig, messages
