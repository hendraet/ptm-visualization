import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

from protein_sequencing import exon_helper


class Plotter:
    def __init__(self, config, plot_config):
        self.plot_config = plot_config

        self.SEQUENCE_BOUNDARIES = {'x0': 0, 'x1': 0, 'y0': 0, 'y1': 0}
        self.PIXELS_PER_AA = 0
        self.SEQUENCE_OFFSET = 0
        self.EXON_1_OFFSET = {'index_start': -1,
                              'index_end': -1,
                              'pixel_start': -1,
                              'pixel_end': -1}
        self.EXON_2_OFFSET = {'index_start': -1,
                              'index_end': -1,
                              'pixel_start': -1,
                              'pixel_end': -1}

        # TODO: maybe lower case these
        self.REGIONS = config.REGIONS
        self.MODIFICATIONS = config.MODIFICATIONS
        self.MODIFICATION_LEGEND_TITLE = config.MODIFICATION_LEGEND_TITLE
        self.FIGURE_ORIENTATION = config.FIGURE_ORIENTATION
        self.FIGURE_WIDTH = config.FIGURE_WIDTH
        self.FIGURE_HEIGHT = config.FIGURE_HEIGHT
        self.FONT_SIZE = config.FONT_SIZE
        self.FONT = config.FONT
        self.EXONS_GAP = config.EXONS_GAP
        self.PTMS_TO_HIGHLIGHT = config.PTMS_TO_HIGHLIGHT
        self.PTM_HIGHLIGHT_LABEL_COLOR = config.PTM_HIGHLIGHT_LABEL_COLOR
        self.INCLUDED_MODIFICATIONS = config.INCLUDED_MODIFICATIONS

        self.regions = config.REGIONS
        self.sequence_region_colors = config.SEQUENCE_REGION_COLORS
        self.min_exon_length = config.MIN_EXON_LENGTH
        self.font = config.FONT
        self.sequence_plot_height = config.SEQUENCE_PLOT_HEIGHT
        self.sequence_plot_font_size = config.SEQUENCE_PLOT_FONT_SIZE

    def get_width(self):
        """Return width of the plot, based on user settings in default_config.py."""
        if self.FIGURE_ORIENTATION == 0:
            return self.FIGURE_WIDTH
        return self.FIGURE_HEIGHT

    def get_height(self):
        """Return height of the plot, based on user settings in default_config.py."""
        if self.FIGURE_ORIENTATION == 0:
            return self.FIGURE_HEIGHT
        return self.FIGURE_WIDTH

    def get_left_margin(self):
        """Return the left margin for the sequence plot.
        Calculated based on the longest label in the legend."""
        longest_text = self.MODIFICATION_LEGEND_TITLE
        for mod in self.MODIFICATIONS:
            if len(self.MODIFICATIONS[mod][0]) > len(longest_text):
                longest_text = self.MODIFICATIONS[mod][0]
        return int((self.get_label_length(longest_text) / self.get_width() * 1.05) * self.get_width())

    def get_top_margin(self):
        """Return the top margin for the sequence plot.
        Calculated based on the number of modifications in the legend."""
        legend_height = (len(self.MODIFICATIONS) + 3) * self.get_label_height()
        return int((legend_height / self.get_height() * 1.05) * self.get_height())

    def get_label_length(self, label):
        """Approximate the length of a label in pixels based on font size and label length."""
        return int(self.FONT_SIZE / 1.5 * len(label))

    def get_label_height(self):
        """Approximate the height of a label in pixels based on font size."""
        return self.FONT_SIZE + self.FONT_SIZE // 5

    def separate_by_group(self, groups_by_position_and_isoform):
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

    def different_possibilities_plot(self, width: int, height: int, different_possibilities: list[int]):
        """Debug option. Plot the different possibilities of the sequence in a heatmap."""
        rectangle = np.zeros((height, width))
        for i, value in enumerate(different_possibilities):
            rectangle[:, i] = value
        fig = go.Figure(data=go.Heatmap(z=rectangle))
        fig.show()

    def clean_up(self, out_dir):
        """Remove all files from the temporary directory."""
        out_path = Path(out_dir)
        for file_path in out_path.iterdir():
            if file_path.is_file():
                file_path.unlink()

    def finalize_plotting(self, fig, output_path, save_plot: bool = True, show_plot: bool = True):
        """Show the plot and save it as a .png and .svg file."""
        # TODO: hardcoded paths -.-
        if save_plot:
            output_svg = f"{output_path}/figure1.svg"
            output_png = f"{output_path}/figure1.png"
            fig.write_image(output_png)
            fig.write_image(output_svg)
        if show_plot:
            fig.show()

        self.clean_up(output_path)

    def get_position_with_offset(self, position, isoform):
        """Return the position in the rendering index based on sequence position and isoform."""
        if isoform == 'exon2':
            position += self.EXON_1_OFFSET['index_end'] - self.EXON_1_OFFSET['index_start'] + 1
        elif position > max(self.EXON_1_OFFSET['index_end'], self.EXON_2_OFFSET['index_end']):
            if isoform != 'general':
                raise ValueError(f"Position {position} is out of range for isoform {isoform}")
            exon_1_length = self.EXON_1_OFFSET['index_end'] - self.EXON_1_OFFSET['index_start'] + 1
            exon_2_length = self.EXON_2_OFFSET['index_end'] - self.EXON_2_OFFSET['index_start'] + 1
            position += max(exon_1_length, exon_2_length)

        return position

    def offset_line_for_exon(self, line_position, aa_position, oritentation):
        """Offset the line position based on the exon boundaries."""
        if aa_position >= self.EXON_1_OFFSET['index_start'] and self.EXON_1_OFFSET['index_start'] != -1:
            if oritentation == 0:
                line_position += self.EXONS_GAP
            else:
                line_position -= self.EXONS_GAP
        if aa_position > self.EXON_1_OFFSET['index_end'] and self.EXON_1_OFFSET['index_start'] != -1:
            if oritentation == 0:
                line_position += self.EXONS_GAP
            else:
                line_position -= self.EXONS_GAP

        return line_position

    @staticmethod
    def get_modifications_from_file(mod_file: str) -> set:
        with open(mod_file, 'r', encoding="utf-8") as f:
            rows = f.readlines()[1:4]
            modification_types = rows[0].strip().split(',')
            present_modifications = set()
            for i, (label) in enumerate(rows[1].strip().split(',')):
                if label == '':
                    continue
                present_modifications.add(modification_types[i])
        return present_modifications

    def _create_plot(
            self,
            input_file: str | os.PathLike,
            present_modifications: list[str] | None,
            groups_missing=None,
            legend_positioning=None,
            out_dir=None,
    ) -> go.Figure:
        """Create the plot with main sequence and all additional information."""

        (
            exon_found,
            exon_start_index,
            _,
            max_exon_length,
            _,
            exon_1_length,
            _,
            exon_2_length,
            _,
            max_sequence_length
        ) = exon_helper.retrieve_exon(Path(input_file), self.min_exon_length, out_dir=Path(out_dir))

        # exon checks
        if exon_found:
            # get exon lengths
            self.EXON_1_OFFSET['index_start'] = exon_start_index
            self.EXON_1_OFFSET['index_end'] = exon_start_index + exon_1_length - 1
            self.EXON_2_OFFSET['index_start'] = exon_start_index
            self.EXON_2_OFFSET['index_end'] = exon_start_index + exon_2_length - 1

            # calculate new max sequence length with exons
            max_sequence_length = max_sequence_length - max_exon_length + exon_1_length + exon_2_length

            # check if exon lengths match with regions
            region_end_matches_exon = False
            for i, region in enumerate(self.regions):
                if region[1] + 1 == exon_start_index:
                    region_end_matches_exon = True
                    if len(self.regions) <= i + 2:
                        raise ValueError(f"Exon start {exon_start_index} matches a region end for region {region}, but "
                                         "there are not enough regions after it, please check your supplied region "
                                         "list.")

                    exon_1_region = self.regions[i + 1]
                    exon_2_region = self.regions[i + 2]
                    if exon_1_region[1] - region[1] != exon_1_length:
                        if exon_1_region[1] - region[1] != exon_2_length:
                            raise ValueError(f"Exon 1 length {exon_1_length} does not match with end for "
                                             f"region {exon_1_region}.")

                        # Swap regions in case they are in the wrong order
                        exon_1_region, exon_2_region = exon_2_region, exon_1_region
                        exon_1_length, exon_2_length = exon_2_length, exon_1_length
                    if exon_2_region[1] - region[1] != exon_2_length:
                        raise ValueError(
                            f"Exon 2 length {exon_2_length} does not match with end for region {exon_2_region}.")
            if not region_end_matches_exon:
                raise ValueError(
                    f"Exon start {exon_start_index} does not match any region end, please check your supplied region "
                    "list - maybe it is missing some regions or it doesn't match the provided fasta sequence.")

        # basis for all pixel calculations
        if self.FIGURE_ORIENTATION == 0:
            max_sequence_length_pixels = (self.get_width() - self.get_left_margin())
            self.PIXELS_PER_AA = int((max_sequence_length_pixels - self.EXONS_GAP * exon_found * 2) // max_sequence_length)
            self.SEQUENCE_OFFSET = self.get_left_margin()
        else:
            max_sequence_length_pixels = self.get_height() - self.get_top_margin()
            self.PIXELS_PER_AA = int(
                (max_sequence_length_pixels - self.EXONS_GAP * exon_found * 2) // max_sequence_length)
            self.SEQUENCE_OFFSET = self.get_top_margin()

        # calculate region boundaries in pixels
        region_boundaries = []
        region_end_pixel = self.SEQUENCE_OFFSET
        region_start = 1
        exon_offset = 0
        region_index = 0
        # 0 = normal region, 1 = region before exon, 2 = region after start exon
        # 3 = end exon/ region after exon, 4 = start exon, 5 = middle exon
        region_plot_type = 0
        while region_index < len(self.regions):
            region_name, region_end, region_group, _ = self.regions[region_index]
            region_start_pixel = region_end_pixel
            region_end_pixel = region_end * self.PIXELS_PER_AA + 1 + self.SEQUENCE_OFFSET
            if exon_found:
                if region_end == exon_1_region[1]:
                    # alter last boundary to include exon
                    if region_index > 0:
                        last_boundary = region_boundaries[-1]
                        region_boundaries[-1] = (last_boundary[0], last_boundary[1], last_boundary[2], last_boundary[3],
                                                 last_boundary[4], last_boundary[5], 1)
                    else:
                        region_plot_type = 4

                    # add current exon
                    # if next region is also last region
                    if region_index + 1 == len(self.regions) - 1:
                        region_plot_type = 3
                    # if next region is not last region and also not start exon
                    elif region_index > 0:
                        region_plot_type = 5
                    first_exon_offset = self.EXONS_GAP
                    self.EXON_1_OFFSET['pixel_start'] = region_start_pixel + first_exon_offset
                    self.EXON_1_OFFSET['pixel_end'] = region_end_pixel + first_exon_offset
                    region_boundaries.append(
                        (region_name, region_start_pixel + first_exon_offset, region_end_pixel + first_exon_offset,
                         self.sequence_region_colors[region_group], region_start, region_end, region_plot_type))
                    exon_offset = exon_1_length * self.PIXELS_PER_AA + self.EXONS_GAP
                    exon_1_region_end = region_end
                    # process next exon
                    region_index += 1
                    region_name, region_end, region_group, _ = self.regions[region_index]
                    region_start_pixel = region_end_pixel + self.EXONS_GAP * 2
                    region_end_pixel = region_end * self.PIXELS_PER_AA + 1 + self.SEQUENCE_OFFSET + exon_offset
                    self.EXON_2_OFFSET['pixel_start'] = region_start_pixel
                    self.EXON_2_OFFSET['pixel_end'] = region_end_pixel
                    region_boundaries.append(
                        (region_name, region_start_pixel, region_end_pixel, self.sequence_region_colors[region_group],
                         region_start, region_end, region_plot_type))
                    region_start = max(exon_1_region_end, region_end) + 1
                    region_index += 1
                    if region_plot_type == 4:
                        region_plot_type = 2
                    else:
                        region_plot_type = 3
                    continue

            region_boundaries.append((region_name, region_start_pixel + exon_offset, region_end_pixel + exon_offset,
                                      self.sequence_region_colors[region_group], region_start, region_end, 0))
            region_start = region_end + 1
            region_index += 1
            region_plot_type = 0

        fig = self._create_sequence_plot(
            region_boundaries,
            present_modifications,
            groups_missing=groups_missing,
            legend_positioning=legend_positioning,
        )
        return fig

    def _create_sequence_plot(
            self,
            region_boundaries: list[tuple[str, int, int, str, int, int]],
            present_modifications,
            groups_missing: str | None,
            legend_positioning: str | None,
    ) -> go.Figure:
        """Create the sequence plot."""
        fig = go.Figure()

        width = self.get_width()
        height = self.get_height()

        # General Layout
        fig.update_layout(
            title="",
            width=width,
            height=height,
            xaxis=dict(range=[0, width], autorange=False),
            yaxis=dict(range=[0, height], autorange=False),
            plot_bgcolor="white",
            font_family=self.font,
            margin=dict(l=0, r=0, t=0, b=0),
        )
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)

        # Legend
        if legend_positioning:
            def sort_key(mod):
                """Sort modifications by length and then by number of specific characters."""
                chars_to_count = ['i', 'l', 't', 'j']
                primary_key = len(mod[0])
                secondary_key = -sum(mod[0].count(char) for char in chars_to_count)
                return (primary_key, secondary_key)

            labels = [self.MODIFICATIONS[mod][0] for mod in present_modifications] + [self.MODIFICATION_LEGEND_TITLE]
            sorted_labels = sorted(labels, key=sort_key)

            if not groups_missing:
                if legend_positioning == 'A':
                    x_legend = 0 if self.FIGURE_ORIENTATION == 0 else width // 2 + self.sequence_plot_height // 2
                    y_legend = (
                        height // 2 -
                        self.sequence_plot_height // 2 +
                        (len(present_modifications) + 1) * self.get_label_height()
                        if self.FIGURE_ORIENTATION == 0 else height
                    )
                else:
                    x_legend = 0 if self.FIGURE_ORIENTATION == 0 else width // 2 - self.sequence_plot_height // 2
                    y_legend = height // 2 + self.sequence_plot_height // 2 if self.FIGURE_ORIENTATION == 0 else height
            else:
                x_legend = 0
                y_legend = 0
                if groups_missing == 'A':
                    if self.FIGURE_ORIENTATION == 1:
                        x_legend = width
                        y_legend = height
                    else:
                        y_legend = height
                if groups_missing == 'B':
                    if self.FIGURE_ORIENTATION == 1:
                        y_legend = height
                    else:
                        y_legend = (len(present_modifications) + 1) * self.get_label_height()

            text_position = "bottom right"
            if legend_positioning == 'A' and self.FIGURE_ORIENTATION == 1:
                text_position = "bottom left"
            fig.add_trace(go.Scatter(x=[x_legend], y=[y_legend],
                                     mode='text',
                                     text=f"<b>{self.MODIFICATION_LEGEND_TITLE}</b>",
                                     textposition=text_position,
                                     showlegend=False, hoverinfo='none',
                                     textfont=dict(size=self.sequence_plot_font_size, color="black")))
            y_legend -= self.get_label_height()

            labels = [self.MODIFICATIONS[mod] for mod in present_modifications]
            sorted_labels = sorted(labels, key=sort_key)
            if groups_missing == 'A' or self.FIGURE_ORIENTATION == 1:
                sorted_labels = sorted_labels[::-1]
            for i, mod in enumerate(sorted_labels):
                fig.add_trace(
                    go.Scatter(x=[x_legend],
                               y=[y_legend - i * self.get_label_height()],
                               mode='text',
                               text=mod[0],
                               textposition=text_position,
                               showlegend=False,
                               hoverinfo='none',
                               textfont=dict(size=self.sequence_plot_font_size, color=mod[1])))
        # Sequence
        fig = self._plot_sequence(fig, region_boundaries, groups_missing)
        return fig

    def _plot_sequence(
            self,
            fig,
            region_boundaries,
            groups_missing,
    ):
        """Plot the sequence with regions and exons."""
        sequence_x0, sequence_y0 = 0, 0
        x0, x1, y0, y1 = 0, 0, 0, 0
        if self.FIGURE_ORIENTATION == 0:
            y0 = self.get_height() // 2 - self.sequence_plot_height // 2
            if groups_missing:
                if groups_missing == 'A':
                    y0 = self.get_height() - self.sequence_plot_height
                if groups_missing == 'B':
                    y0 = 0
            y1 = y0 + self.sequence_plot_height
        else:
            x0 = self.get_width() // 2 - self.sequence_plot_height // 2
            if groups_missing:
                if groups_missing == 'B':
                    x0 = 0
                if groups_missing == 'A':
                    x0 = self.get_width() - self.sequence_plot_height
            x1 = x0 + self.sequence_plot_height
        last_region_end = 0
        last_i = 0
        for i, (region_name, region_start_pixel, region_end_pixel, region_color, region_start, region_end,
                exon_type) in enumerate(region_boundaries):
            if self.FIGURE_ORIENTATION == 0:
                x0 = region_start_pixel
                x1 = region_end_pixel
            else:
                y0 = self.get_height() - region_start_pixel
                y1 = self.get_height() - region_end_pixel

            # 0 = normal region, 1 = region before exon, 2 = region after start exon
            # 3 = end exon/ region after exon, 4 = start exon, 5 = middle exon
            if exon_type == 0:
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
            elif exon_type == 1:
                if self.FIGURE_ORIENTATION == 0:
                    x = [x0, x1, x1 + self.EXONS_GAP // 2, x0, x0]
                    y = [y0, y0, y1, y1, y0]
                else:
                    x = [x0, x1, x1, x0, x0]
                    y = [y0, y0, y1, y1 - self.EXONS_GAP // 2, y0]
            elif exon_type == 2:
                if self.FIGURE_ORIENTATION == 0:
                    x = [x0, x1, x1, x0 - self.EXONS_GAP // 2, x0]
                    y = [y0, y0, y1, y1, y0]
                else:
                    x = [x0, x1, x1, x0, x0]
                    y = [y0 + self.EXONS_GAP // 2, y0, y1, y1, y0 + self.EXONS_GAP // 2]
            elif exon_type == 3:
                if self.FIGURE_ORIENTATION == 0:
                    x = [x0 - self.EXONS_GAP // 2, x1, x1, x0, x0 - self.EXONS_GAP // 2]
                    y = [y0, y0, y1, y1, y0]
                else:
                    x = [x0, x1, x1, x0, x0]
                    y = [y0, y0 + self.EXONS_GAP // 2, y1, y1, y0]
            elif exon_type == 4:
                if self.FIGURE_ORIENTATION == 0:
                    x = [x0, x1 + self.EXONS_GAP // 2, x1, x0, x0]
                    y = [y0, y0, y1, y1, y0]
                else:
                    x = [x0, x1, x1, x0, x0]
                    y = [y0, y0, y1 - self.EXONS_GAP // 2, y1, y0]
            elif exon_type == 5:
                if self.FIGURE_ORIENTATION == 0:
                    x = [x0 - self.EXONS_GAP // 2, x1, x1 + self.EXONS_GAP // 2, x0, x0 - self.EXONS_GAP // 2]
                    y = [y0, y0, y1, y1, y0]
                else:
                    x = [x0, x1, x1, x0, x0]
                    y = [y0, y0 + self.EXONS_GAP // 2, y1, y1 - self.EXONS_GAP // 2, y0]
            if exon_type != 0:
                fig.add_trace(go.Scatter(x=x,
                                         y=y,
                                         mode='lines',
                                         fillcolor=region_color,
                                         fill='toself',
                                         line=dict(color="darkgrey", width=2), showlegend=False, hoverinfo='none'))
            # Labels
            x_label = (x0 + x1) / 2
            y_label = (y0 + y1) / 2
            fig.add_annotation(
                x=x_label,
                y=y_label,
                text=region_name,
                showarrow=False,
                font=dict(size=self.sequence_plot_font_size, color="black"),
                textangle=90 if self.FIGURE_ORIENTATION == 1 else 0
            )
            if i == 0:
                if self.FIGURE_ORIENTATION == 0:
                    x = x0 - self.get_label_length(str(region_start))
                    y = y_label
                else:
                    x = x_label
                    y = y0 + self.get_label_height()
                fig.add_annotation(
                    x=x,
                    y=y,
                    text='1',
                    showarrow=False,
                    font=dict(size=self.sequence_plot_font_size, color="gray"),
                    textangle=0
                )
                sequence_x0, sequence_y0 = x0, y0
            last_i = i
            last_region_end = region_end
        if self.FIGURE_ORIENTATION == 0:
            x = x1 + self.get_label_length(str(last_region_end))
            y = y_label
        else:
            x = x_label
        fig.add_annotation(
            x=x,
            y=y,
            text=max(last_region_end, region_boundaries[last_i - 1][5]),
            showarrow=False,
            font=dict(size=self.sequence_plot_font_size, color="gray"),
            textangle=0
        )

        self.SEQUENCE_BOUNDARIES = {'x0': sequence_x0, 'x1': x1, 'y0': sequence_y0, 'y1': y1}
        return fig
