"""Module for creating main sequence plot."""

import os
import importlib
from pathlib import Path

import plotly.graph_objects as go
from protein_sequencing import utils, exon_helper

CONFIG = importlib.import_module('configs.default_config', 'configs')


def create_plot(
        input_file: str | os.PathLike,
        present_modifications,
        groups_missing=None,
        legend_positioning=None,
        out_dir=None
) -> go.Figure:
    """Create the plot with main sequence and all addiational information."""
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
    ) = exon_helper.retrieve_exon(input_file, CONFIG.MIN_EXON_LENGTH, out_dir=Path(out_dir))

    # exon checks
    if exon_found:
        # get exon lengths
        utils.EXON_1_OFFSET['index_start'] = exon_start_index
        utils.EXON_1_OFFSET['index_end'] = exon_start_index + exon_1_length - 1
        utils.EXON_2_OFFSET['index_start'] = exon_start_index
        utils.EXON_2_OFFSET['index_end'] = exon_start_index + exon_2_length - 1

        # calculate new max sequence length with exons
        max_sequence_length = max_sequence_length - max_exon_length + exon_1_length + exon_2_length

        # check if exon lengths match with regions
        region_end_matches_exon = False
        for i, region in enumerate(CONFIG.REGIONS):
            if region[1] + 1 == exon_start_index:
                region_end_matches_exon = True
                if len(CONFIG.REGIONS) < i + 2:
                    raise ValueError(f"Exon start {exon_start_index} matches a region end for region {region}, \
                                     but there are not enough regions after it, please check your supplied region list.")

                exon_1_region = CONFIG.REGIONS[i + 1]
                exon_2_region = CONFIG.REGIONS[i + 2]
                if exon_1_region[1] - region[1] != exon_1_length:
                    if exon_1_region[1] - region[1] != exon_2_length:
                        raise ValueError(
                            f"Exon 1 length {exon_1_length} does not match with end for region {exon_1_region}.")

                    # Swap regions in case they are in the wrong order
                    exon_1_region, exon_2_region = exon_2_region, exon_1_region
                    exon_1_length, exon_2_length = exon_2_length, exon_1_length
                if exon_2_region[1] - region[1] != exon_2_length:
                    raise ValueError(
                        f"Exon 2 length {exon_2_length} does not match with end for region {exon_2_region}.")
        if not region_end_matches_exon:
            raise ValueError(
                f"Exon start {exon_start_index} does not match any region end, please check your supplied region list.")

    # basis for all pixel calculations
    if CONFIG.FIGURE_ORIENTATION == 0:
        max_sequence_length_pixels = utils.get_width() - utils.get_left_margin()
        utils.PIXELS_PER_AA = int(
            (max_sequence_length_pixels - CONFIG.EXONS_GAP * exon_found * 2) // max_sequence_length)
        utils.SEQUENCE_OFFSET = utils.get_left_margin()
    else:
        max_sequence_length_pixels = utils.get_height() - utils.get_top_margin()
        utils.PIXELS_PER_AA = int(
            (max_sequence_length_pixels - CONFIG.EXONS_GAP * exon_found * 2) // max_sequence_length)
        utils.SEQUENCE_OFFSET = utils.get_top_margin()

    # calculate region boundaries in pixels
    region_boundaries = []
    region_end_pixel = utils.SEQUENCE_OFFSET
    region_start = 1
    exon_offset = 0
    region_index = 0
    # 0 = normal region, 1 = region before exon, 2 = region after start exon
    # 3 = end exon/ region after exon, 4 = start exon, 5 = middle exon
    region_plot_type = 0
    while region_index < len(CONFIG.REGIONS):
        region_name, region_end, region_group, _ = CONFIG.REGIONS[region_index]
        region_start_pixel = region_end_pixel
        region_end_pixel = region_end * utils.PIXELS_PER_AA + 1 + utils.SEQUENCE_OFFSET
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
                if region_index + 1 == len(CONFIG.REGIONS) - 1:
                    region_plot_type = 3
                # if next region is not last region and also not start exon
                elif region_index > 0:
                    region_plot_type = 5
                first_exon_offset = CONFIG.EXONS_GAP
                utils.EXON_1_OFFSET['pixel_start'] = region_start_pixel + first_exon_offset
                utils.EXON_1_OFFSET['pixel_end'] = region_end_pixel + first_exon_offset
                region_boundaries.append(
                    (region_name, region_start_pixel + first_exon_offset, region_end_pixel + first_exon_offset,
                     CONFIG.SEQUENCE_REGION_COLORS[region_group], region_start, region_end, region_plot_type))
                exon_offset = exon_1_length * utils.PIXELS_PER_AA + CONFIG.EXONS_GAP
                exon_1_region_end = region_end
                # process next exon
                region_index += 1
                region_name, region_end, region_group, _ = CONFIG.REGIONS[region_index]
                region_start_pixel = region_end_pixel + CONFIG.EXONS_GAP * 2
                region_end_pixel = region_end * utils.PIXELS_PER_AA + 1 + utils.SEQUENCE_OFFSET + exon_offset
                utils.EXON_2_OFFSET['pixel_start'] = region_start_pixel
                utils.EXON_2_OFFSET['pixel_end'] = region_end_pixel
                region_boundaries.append(
                    (region_name, region_start_pixel, region_end_pixel, CONFIG.SEQUENCE_REGION_COLORS[region_group],
                     region_start, region_end, region_plot_type))
                region_start = max(exon_1_region_end, region_end) + 1
                region_index += 1
                if region_plot_type == 4:
                    region_plot_type = 2
                else:
                    region_plot_type = 3
                continue

        region_boundaries.append((region_name, region_start_pixel + exon_offset, region_end_pixel + exon_offset,
                                  CONFIG.SEQUENCE_REGION_COLORS[region_group], region_start, region_end, 0))
        region_start = region_end + 1
        region_index += 1
        region_plot_type = 0

    fig = create_sequence_plot(region_boundaries, present_modifications, groups_missing, legend_positioning)

    return fig


def create_sequence_plot(region_boundaries: list[tuple[str, int, int, str, int, int]], present_modifications,
                         groups_missing: str | None, legend_positioning: str | None) -> go.Figure:
    """Create the sequence plot."""
    fig = go.Figure()

    width = utils.get_width()
    height = utils.get_height()

    # General Layout
    fig.update_layout(
        title="",
        width=width,
        height=height,
        xaxis=dict(range=[0, width], autorange=False),
        yaxis=dict(range=[0, height], autorange=False),
        plot_bgcolor="white",
        font_family=CONFIG.FONT,
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

        labels = [CONFIG.MODIFICATIONS[mod][0] for mod in present_modifications] + [CONFIG.MODIFICATION_LEGEND_TITLE]
        sorted_labels = sorted(labels, key=sort_key)

        if not groups_missing:
            if legend_positioning == 'A':
                x_legend = 0 if CONFIG.FIGURE_ORIENTATION == 0 else width // 2 + CONFIG.SEQUENCE_PLOT_HEIGHT // 2
                y_legend = height // 2 - CONFIG.SEQUENCE_PLOT_HEIGHT // 2 + (
                            len(present_modifications) + 1) * utils.get_label_height() if CONFIG.FIGURE_ORIENTATION == 0 else height
            else:
                x_legend = 0 if CONFIG.FIGURE_ORIENTATION == 0 else width // 2 - CONFIG.SEQUENCE_PLOT_HEIGHT // 2
                y_legend = height // 2 + CONFIG.SEQUENCE_PLOT_HEIGHT // 2 if CONFIG.FIGURE_ORIENTATION == 0 else height
        else:
            x_legend = 0
            y_legend = 0
            if groups_missing == 'A':
                if CONFIG.FIGURE_ORIENTATION == 1:
                    x_legend = width
                    y_legend = height
                else:
                    y_legend = height
            if groups_missing == 'B':
                if CONFIG.FIGURE_ORIENTATION == 1:
                    y_legend = height
                else:
                    y_legend = (len(present_modifications) + 1) * utils.get_label_height()

        text_position = "bottom right"
        if legend_positioning == 'A' and CONFIG.FIGURE_ORIENTATION == 1:
            text_position = "bottom left"
        fig.add_trace(go.Scatter(x=[x_legend], y=[y_legend],
                                 mode='text',
                                 text=f"<b>{CONFIG.MODIFICATION_LEGEND_TITLE}</b>",
                                 textposition=text_position,
                                 showlegend=False, hoverinfo='none',
                                 textfont=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                                               color="black")))
        y_legend -= utils.get_label_height()

        labels = [CONFIG.MODIFICATIONS[mod] for mod in present_modifications]
        sorted_labels = sorted(labels, key=sort_key)
        if groups_missing == 'A' or CONFIG.FIGURE_ORIENTATION == 1:
            sorted_labels = sorted_labels[::-1]
        for i, mod in enumerate(sorted_labels):
            fig.add_trace(
                go.Scatter(x=[x_legend],
                           y=[y_legend - i * utils.get_label_height()],
                           mode='text',
                           text=mod[0],
                           textposition=text_position,
                           showlegend=False,
                           hoverinfo='none',
                           textfont=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE, color=mod[1])))

    # Sequence
    fig = plot_sequence(fig, region_boundaries, groups_missing)

    return fig


def plot_sequence(fig, region_boundaries, groups_missing):
    """Plot the sequence with regions and exons."""
    sequence_x0, sequence_y0 = 0, 0
    x0, x1, y0, y1 = 0, 0, 0, 0
    if CONFIG.FIGURE_ORIENTATION == 0:
        y0 = utils.get_height() // 2 - CONFIG.SEQUENCE_PLOT_HEIGHT // 2
        if groups_missing:
            if groups_missing == 'A':
                y0 = utils.get_height() - CONFIG.SEQUENCE_PLOT_HEIGHT
            if groups_missing == 'B':
                y0 = 0
        y1 = y0 + CONFIG.SEQUENCE_PLOT_HEIGHT
    else:
        x0 = utils.get_width() // 2 - CONFIG.SEQUENCE_PLOT_HEIGHT // 2
        if groups_missing:
            if groups_missing == 'B':
                x0 = 0
            if groups_missing == 'A':
                x0 = utils.get_width() - CONFIG.SEQUENCE_PLOT_HEIGHT
        x1 = x0 + CONFIG.SEQUENCE_PLOT_HEIGHT
    last_region_end = 0
    last_i = 0
    for i, (region_name, region_start_pixel, region_end_pixel, region_color, region_start, region_end,
            exon_type) in enumerate(region_boundaries):
        if CONFIG.FIGURE_ORIENTATION == 0:
            x0 = region_start_pixel
            x1 = region_end_pixel
        else:
            y0 = utils.get_height() - region_start_pixel
            y1 = utils.get_height() - region_end_pixel

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
            if CONFIG.FIGURE_ORIENTATION == 0:
                x = [x0, x1, x1 + CONFIG.EXONS_GAP // 2, x0, x0]
                y = [y0, y0, y1, y1, y0]
            else:
                x = [x0, x1, x1, x0, x0]
                y = [y0, y0, y1, y1 - CONFIG.EXONS_GAP // 2, y0]
        elif exon_type == 2:
            if CONFIG.FIGURE_ORIENTATION == 0:
                x = [x0, x1, x1, x0 - CONFIG.EXONS_GAP // 2, x0]
                y = [y0, y0, y1, y1, y0]
            else:
                x = [x0, x1, x1, x0, x0]
                y = [y0 + CONFIG.EXONS_GAP // 2, y0, y1, y1, y0 + CONFIG.EXONS_GAP // 2]
        elif exon_type == 3:
            if CONFIG.FIGURE_ORIENTATION == 0:
                x = [x0 - CONFIG.EXONS_GAP // 2, x1, x1, x0, x0 - CONFIG.EXONS_GAP // 2]
                y = [y0, y0, y1, y1, y0]
            else:
                x = [x0, x1, x1, x0, x0]
                y = [y0, y0 + CONFIG.EXONS_GAP // 2, y1, y1, y0]
        elif exon_type == 4:
            if CONFIG.FIGURE_ORIENTATION == 0:
                x = [x0, x1 + CONFIG.EXON_GAP // 2, x1, x0, x0]
                y = [y0, y0, y1, y1, y0]
            else:
                x = [x0, x1, x1, x0, x0]
                y = [y0, y0, y1 - CONFIG.EXON_GAP // 2, y1, y0]
        elif exon_type == 5:
            if CONFIG.FIGURE_ORIENTATION == 0:
                x = [x0 - CONFIG.EXON_GAP // 2, x1, x1 + CONFIG.EXON_GAP // 2, x0, x0 - CONFIG.EXON_GAP // 2]
                y = [y0, y0, y1, y1, y0]
            else:
                x = [x0, x1, x1, x0, x0]
                y = [y0, y0 + CONFIG.EXON_GAP // 2, y1, y1 - CONFIG.EXON_GAP // 2, y0]
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
            font=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE, color="black"),
            textangle=90 if CONFIG.FIGURE_ORIENTATION == 1 else 0
        )
        if i == 0:
            if CONFIG.FIGURE_ORIENTATION == 0:
                x = x0 - utils.get_label_length(str(region_start))
                y = y_label
            else:
                x = x_label
                y = y0 + utils.get_label_height()
            fig.add_annotation(
                x=x,
                y=y,
                text='1',
                showarrow=False,
                font=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE, color="gray"),
                textangle=0
            )
            sequence_x0, sequence_y0 = x0, y0
        last_i = i
        last_region_end = region_end
    if CONFIG.FIGURE_ORIENTATION == 0:
        x = x1 + utils.get_label_length(str(last_region_end))
        y = y_label
    else:
        x = x_label
        y = y1 - utils.get_label_height()
    fig.add_annotation(
        x=x,
        y=y,
        text=max(last_region_end, region_boundaries[last_i - 1][5]),
        showarrow=False,
        font=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE, color="gray"),
        textangle=0
    )

    utils.SEQUENCE_BOUNDARIES = {'x0': sequence_x0, 'x1': x1, 'y0': sequence_y0, 'y1': y1}
    return fig
