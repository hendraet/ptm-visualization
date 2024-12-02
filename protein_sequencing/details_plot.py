from collections import defaultdict
import math
import os
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import importlib

from protein_sequencing import utils, sequence_plot

CONFIG = importlib.import_module('configs.default_config', 'configs')
PLOT_CONFIG = importlib.import_module('configs.default_details', 'configs')

def get_present_regions(positions, isoforms):
    ranges = []
    for range in positions:
        if '-' in str(range):
            start, end = map(int, range.split('-'))
            ranges.append((start, end))
        else:
            start = end = int(range)
            ranges.append((start, end))

    region_ranges = []
    region_start = 1
    for _, region_end, _, _ in CONFIG.REGIONS:
        region_ranges.append((region_start, region_end))
        region_start = region_end + 1

    regions_present = [False] * len(region_ranges)
    region_index = 0
    for i, range in enumerate(ranges):
        if isoforms[i] == 'exon1':
            index = next((index for index, region in enumerate(CONFIG.REGIONS) if region[1] == utils.EXON_1_OFFSET["index_end"]), None)
            if index:
                regions_present[index] = True
        elif isoforms[i] == 'exon2':
            index = next((index for index, region in enumerate(CONFIG.REGIONS) if region[1] == utils.EXON_2_OFFSET["index_end"]), None)
            if index:
                regions_present[index] = True
        while range[0] > region_ranges[region_index][1]:
            region_index += 1
        regions_present[region_index] = True
    return regions_present

def get_present_regions_cleavage(cleavage_df: pd.DataFrame):
    cleavages = cleavage_df.iloc[1:2,2:].values[0].tolist()
    isoforms = cleavage_df.iloc[2:3,2:].values[0].tolist()
    return get_present_regions(cleavages, isoforms)

def get_present_regions_ptm(ptm_df: pd.DataFrame):
    ptms = ptm_df.iloc[1:2,2:].values[0].tolist()
    ptms = [ptm[1:] for ptm in ptms]
    isoforms = ptm_df.iloc[2:3,2:].values[0].tolist()
    return get_present_regions(ptms, isoforms)
    

def plot_line_with_label_horizontal(fig: go.Figure, x_0: int, x_1: int, y_0: int, y_1: int, y_2: int, y_3: int, y_label: int, label: str, ptm: bool, ptm_color: str | None = None, ptm_modification: str | None = None):
    line_color = "black"
    if ptm:
        line_color = ptm_color
    fig.add_trace(go.Scatter(x=[x_0, x_0, x_1, x_1],
                            y=[y_0, y_1, y_2, y_3],
                            mode='lines',
                            line=dict(color=line_color, width=1), showlegend=False, hoverinfo='none'))
    if ptm:
        color=ptm_color
        if f'{ptm_modification}({label[0]})@{label[1:]}' in CONFIG.PTMS_TO_HIGHLIGHT:
            fig.add_shape(type='rect',
                            x0 = x_1-utils.get_label_height()//2-1,
                            x1 = x_1+utils.get_label_height()//2+1,
                            y0 = y_label-utils.get_label_length(label)//2-3,
                            y1 = y_label+utils.get_label_length(label)//2+3,
                            line=dict(width=0),
                            fillcolor=CONFIG.PTM_HIGHLIGHT_LABEL_COLOR,
                            showlegend=False,)
    else:
        color = PLOT_CONFIG.CLEAVAGE_LABEL_COLOR
        if label in PLOT_CONFIG.CLEAVAGES_TO_HIGHLIGHT:
            color = PLOT_CONFIG.CLEAVAGE_HIGHLIGHT_COLOR
    fig.add_annotation(x=x_1, y=y_label,
                        text=label,
                        showarrow=False,
                        textangle=-90,
                        font=dict(
                            family=CONFIG.FONT,
                            size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                            color=color,
                            ))
    return fig

def plot_line_with_label_vertical(fig: go.Figure, x_0: int, x_1: int, x_2: int, x_3: int, y_0: int, y_1: int, x_label: int, label: str, ptm: bool, ptm_color: str | None = None, ptm_modification: str | None = None):
    line_color = "black"
    if ptm:
        line_color = ptm_color
    fig.add_trace(go.Scatter(x=[x_0, x_1, x_2, x_3],
                            y=[y_0, y_0, y_1, y_1],
                            mode='lines',
                            line=dict(color=line_color, width=1), showlegend=False, hoverinfo='none'))
    if ptm:
        color=ptm_color
        if f'{ptm_modification}({label[0]})@{label[1:]}' in CONFIG.PTMS_TO_HIGHLIGHT:
            fig.add_shape(type='rect',
                            x0 = x_label-utils.get_label_length(label)//2-3,
                            x1 = x_label+utils.get_label_length(label)//2+3,
                            y0 = y_1-utils.get_label_height()//2-1,
                            y1 = y_1+utils.get_label_height()//2+1,
                            line=dict(width=0),
                            fillcolor=CONFIG.PTM_HIGHLIGHT_LABEL_COLOR,
                            showlegend=False,)
    else:
        color = PLOT_CONFIG.CLEAVAGE_LABEL_COLOR
        if label in PLOT_CONFIG.CLEAVAGES_TO_HIGHLIGHT:
            color = PLOT_CONFIG.CLEAVAGE_HIGHLIGHT_COLOR
    fig.add_annotation(x=x_label, y=y_1,
                        text=label,
                        showarrow=False,
                        font=dict(
                            family=CONFIG.FONT,
                            size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                            color=color,
                            ))
    return fig

def plot_range_with_label_horizontal(fig: go.Figure, x_0_start: int, x_0_end: int, x_1: int, y_0: int, y_1: int, y_2: int, y_3: int, y_label: int, label: str):
    fig.add_trace(go.Scatter(x=[x_0_start, x_0_start, x_1, x_1, x_1, x_0_end, x_0_end],
                            y=[y_0, y_1, y_2, y_3, y_2, y_1, y_0],
                            mode='lines',
                            fill='toself',
                            line=dict(color="black", width=1), showlegend=False, hoverinfo='none'))
    
    color = PLOT_CONFIG.CLEAVAGE_LABEL_COLOR
    if label in PLOT_CONFIG.CLEAVAGES_TO_HIGHLIGHT:
        color = PLOT_CONFIG.CLEAVAGE_HIGHLIGHT_COLOR
    fig.add_annotation(x=x_1, y=y_label,
                        text=label,
                        showarrow=False,
                        textangle=-90,
                        font=dict(
                            family=CONFIG.FONT,
                            size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                            color=color,
                            ))
    return fig

def plot_range_with_label_vertical(fig: go.Figure, x_0: int, x_1: int, x_2: int, x_3: int, y_0_start: int, y_0_end: int, y_1: int, x_label: int, label: str):
    fig.add_trace(go.Scatter(x=[x_0, x_1, x_2, x_3, x_2, x_1, x_0],
                            y=[y_0_start, y_0_start, y_1, y_1, y_1, y_0_end, y_0_end],
                            mode='lines',
                            fill='toself',
                            line=dict(color="black", width=1), showlegend=False, hoverinfo='none'))
    color = PLOT_CONFIG.CLEAVAGE_LABEL_COLOR
    if label in PLOT_CONFIG.CLEAVAGES_TO_HIGHLIGHT:
        color = PLOT_CONFIG.CLEAVAGE_HIGHLIGHT_COLOR
    fig.add_annotation(x=x_label, y=y_1,
                        text=label,
                        showarrow=False,
                        font=dict(
                            family=CONFIG.FONT,
                            size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                            color=color,
                            ))
    return fig

def plot_groups_horizontal(fig: go.Figure, df: pd.DataFrame, x_0_groups: int, y_0_groups: int, dx: int, dy: int, x_label: int, y_label: int, last_region: int, group_dircetion: int, ptm: bool):
    x_margin = 0
    if dx % 2 != 0:
        x_margin = 1
    color_low = PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_LOW
    color_mid = PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_MID
    color_high = PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_HIGH
    if ptm:
        color_low = PLOT_CONFIG.PTM_SCALE_COLOR_LOW
        color_mid = PLOT_CONFIG.PTM_SCALE_COLOR_MID
        color_high = PLOT_CONFIG.PTM_SCALE_COLOR_HIGH
    fig.add_shape(type='rect',
                    x0=x_0_groups - dx//2 - x_margin,
                    y0=y_0_groups,
                    x1=x_0_groups + dx * len(df.iloc[0:1,:].columns) - dx//2,
                    y1=y_0_groups + dy * len(df.index) + 1,
                    fillcolor='grey',
                    line=dict(color='grey', width=1),
                    showlegend=False,
                    layer='below',)
    fig.add_trace(go.Heatmap(z=df,
                    x0=x_0_groups,
                    y0=y_0_groups+dy//2,
                    dx=dx, dy=dy,
                    showscale=False, hoverinfo='none',
                    xgap=1, ygap=1,
                    zmin=0,
                    zmax=1,
                    zmid=0.5,
                    colorscale=[[0, color_low], [0.5, color_mid], [1, color_high]]))
    yanchor = 'bottom'
    xanchor = 'left'
    if group_dircetion == -1:
        yanchor = 'top'
        xanchor = 'right'
    fig.add_annotation(x=x_label-utils.get_label_height()*group_dircetion, y=y_label-int((8/offset_region_label_from_angle())*group_dircetion),    
            text=CONFIG.REGIONS[last_region][3],
            showarrow=False,
            textangle=-PLOT_CONFIG.REGION_LABEL_ANGLE_GROUPS,
            xanchor=xanchor,
            yanchor=yanchor,
            font=dict(
                family=CONFIG.FONT,
                size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                color='black',
                ))
    return fig

def plot_groups_vertical(fig: go.Figure, df: pd.DataFrame, x_0_groups: int, y_0_groups: int, dx: int, dy: int, x_label: int, y_label: int, last_region: int, group_dircetion: int, ptm: bool):
    y_margin = 0
    if dy % 2 != 0:
        y_margin = 1
    color_low = PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_LOW
    color_mid = PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_MID
    color_high = PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_HIGH
    if ptm:
        color_low = PLOT_CONFIG.PTM_SCALE_COLOR_LOW
        color_mid = PLOT_CONFIG.PTM_SCALE_COLOR_MID
        color_high = PLOT_CONFIG.PTM_SCALE_COLOR_HIGH
    fig.add_shape(type='rect',
                    x0=x_0_groups,
                    y0=y_0_groups + dy//2 + y_margin,
                    x1=x_0_groups + dx * len(df.index) + 1,
                    y1=y_0_groups - dy * len(df.iloc[0:1,:].columns) + dy//2,
                    fillcolor='grey',
                    line=dict(color='grey', width=1),
                    showlegend=False,
                    layer='below',)
    
    fig.add_trace(go.Heatmap(z=df.T,
                    x0=x_0_groups+dx//2,
                    y0=y_0_groups,
                    dx=dx, dy=-dy,
                    showscale=False, hoverinfo='none',
                    xgap=1, ygap=1,
                    zmin=0,
                    zmax=1,
                    zmid=0.5,
                    colorscale=[[0, color_low], [0.5, color_mid], [1, color_high]]))
    xanchor = 'left'
    if group_dircetion == -1:
        xanchor = 'right'
    fig.add_annotation(x=x_label, y=y_label,
            text=CONFIG.REGIONS[last_region][3],
            showarrow=False,
            textangle=-PLOT_CONFIG.REGION_LABEL_ANGLE_GROUPS+90,
            xanchor=xanchor,
            font=dict(
                family=CONFIG.FONT,
                size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                color='black',
                ))
    return fig

def plot_group_labels_horizontal(fig: go.Figure, mean_values: pd.DataFrame, y_0_groups: int, dy: int, dx: int):
    for i, group in enumerate(mean_values.index):
        y_0_rect = y_0_groups + i*dy
        x_1_rect = utils.SEQUENCE_OFFSET - dx//2
        fig.add_shape(type='rect',
                      x0 = 0,
                      x1 = x_1_rect,
                      y0 = y_0_rect,
                      y1 = y_0_rect + dy,
                      fillcolor=PLOT_CONFIG.GROUPS[group][1],
                      line=dict(width=0),
                      showlegend=False,
                      layer='below',)
        color = get_label_color(group)

        fig.add_annotation(x=x_1_rect//2, y=y_0_rect + dy//2,
            text=group,
            showarrow=False,
            align='center',
            font=dict(
                family=CONFIG.FONT,
                size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                color=color))
def get_label_color(group: str):
    # based on https://stackoverflow.com/questions/3942878/
    red, green, blue = tuple(int(PLOT_CONFIG.GROUPS[group][1][i:i+2], 16) for i in (1, 3, 5))
    return '#000000' if red*0.299 + green*0.587 + blue*0.114 > 130 else '#ffffff'

def plot_group_labels_vertical(fig: go.Figure, mean_values: pd.DataFrame, x_0_groups: int, dx: int, dy: int):
    for i, group in enumerate(mean_values.index):
        x_0_rect = x_0_groups + i*dx
        y_0_rect = utils.get_height()
        y_rect = utils.SEQUENCE_OFFSET - dy//2
        fig.add_shape(type='rect',
                      x0 = x_0_rect,
                      x1 = x_0_rect + dx,
                      y0 = y_0_rect,
                      y1 = y_0_rect - y_rect,
                      fillcolor=PLOT_CONFIG.GROUPS[group][1],
                      line=dict(width=0),
                      showlegend=False,
                      layer='below',)
        
        color = get_label_color(group)

        fig.add_annotation(x=x_0_rect + dx//2, y=y_0_rect - y_rect//2,
            text=group,
            showarrow=False,
            align='center',
            textangle=90,
            font=dict(
                family=CONFIG.FONT,
                size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                color=color))
        
def preprocess_groups(df: pd.DataFrame):
    df.columns = df.iloc[0]
    labels = df.iloc[1:2,2:].values.flatten().tolist()
    df = df.iloc[3:]

    reverse_group_mapping = {}
    for key, list in PLOT_CONFIG.GROUPS.items():
        values = list[0]
        for value in values:
            reverse_group_mapping[value] = key
    df.iloc[:,1] = df.iloc[:,1].map(reverse_group_mapping)
    mean_values = df.iloc[:,2:].astype(float).groupby(df.iloc[:,1]).mean()
    mean_values = mean_values.reindex([*PLOT_CONFIG.GROUPS])

    return mean_values, labels

def offset_region_label_from_angle():
    longest_label = ''
    for i, (_, _, _, region_label_short) in enumerate(CONFIG.REGIONS):
        if utils.get_label_length(region_label_short) > utils.get_label_length(longest_label):
            longest_label = region_label_short
        
    length = utils.get_label_length(longest_label)
    height = utils.get_label_height()

    angle_radians = math.radians(-PLOT_CONFIG.REGION_LABEL_ANGLE_GROUPS)
    dy = abs((length / 2) * math.sin(angle_radians)) + abs((height / 2) * math.cos(angle_radians))
    
    return int(dy)+10


def plot_cleavages(fig: go.Figure, cleavage_df: pd.DataFrame, pixels_per_cleavage: int, label_plot_height: int, above: str):
    # prepare data:
    mean_values, cleavages = preprocess_groups(cleavage_df)
    isoforms = cleavage_df.iloc[2:3,2:].values.flatten().tolist()
    if above == 'B':
        mean_values = mean_values.iloc[::-1]

    longest_label = ''
    for cleavage in cleavages[::-1]:
        if utils.get_label_length(str(cleavage)) > utils.get_label_length(longest_label):
            longest_label = str(cleavage)

    group_direction = 1 if above == 'A' else -1
    first_cleavage_in_region = 0
    cleavage_idx = 0
    last_end = CONFIG.REGIONS[0][1]
    last_region = 0

    if CONFIG.FIGURE_ORIENTATION == 0:
        y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if above == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
        y_1_line = y_0_line + 10 * group_direction
        y_2_line = y_0_line + (label_plot_height - utils.get_label_length(longest_label) - 10) * group_direction

        y_0_groups = y_0_line + (label_plot_height + 10) * group_direction
        vertical_space_left = utils.get_height() - y_0_groups if above == 'A' else y_0_groups
        # offset for border around heatmap
        vertical_space_left -= 2
        # offset for label for region
        dy_label = offset_region_label_from_angle()
        vertical_space_left -= dy_label*2
        dx = pixels_per_cleavage
        dy = vertical_space_left//len(mean_values.index)*group_direction

        plot_group_labels_horizontal(fig, mean_values, y_0_groups, dy, dx)
    else:
        x_0_line = utils.SEQUENCE_BOUNDARIES['x1'] if above == 'A' else utils.SEQUENCE_BOUNDARIES['x0']
        x_1_line = x_0_line + 10 * group_direction
        x_2_line = x_0_line + (label_plot_height - utils.get_label_length(longest_label) - 10) * group_direction

        x_0_groups = x_0_line + (label_plot_height + 10) * group_direction
        horizontal_space_left = utils.get_width() - x_0_groups if above == 'A' else x_0_groups
        # offset for border around heatmap
        horizontal_space_left -= 2
        # offset for label for region
        dx_label = offset_region_label_from_angle()
        horizontal_space_left -= dx_label*2
        dy = pixels_per_cleavage
        dx = horizontal_space_left//len(mean_values.index)*group_direction

        plot_group_labels_vertical(fig, mean_values, x_0_groups, dx, dy)
        
    previous_index = 0
    for i, cleavage in enumerate(cleavages):
        if '-' in str(cleavage):
            start, end = map(int, cleavage.split('-'))
        else:
            start = end = int(cleavage)
        if start > last_end or start < previous_index:
            if CONFIG.FIGURE_ORIENTATION == 0:
                start_idx = cleavage_idx - (i - first_cleavage_in_region)
                x_0_groups = start_idx * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                x_divider = cleavage_idx * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                x_label = x_0_groups + (x_divider-x_0_groups)//2 - dx//2
                y_label = y_0_groups + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction

                plot_groups_horizontal(fig, mean_values.iloc[:,first_cleavage_in_region:i], x_0_groups, y_0_groups, dx, dy, x_label, y_label, last_region, group_direction, False)                
                
                fig.add_trace(go.Scatter(x=[x_divider,x_divider],
                            y=[y_0_groups, y_0_groups+len(mean_values.index)*dy],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
            else:
                start_idx = cleavage_idx - (i - first_cleavage_in_region)
                y_0_groups = utils.get_height() - start_idx * pixels_per_cleavage - utils.SEQUENCE_OFFSET
                y_divider = utils.get_height() - cleavage_idx * pixels_per_cleavage - utils.SEQUENCE_OFFSET
                y_label = y_0_groups - (y_0_groups - y_divider)//2 + dy//2
                x_label = x_0_groups + len(mean_values.index)*dx + (5+utils.get_label_height()//2) * group_direction

                plot_groups_vertical(fig, mean_values.iloc[:,first_cleavage_in_region:i], x_0_groups, y_0_groups, dx, dy, x_label, y_label, last_region, group_direction, False)
                
                fig.add_trace(go.Scatter(x=[x_0_groups, x_0_groups+len(mean_values.index)*dx],
                            y=[y_divider, y_divider],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
            if start < previous_index:
                last_region += 1
                last_end = CONFIG.REGIONS[last_region][1]
            else:    
                while start > last_end:
                    last_region += 1
                    last_end = CONFIG.REGIONS[last_region][1]
            cleavage_idx += 1
            first_cleavage_in_region = i
        if CONFIG.FIGURE_ORIENTATION == 0:
            if start == end:
                label = str(start)
                position = utils.get_position_with_offset(start, isoforms[i])
                x_0_line = position * utils.PIXELS_PER_AA + utils.SEQUENCE_OFFSET
                x_0_line = utils.offset_line_for_exon(x_0_line, start, CONFIG.FIGURE_ORIENTATION)
                x_1_line = cleavage_idx * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                y_3_line = y_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction
                
                plot_line_with_label_horizontal(fig,
                                 x_0_line, x_1_line,
                                 y_0_line, y_1_line, y_2_line, y_3_line,
                                 y_label,
                                 label, False, None, None)
            else:
                label = f'{start}-{end}'
                start_position = utils.get_position_with_offset(start, isoforms[i])
                end_position = utils.get_position_with_offset(end, isoforms[i])
                x_0_start_line = start_position * utils.PIXELS_PER_AA + utils.SEQUENCE_OFFSET
                x_0_end_line = end_position * utils.PIXELS_PER_AA + utils.SEQUENCE_OFFSET
                x_0_start_line = utils.offset_line_for_exon(x_0_start_line, start, CONFIG.FIGURE_ORIENTATION)
                x_0_end_line = utils.offset_line_for_exon(x_0_end_line, end, CONFIG.FIGURE_ORIENTATION)
                x_1_line = cleavage_idx * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                y_3_line = y_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction

                plot_range_with_label_horizontal(fig,
                                    x_0_start_line, x_0_end_line, x_1_line,
                                    y_0_line, y_1_line, y_2_line, y_3_line,
                                    y_label,
                                    label)
        else:
            if start == end:
                label = str(start)
                position = utils.get_position_with_offset(start, isoforms[i])
                y_0_line = utils.get_height() - position * utils.PIXELS_PER_AA - utils.SEQUENCE_OFFSET
                y_0_line = utils.offset_line_for_exon(y_0_line, start, CONFIG.FIGURE_ORIENTATION)
                y_1_line = utils.get_height() - cleavage_idx * pixels_per_cleavage - utils.SEQUENCE_OFFSET
                x_3_line = x_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                x_label = x_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction

                plot_line_with_label_vertical(fig,
                                 x_0_line, x_1_line, x_2_line, x_3_line,
                                 y_0_line, y_1_line,
                                 x_label,
                                 label, False, None, None)
            else:
                label = f'{start}-{end}'
                start_position = utils.get_position_with_offset(start, isoforms[i])
                end_position = utils.get_position_with_offset(end, isoforms[i])
                y_0_start_line = utils.get_height() - start_position * utils.PIXELS_PER_AA - utils.SEQUENCE_OFFSET
                y_0_end_line = utils.get_height() - end_position * utils.PIXELS_PER_AA - utils.SEQUENCE_OFFSET
                y_0_start_line = utils.offset_line_for_exon(y_0_start_line, start, CONFIG.FIGURE_ORIENTATION)
                y_0_end_line = utils.offset_line_for_exon(y_0_end_line, end, CONFIG.FIGURE_ORIENTATION)
                y_1_line = utils.get_height() - cleavage_idx * pixels_per_cleavage - utils.SEQUENCE_OFFSET
                x_3_line = x_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                x_label = x_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction

                plot_range_with_label_vertical(fig,
                                    x_0_line, x_1_line, x_2_line, x_3_line,
                                    y_0_start_line, y_0_end_line,
                                    y_1_line,
                                    x_label,
                                    label)
        cleavage_idx += 1
        previous_index = start
      
    while start > last_end:
        last_region += 1
        last_end = CONFIG.REGIONS[last_region][1]

    if isoforms[first_cleavage_in_region] == 'exon2' and isoforms[first_cleavage_in_region-1] != 'exon1':
        last_region += 1
        
    # plot groups for last region
    if CONFIG.FIGURE_ORIENTATION == 0:
        start_idx = cleavage_idx - (i - first_cleavage_in_region)-1
        x_0_groups = start_idx * pixels_per_cleavage + utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_cleavage_in_region:].columns)
        x_label = x_0_groups + (region_length * pixels_per_cleavage)//2 - dx//2
        y_label = y_0_groups + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction
        plot_groups_horizontal(fig, mean_values.iloc[:,first_cleavage_in_region:], x_0_groups, y_0_groups, dx, dy, x_label, y_label, last_region, group_direction, False)

        create_custome_colorscale(fig, vertical_space_left, group_direction, x_0_groups, y_0_groups, region_length, pixels_per_cleavage, False)
    else:
        start_idx = cleavage_idx - (i - first_cleavage_in_region)-1
        y_0_groups = utils.get_height() - start_idx * pixels_per_cleavage - utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_cleavage_in_region:].columns)
        y_label = y_0_groups - (region_length * pixels_per_cleavage)//2 + dy//2
        x_label = x_0_groups + len(mean_values.index)*dx + (5+utils.get_label_height()//2) * group_direction
        plot_groups_vertical(fig, mean_values.iloc[:,first_cleavage_in_region:], x_0_groups, y_0_groups, dx, dy, x_label, y_label, last_region, group_direction, False)

        create_custome_colorscale(fig, horizontal_space_left, group_direction, x_0_groups, y_0_groups, region_length, pixels_per_cleavage, False)

def plot_ptms(fig: go.Figure, ptm_df: pd.DataFrame, pixels_per_ptm: int, label_plot_height: int, above: str, second_row: bool):
    group_direction = 1 if above == 'A' else -1
    mean_values, ptms = preprocess_groups(ptm_df)
    isoforms = ptm_df.iloc[2:3,2:].values.flatten().tolist()
    # For debugging purposes
    #new_row = pd.DataFrame([ptms], columns=mean_values.columns, index=['ptm_position'])
    #mean_values = pd.concat([new_row, mean_values])
    #mean_values.to_csv('plotting_data_ptms.csv', sep=',')

    label_length = utils.get_label_length(ptms[-1])
    # inverse index for group B
    if above == 'B':
        mean_values = mean_values.iloc[::-1]

    if CONFIG.FIGURE_ORIENTATION == 0:
        y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if above == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
        y_1_line = y_0_line + 10 * group_direction
        y_2_line = y_0_line + (label_plot_height - label_length - 10 - PLOT_CONFIG.PTM_RECT_LENGTH - 10) * group_direction
        if second_row:
            y_2_line = y_0_line + (label_plot_height - 2*(label_length + 10) - PLOT_CONFIG.PTM_RECT_LENGTH - 5) * group_direction
    else:
        x_0_line = utils.SEQUENCE_BOUNDARIES['x1'] if above == 'A' else utils.SEQUENCE_BOUNDARIES['x0']
        x_1_line = x_0_line + 10 * group_direction
        x_2_line = x_0_line + (label_plot_height - label_length - 10 - PLOT_CONFIG.PTM_RECT_LENGTH - 10) * group_direction
        if second_row:
            x_2_line = x_0_line + (label_plot_height - 2*(label_length + 10) - PLOT_CONFIG.PTM_RECT_LENGTH - 5) * group_direction

    last_end = CONFIG.REGIONS[0][1]
    first_ptm_in_region = 0
    ptm_idx = 0
    last_region = 0

    if CONFIG.FIGURE_ORIENTATION == 0:
        dx = pixels_per_ptm
        y_0_groups = y_0_line + (label_plot_height + 10) * group_direction
        vertical_space_left = utils.get_height() - y_0_groups if above == 'A' else y_0_groups
        # offset for region label
        dy_label = offset_region_label_from_angle()
        vertical_space_left -= dy_label*2
        dy = vertical_space_left//len(mean_values.index)*group_direction

        plot_group_labels_horizontal(fig, mean_values, y_0_groups, dy, dx)
    else:
        dy = pixels_per_ptm
        x_0_groups = x_0_line + (label_plot_height + 10) * group_direction
        horizontal_space_left = utils.get_width() - x_0_groups if above == 'A' else x_0_groups
        # offset for region label
        dx_label = offset_region_label_from_angle()
        horizontal_space_left -= dx_label*2
        dx = horizontal_space_left//len(mean_values.index)*group_direction

        plot_group_labels_vertical(fig, mean_values, x_0_groups, dx, dy)
    
    previous_ptm = 0
    for i, ptm in enumerate(ptms):
        ptm_position = int(ptm[1:])
        if ptm_position > last_end or ptm_position < previous_ptm:
            if CONFIG.FIGURE_ORIENTATION == 0:
                start_idx = ptm_idx - (i - first_ptm_in_region)
                x_0_groups = start_idx * pixels_per_ptm + utils.SEQUENCE_OFFSET
                x_divider = ptm_idx * pixels_per_ptm + utils.SEQUENCE_OFFSET                  
                x_label = x_0_groups + (x_divider-x_0_groups)//2 - dx//2
                y_label = y_0_groups + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction
                
                plot_groups_horizontal(fig, mean_values.iloc[:,first_ptm_in_region:i], x_0_groups, y_0_groups, dx, dy, x_label, y_label, last_region, group_direction, True)

                fig.add_trace(go.Scatter(x=[x_divider,x_divider],
                            y=[y_0_groups, y_0_groups+len(mean_values.index)*dy],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
            else:
                start_idx = ptm_idx - (i - first_ptm_in_region)
                y_0_groups = utils.get_height() - start_idx * pixels_per_ptm - utils.SEQUENCE_OFFSET
                y_divider = utils.get_height() - ptm_idx * pixels_per_ptm - utils.SEQUENCE_OFFSET
                y_label = y_0_groups - (y_0_groups - y_divider)//2 + dy//2
                x_label = x_0_groups + len(mean_values.index)*dx + (5+utils.get_label_height()//2) * group_direction

                plot_groups_vertical(fig, mean_values.iloc[:,first_ptm_in_region:i], x_0_groups, y_0_groups, dx, dy, x_label, y_label, last_region, group_direction, True)

                fig.add_trace(go.Scatter(x=[x_0_groups, x_0_groups+len(mean_values.index)*dx],
                            y=[y_divider, y_divider],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
            if ptm_position < previous_ptm:
                last_region += 1
                last_end = CONFIG.REGIONS[last_region][1]
            else:
                while ptm_position > last_end:
                    last_region += 1
                    last_end = CONFIG.REGIONS[last_region][1]
            ptm_idx += 1
            first_ptm_in_region = i
        if CONFIG.FIGURE_ORIENTATION == 0:
            position = utils.get_position_with_offset(ptm_position, isoforms[i])
            x_0_line = position * utils.PIXELS_PER_AA + utils.SEQUENCE_OFFSET
            x_0_line = utils.offset_line_for_exon(x_0_line, ptm_position, CONFIG.FIGURE_ORIENTATION)
            x_1_line = ptm_idx * pixels_per_ptm + utils.SEQUENCE_OFFSET
            y_3_line = y_2_line + 10 * group_direction
            if second_row and i % 2 == 1:
                x_1_line = ptm_idx * pixels_per_ptm + utils.SEQUENCE_OFFSET
                y_3_line = y_2_line + (label_length + 10 + 5) * group_direction
            y_label = y_3_line + (utils.get_label_length(ptm)+10) // 2 * group_direction
            text_color = CONFIG.MODIFICATIONS[str(ptm_df.iloc[0,i+2])][1]
            plot_line_with_label_horizontal(fig, x_0_line, x_1_line, y_0_line, y_1_line, y_2_line, y_3_line, y_label, ptm, True, text_color, str(ptm_df.iloc[0,i+2]))
            x_0_rect = x_1_line - dx//2
            fig.add_shape(type='rect',
                      x0 = x_0_rect,
                      x1 = x_0_rect + dx,
                      y0 = y_0_line + (label_plot_height-PLOT_CONFIG.PTM_RECT_LENGTH)*group_direction,
                      y1 = y_0_line + label_plot_height*group_direction,
                      fillcolor=text_color,
                      line=dict(width=1, color='grey'),
                      showlegend=False,)
        else:
            position = utils.get_position_with_offset(ptm_position, isoforms[i])
            y_0_line = utils.get_height() - position * utils.PIXELS_PER_AA - utils.SEQUENCE_OFFSET
            y_0_line = utils.offset_line_for_exon(y_0_line, ptm_position, CONFIG.FIGURE_ORIENTATION)
            y_1_line = utils.get_height() - ptm_idx * pixels_per_ptm - utils.SEQUENCE_OFFSET
            x_3_line = x_2_line + 10 * group_direction
            if second_row and i % 2 == 1:
                y_1_line = utils.get_height() - ptm_idx * pixels_per_ptm - utils.SEQUENCE_OFFSET
                x_3_line = x_2_line + (label_length + 10 + 5) * group_direction
            x_label = x_3_line + (utils.get_label_length(ptm)+10) // 2 * group_direction
            text_color = CONFIG.MODIFICATIONS[str(ptm_df.iloc[0,i+2])][1]
            plot_line_with_label_vertical(fig, x_0_line, x_1_line, x_2_line, x_3_line, y_0_line, y_1_line, x_label, ptm, True, text_color, str(ptm_df.iloc[0,i+2]))
            y_0_rect = y_1_line - dy//2
            fig.add_shape(type='rect',
                      x0 = x_0_line + (label_plot_height-PLOT_CONFIG.PTM_RECT_LENGTH)*group_direction,
                      x1 = x_0_line + label_plot_height*group_direction,
                      y0 = y_0_rect,
                      y1 = y_0_rect + dy,
                      fillcolor=text_color,
                      line=dict(width=1, color='grey'),
                      showlegend=False,)
        ptm_idx += 1
        previous_ptm = ptm_position

    
    while ptm_position > last_end:
        last_region += 1
        last_end = CONFIG.REGIONS[last_region][1]

    if isoforms[first_ptm_in_region] == 'exon2' and isoforms[first_ptm_in_region-1] != 'exon1':
        last_region += 1

    # plot groups for last region
    if CONFIG.FIGURE_ORIENTATION == 0:
        start_idx = ptm_idx - (i - first_ptm_in_region)-1
        x_0_groups = start_idx * pixels_per_ptm + utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_ptm_in_region:].columns)
        x_label = x_0_groups + (region_length * pixels_per_ptm)//2 - dx//2
        y_label = y_0_groups + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction
        plot_groups_horizontal(fig, mean_values.iloc[:,first_ptm_in_region:], x_0_groups, y_0_groups, dx, dy, x_label, y_label, last_region, group_direction, True)
        
        create_custome_colorscale(fig, vertical_space_left, group_direction, x_0_groups, y_0_groups, region_length, pixels_per_ptm, True)
        
    else:
        start_idx = ptm_idx - (i - first_ptm_in_region)-1
        y_0_groups = utils.get_height() - start_idx * pixels_per_ptm - utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_ptm_in_region:].columns)
        y_label = y_0_groups - (region_length * pixels_per_ptm)//2 + dy//2
        x_label = x_0_groups + len(mean_values.index)*dx + (5+utils.get_label_height()//2) * group_direction
        plot_groups_vertical(fig, mean_values.iloc[:,first_ptm_in_region:], x_0_groups, y_0_groups, dx, dy, x_label, y_label, last_region, group_direction, True)

        create_custome_colorscale(fig, horizontal_space_left, group_direction, x_0_groups, y_0_groups, region_length, pixels_per_ptm, True)
        
def create_custome_colorscale(fig: go.Figure, vertical_space_left: int, group_direction: int, x_0_groups: int, y_0_groups: int, region_length: int, pixels_per_step: int, ptm: bool):
    if ptm:
        colorscale = [
            [0.0, PLOT_CONFIG.PTM_SCALE_COLOR_LOW],
            [0.5, PLOT_CONFIG.PTM_SCALE_COLOR_MID],
            [1.0, PLOT_CONFIG.PTM_SCALE_COLOR_HIGH] 
        ]
        label = PLOT_CONFIG.PTM_LEGEND_TITLE
    else:
        colorscale = [
            [0.0, PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_LOW],
            [0.5, PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_MID],
            [1.0, PLOT_CONFIG.CLEAVAGE_SCALE_COLOR_HIGH] 
        ]
        label = PLOT_CONFIG.CLEAVAGE_LEGEND_TITLE
    # Create a heatmap
    z = np.linspace(0, 1, 100).reshape(100, 1)

    if CONFIG.FIGURE_ORIENTATION == 0:
        dx = 15
        dy = 1
        scale_height = dy * 100 + 10 + utils.get_label_height() * label.count('<br>')
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
    fig.add_trace(go.Heatmap(
        x0=x_bar,
        y0=y_bar,
        z=z,
        dx=dx,
        dy=dy,
        colorscale=colorscale,
        showscale=False,
        hoverinfo='none',
    ))
    for i in range(3):
        percentage_label = f'{i*50}%'
        if CONFIG.FIGURE_ORIENTATION == 0:
            x_scale = x_bar + 15 + utils.get_label_length(percentage_label)//2
            y_scale = y_bar + i*100*dy/2
        else:
            x_scale = x_bar + i*100*dx/2
            y_scale = y_bar - utils.get_label_height()
        fig.add_annotation(x=x_scale,
                            y=y_scale,
                            text=percentage_label,
                            showarrow=False,
                            font=dict(
                                family=CONFIG.FONT,
                                size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                                color='black',
                                ))
    longest_label = ''
    for string in label.split('<br>'):
        if utils.get_label_length(string) > utils.get_label_length(longest_label):
            longest_label = string
    if CONFIG.FIGURE_ORIENTATION == 0:
        x_legend_title = x_bar + utils.get_label_length(longest_label)//2 - 15
        y_legend_title = y_bar + scale_height
    else:
        x_legend_title = x_bar + dx * 50
        y_legend_title = y_scale - utils.get_label_height() * (label.count('<br>')+1)
    fig.add_annotation(x=x_legend_title,
                        y=y_legend_title,
                        text=label,
                        showarrow=False,
                        font=dict(
                            family=CONFIG.FONT,
                            size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                            color='black',
                            ))

def filter_relevant_modification_sights(ptm_file: str, threshold: int):
    df = pd.read_csv(ptm_file)
    columns_to_keep = []
    for col in df.columns:
        if CONFIG.INCLUDED_MODIFICATIONS.get(df[col].iloc[0]):
            if df[col].iloc[1][:1] not in CONFIG.INCLUDED_MODIFICATIONS.get(df[col].iloc[0]):
                continue
            if df[col].iloc[0] not in CONFIG.MODIFICATIONS:
                continue
            if df[col].iloc[1][:1] == 'R' and df[col].iloc[0] == 'Deamidation':
                df[col].iloc[0] = 'Citrullination'
            columns_to_keep.append(col)
    df_filtered = df[columns_to_keep]
    df_values = df_filtered.iloc[3:].astype(int)
    sums = df_values.sum()
    filtered_columns = sums[sums >= threshold].index
    # all filter options result in an empty dataframe, check relevant modifications and threshold to keep more columns
    assert len(filtered_columns) > 0
    filtered_df = df[filtered_columns]

    result_df = pd.concat([df.iloc[:, :2], filtered_df], axis=1)
    
    return result_df

def calculate_legend_space(ptm: bool):
    if CONFIG.FIGURE_ORIENTATION == 0:
        longest_label = ''
        if ptm:
            for string in PLOT_CONFIG.PTM_LEGEND_TITLE.split('<br>'):
                if utils.get_label_length(string) > utils.get_label_length(longest_label):
                    longest_label = string
        else:
            for string in PLOT_CONFIG.CLEAVAGE_LEGEND_TITLE.split('<br>'):
                if utils.get_label_length(string) > utils.get_label_length(longest_label):
                    longest_label = string
        if utils.get_label_length('100%') + 10 > utils.get_label_length(longest_label):
            return utils.get_label_length('100%') + 10
        return utils.get_label_length(longest_label)
    else:
        if ptm:
            title_height = utils.get_label_height() * (PLOT_CONFIG.PTM_LEGEND_TITLE.count('<br')+1)
        else:
            title_height = utils.get_label_height() * (PLOT_CONFIG.CLEAVAGE_LEGEND_TITLE.count('<br')+1)
        return utils.get_label_height() + title_height + 10
    
def create_details_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
    legend = None
    if not 'A' in PLOT_CONFIG.INPUT_FILES.keys():
        if PLOT_CONFIG.INPUT_FILES['B'][0] == 'PTM':
            legend = 'B'
        fig = sequence_plot.create_plot(input_file, 'A', legend)
    elif not 'B' in PLOT_CONFIG.INPUT_FILES.keys():
        if PLOT_CONFIG.INPUT_FILES['A'][0] == 'PTM':
            legend = 'A'
        fig = sequence_plot.create_plot(input_file, 'B', legend)
    else:
        if PLOT_CONFIG.INPUT_FILES['A'][0] == 'PTM':
            legend = 'A'
        if PLOT_CONFIG.INPUT_FILES['B'][0] == 'PTM':
            legend = 'B'
        fig = sequence_plot.create_plot(input_file, None, legend)
    cleavage_file_path = None
    ptm_file_path = None
    for above in PLOT_CONFIG.INPUT_FILES.keys():
        match PLOT_CONFIG.INPUT_FILES[above][0]:
            case 'Cleavage':
                cleavage_file_path = PLOT_CONFIG.INPUT_FILES[above][1]
                cleavage_above = above
            case 'PTM':
                ptm_file_path = PLOT_CONFIG.INPUT_FILES[above][1]
                ptm_above = above
    
    if CONFIG.FIGURE_ORIENTATION == 0:
        plot_space = utils.get_width()-utils.SEQUENCE_BOUNDARIES['x0']
    else:
        # first we calculate the missing space above the sequence and then subtract it from the total height
        plot_space = utils.get_height() - (utils.get_height()-utils.SEQUENCE_BOUNDARIES['y0'])

    label_plot_height = 150
    
    if cleavage_file_path:
        cleavage_df = pd.read_csv(cleavage_file_path)
        present_regions = get_present_regions_cleavage(cleavage_df)
        number_of_cleavages = len(cleavage_df.columns)
        number_of_dividers = present_regions.count(True)-1
        cleavage_space = plot_space - calculate_legend_space(False)
        pixels_per_cleavage = cleavage_space // (number_of_cleavages + number_of_dividers)
        assert(pixels_per_cleavage >= CONFIG.FONT_SIZE)

        plot_cleavages(fig, cleavage_df, pixels_per_cleavage, label_plot_height, cleavage_above)

    if ptm_file_path:
        ptm_df = filter_relevant_modification_sights(ptm_file_path, PLOT_CONFIG.MODIFICATION_THRESHOLD)
        present_regions = get_present_regions_ptm(ptm_df)
        number_of_ptms = len(ptm_df.columns)
        number_of_dividers = present_regions.count(True)-1
        second_row = False
        ptm_space = plot_space - calculate_legend_space(True)
        pixels_per_ptm = ptm_space // (number_of_ptms + number_of_dividers)
        if (number_of_ptms + number_of_dividers) * utils.get_label_height() > 2*ptm_space:
            raise ValueError('Too many PTMs to fit in plot')
        elif (number_of_ptms + 2*number_of_dividers) * utils.get_label_height() > ptm_space:
            second_row = True
        
        plot_ptms(fig, ptm_df, pixels_per_ptm, label_plot_height, ptm_above, second_row)
    
    utils.show_plot(fig, output_path)
