from collections import defaultdict
import math
import os
import plotly.graph_objects as go
import pandas as pd
import numpy as np

from protein_sequencing import parameters, utils, sequence_plot, constants

def get_present_regions(positions):
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
    for _, region_end, _, _ in parameters.REGIONS:
        region_ranges.append((region_start, region_end))
        region_start = region_end + 1

    regions_present = [False] * len(region_ranges)
    region_index = 0
    for range in ranges:
        while range[0] > region_ranges[region_index][1]:
            region_index += 1
        regions_present[region_index] = True
    return regions_present

def get_present_regions_cleavage(cleavage_df: pd.DataFrame):
    cleavages = cleavage_df.iloc[0:1,2:].values[0].tolist()
    return get_present_regions(cleavages)

def get_present_regions_ptm(ptm_df: pd.DataFrame):
    ptms = ptm_df.iloc[1:2,2:].values[0].tolist()
    ptms = [ptm[1:] for ptm in ptms]
    return get_present_regions(ptms)
    

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
        if f'{ptm_modification}({label[0]})@{label[1:]}' in parameters.PTMS_TO_HIGHLIGHT:
            fig.add_shape(type='rect',
                            x0 = x_1-utils.get_label_height()//2-1,
                            x1 = x_1+utils.get_label_height()//2+1,
                            y0 = y_label-utils.get_label_length(label)//2-3,
                            y1 = y_label+utils.get_label_length(label)//2+3,
                            line=dict(width=0),
                            fillcolor=parameters.PTM_HIGHLIGHT_LABEL_COLOR,
                            showlegend=False,)
    else:
        color = parameters.CLEAVAGE_LABEL_COLOR
        if label in parameters.CLEAVAGES_TO_HIGHLIGHT:
            color = parameters.CLEAVAGE_HIGHLIGHT_LABEL_COLOR
    fig.add_annotation(x=x_1, y=y_label,
                        text=label,
                        showarrow=False,
                        textangle=-90,
                        font=dict(
                            family=parameters.FONT,
                            size=parameters.SEQUENCE_PLOT_FONT_SIZE,
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
        if f'{ptm_modification}({label[0]})@{label[1:]}' in parameters.PTMS_TO_HIGHLIGHT:
            fig.add_shape(type='rect',
                            x0 = x_label-utils.get_label_length(label)//2-3,
                            x1 = x_label+utils.get_label_length(label)//2+3,
                            y0 = y_1-utils.get_label_height()//2-1,
                            y1 = y_1+utils.get_label_height()//2+1,
                            line=dict(width=0),
                            fillcolor=parameters.PTM_HIGHLIGHT_LABEL_COLOR,
                            showlegend=False,)
    else:
        color = parameters.CLEAVAGE_LABEL_COLOR
        if label in parameters.CLEAVAGES_TO_HIGHLIGHT:
            color = parameters.CLEAVAGE_HIGHLIGHT_LABEL_COLOR
    fig.add_annotation(x=x_label, y=y_1,
                        text=label,
                        showarrow=False,
                        font=dict(
                            family=parameters.FONT,
                            size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                            color=color,
                            ))
    return fig

def plot_range_with_label_horizontal(fig: go.Figure, x_0_start: int, x_0_end: int, x_1: int, y_0: int, y_1: int, y_2: int, y_3: int, y_label: int, label: str):
    fig.add_trace(go.Scatter(x=[x_0_start, x_0_start, x_1, x_1, x_1, x_0_end, x_0_end],
                            y=[y_0, y_1, y_2, y_3, y_2, y_1, y_0],
                            mode='lines',
                            fill='toself',
                            line=dict(color="black", width=1), showlegend=False, hoverinfo='none'))
    
    color = parameters.CLEAVAGE_LABEL_COLOR
    if label in parameters.CLEAVAGES_TO_HIGHLIGHT:
        color = parameters.CLEAVAGE_HIGHLIGHT_LABEL_COLOR
    fig.add_annotation(x=x_1, y=y_label,
                        text=label,
                        showarrow=False,
                        textangle=-90,
                        font=dict(
                            family=parameters.FONT,
                            size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                            color=color,
                            ))
    return fig

def plot_range_with_label_vertical(fig: go.Figure, x_0: int, x_1: int, x_2: int, x_3: int, y_0_start: int, y_0_end: int, y_1: int, x_label: int, label: str):
    fig.add_trace(go.Scatter(x=[x_0, x_1, x_2, x_3, x_2, x_1, x_0],
                            y=[y_0_start, y_0_start, y_1, y_1, y_1, y_0_end, y_0_end],
                            mode='lines',
                            fill='toself',
                            line=dict(color="black", width=1), showlegend=False, hoverinfo='none'))
    color = parameters.CLEAVAGE_LABEL_COLOR
    if label in parameters.CLEAVAGES_TO_HIGHLIGHT:
        color = parameters.CLEAVAGE_HIGHLIGHT_LABEL_COLOR
    fig.add_annotation(x=x_label, y=y_1,
                        text=label,
                        showarrow=False,
                        font=dict(
                            family=parameters.FONT,
                            size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                            color=color,
                            ))
    return fig

def plot_neuropathologies_horizontal(fig: go.Figure, df: pd.DataFrame, x_0_neuropathologies: int, y_0_neuropathologies: int, dx: int, dy: int, x_label: int, y_label: int, last_region: int, group_dircetion: int, ptm: bool):
    x_margin = 0
    if dx % 2 != 0:
        x_margin = 1
    color_low = parameters.CLEAVAGE_SCALE_COLOR_LOW
    color_mid = parameters.CLEAVAGE_SCALE_COLOR_MID
    color_high = parameters.CLEAVAGE_SCALE_COLOR_HIGH
    if ptm:
        color_low = parameters.PTM_SCALE_COLOR_LOW
        color_mid = parameters.PTM_SCALE_COLOR_MID
        color_high = parameters.PTM_SCALE_COLOR_HIGH
    fig.add_shape(type='rect',
                    x0=x_0_neuropathologies - dx//2 - x_margin,
                    y0=y_0_neuropathologies,
                    x1=x_0_neuropathologies + dx * len(df.iloc[0:1,:].columns) - dx//2,
                    y1=y_0_neuropathologies + dy * len(df.index) + 1,
                    fillcolor='grey',
                    line=dict(color='grey', width=1),
                    showlegend=False,
                    layer='below',)
    fig.add_trace(go.Heatmap(z=df,
                    x0=x_0_neuropathologies,
                    y0=y_0_neuropathologies+dy//2,
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
    fig.add_annotation(x=x_label-utils.get_label_height()*group_dircetion, y=y_label-int((7/offset_region_label_from_angle())*group_dircetion),    
            text=parameters.REGIONS[last_region][3],
            showarrow=False,
            textangle=-parameters.REGION_LABEL_ANGLE_NEUROPATHOLOGIES,
            xanchor=xanchor,
            yanchor=yanchor,
            font=dict(
                family=parameters.FONT,
                size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                color='black',
                ))
    return fig

def plot_neurophatologies_vertical(fig: go.Figure, df: pd.DataFrame, x_0_neuropathologies: int, y_0_neuropathologies: int, dx: int, dy: int, x_label: int, y_label: int, last_region: int, group_dircetion: int, ptm: bool):
    y_margin = 0
    if dy % 2 != 0:
        y_margin = 1
    color_low = parameters.CLEAVAGE_SCALE_COLOR_LOW
    color_mid = parameters.CLEAVAGE_SCALE_COLOR_MID
    color_high = parameters.CLEAVAGE_SCALE_COLOR_HIGH
    if ptm:
        color_low = parameters.PTM_SCALE_COLOR_LOW
        color_mid = parameters.PTM_SCALE_COLOR_MID
        color_high = parameters.PTM_SCALE_COLOR_HIGH
    fig.add_shape(type='rect',
                    x0=x_0_neuropathologies,
                    y0=y_0_neuropathologies + dy//2 + y_margin,
                    x1=x_0_neuropathologies + dx * len(df.index) + 1,
                    y1=y_0_neuropathologies - dy * len(df.iloc[0:1,:].columns) + dy//2,
                    fillcolor='grey',
                    line=dict(color='grey', width=1),
                    showlegend=False,
                    layer='below',)
    
    fig.add_trace(go.Heatmap(z=df.T,
                    x0=x_0_neuropathologies+dx//2,
                    y0=y_0_neuropathologies,
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
            text=parameters.REGIONS[last_region][3],
            showarrow=False,
            textangle=-parameters.REGION_LABEL_ANGLE_NEUROPATHOLOGIES+90,
            xanchor=xanchor,
            font=dict(
                family=parameters.FONT,
                size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                color='black',
                ))
    return fig

def plot_neuropathology_labels_horizontal(fig: go.Figure, mean_values: pd.DataFrame, y_0_neuropathologies: int, dy: int, dx: int):
    for i, neuropathology in enumerate(mean_values.index):
        y_0_rect = y_0_neuropathologies + i*dy
        x_1_rect = utils.SEQUENCE_OFFSET - dx//2
        fig.add_shape(type='rect',
                      x0 = 0,
                      x1 = x_1_rect,
                      y0 = y_0_rect,
                      y1 = y_0_rect + dy,
                      fillcolor=parameters.NEUROPATHOLOGIES[neuropathology][1],
                      line=dict(width=0),
                      showlegend=False,
                      layer='below',)
        color = utils.get_label_color(neuropathology)

        fig.add_annotation(x=x_1_rect//2, y=y_0_rect + dy//2,
            text=neuropathology,
            showarrow=False,
            align='center',
            font=dict(
                family=parameters.FONT,
                size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                color=color))

def plot_neuropathology_labels_vertical(fig: go.Figure, mean_values: pd.DataFrame, x_0_neuropathologies: int, dx: int, dy: int):
    for i, neuropathology in enumerate(mean_values.index):
        x_0_rect = x_0_neuropathologies + i*dx
        y_0_rect = utils.get_height()
        y_rect = utils.SEQUENCE_OFFSET - dy//2
        fig.add_shape(type='rect',
                      x0 = x_0_rect,
                      x1 = x_0_rect + dx,
                      y0 = y_0_rect,
                      y1 = y_0_rect - y_rect,
                      fillcolor=parameters.NEUROPATHOLOGIES[neuropathology][1],
                      line=dict(width=0),
                      showlegend=False,
                      layer='below',)
        
        color = utils.get_label_color(neuropathology)

        fig.add_annotation(x=x_0_rect + dx//2, y=y_0_rect - y_rect//2,
            text=neuropathology,
            showarrow=False,
            align='center',
            textangle=90,
            font=dict(
                family=parameters.FONT,
                size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                color=color))
        
def preprocess_neuropathologies(df: pd.DataFrame, ptm: bool):
    df.columns = df.iloc[0]
    if ptm:
        labels = df.iloc[1:2,2:].values.flatten().tolist()
        df = df.iloc[2:]
    else:
        labels = df.iloc[0:1,2:].values.flatten().tolist()
        df = df.iloc[1:]

    reverse_neuropathology_mapping = {}
    for key, list in parameters.NEUROPATHOLOGIES.items():
        values = list[0]
        for value in values:
            reverse_neuropathology_mapping[value] = key
    df.iloc[:,1] = df.iloc[:,1].map(reverse_neuropathology_mapping)
    mean_values = df.iloc[:,2:].astype(float).groupby(df.iloc[:,1]).mean()
    mean_values = mean_values.reindex([*parameters.NEUROPATHOLOGIES])

    return mean_values, labels

def offset_region_label_from_angle():
    longest_label = ''
    for i, (_, _, _, region_label_short) in enumerate(parameters.REGIONS):
        if utils.get_label_length(region_label_short) > utils.get_label_length(longest_label):
            longest_label = region_label_short
        
    length = utils.get_label_length(longest_label)
    height = utils.get_label_height()

    angle_radians = math.radians(-parameters.REGION_LABEL_ANGLE_NEUROPATHOLOGIES)
    dy = abs((length / 2) * math.sin(angle_radians)) + abs((height / 2) * math.cos(angle_radians))
    
    return int(dy)


def plot_cleavages(fig: go.Figure, cleavage_df: pd.DataFrame, pixels_per_cleavage: int, label_plot_height: int, group: str):
    # prepare data:
    mean_values, cleavages = preprocess_neuropathologies(cleavage_df, False)
    if group == 'B':
        mean_values = mean_values.iloc[::-1]

    longest_label = ''
    for cleavage in cleavages[::-1]:
        if utils.get_label_length(str(cleavage)) > utils.get_label_length(longest_label):
            longest_label = str(cleavage)

    group_direction = 1 if group == 'A' else -1
    first_cleavage_in_region = 0
    cleavages_visited = 0
    last_end = parameters.REGIONS[0][1]
    last_region = 0

    if parameters.FIGURE_ORIENTATION == 0:
        y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
        y_1_line = y_0_line + 10 * group_direction
        y_2_line = y_0_line + (label_plot_height - utils.get_label_length(longest_label) - 10) * group_direction

        y_0_neuropathologies = y_0_line + (label_plot_height + 10) * group_direction
        vertical_space_left = utils.get_height() - y_0_neuropathologies if group == 'A' else y_0_neuropathologies
        # offset for border around heatmap
        vertical_space_left -= 2
        # offset for label for region
        dy_label = offset_region_label_from_angle()
        vertical_space_left -= dy_label*2 + 5
        dx = pixels_per_cleavage
        dy = vertical_space_left//len(mean_values.index)*group_direction

        plot_neuropathology_labels_horizontal(fig, mean_values, y_0_neuropathologies, dy, dx)
    else:
        x_0_line = utils.SEQUENCE_BOUNDARIES['x1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['x0']
        x_1_line = x_0_line + 10 * group_direction
        x_2_line = x_0_line + (label_plot_height - utils.get_label_length(longest_label) - 10) * group_direction

        x_0_neuropathologies = x_0_line + (label_plot_height + 10) * group_direction
        horizontal_space_left = utils.get_width() - x_0_neuropathologies if group == 'A' else x_0_neuropathologies
        # offset for border around heatmap
        horizontal_space_left -= 2
        # offset for label for region
        dx_label = offset_region_label_from_angle()
        horizontal_space_left -= dx_label*2 + 5
        dy = pixels_per_cleavage
        dx = horizontal_space_left//len(mean_values.index)*group_direction

        plot_neuropathology_labels_vertical(fig, mean_values, x_0_neuropathologies, dx, dy)

    for cleavage in cleavages:
        if '-' in str(cleavage):
            start, end = map(int, cleavage.split('-'))
        else:
            start = end = int(cleavage)
        if start > last_end:
            if parameters.FIGURE_ORIENTATION == 0:
                x_0_neuropathologies = first_cleavage_in_region * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                x_divider = cleavages_visited * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                x_label = x_0_neuropathologies + (x_divider-x_0_neuropathologies)//2 - dx//2
                y_label = y_0_neuropathologies + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction

                plot_neuropathologies_horizontal(fig, mean_values.iloc[:,first_cleavage_in_region-last_region:cleavages_visited-last_region], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region, group_direction, False)                
                
                fig.add_trace(go.Scatter(x=[x_divider,x_divider],
                            y=[y_0_neuropathologies, y_0_neuropathologies+len(mean_values.index)*dy],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
            else:
                y_0_neuropathologies = utils.get_height() - first_cleavage_in_region * pixels_per_cleavage - utils.SEQUENCE_OFFSET
                y_divider = utils.get_height() - cleavages_visited * pixels_per_cleavage - utils.SEQUENCE_OFFSET
                y_label = y_0_neuropathologies - (y_0_neuropathologies - y_divider)//2 + dy//2
                x_label = x_0_neuropathologies + len(mean_values.index)*dx + (5+utils.get_label_height()//2) * group_direction

                plot_neurophatologies_vertical(fig, mean_values.iloc[:,first_cleavage_in_region-last_region:cleavages_visited-last_region], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region, group_direction, False)
                
                fig.add_trace(go.Scatter(x=[x_0_neuropathologies, x_0_neuropathologies+len(mean_values.index)*dx],
                            y=[y_divider, y_divider],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
                
            while start > last_end:
                last_region += 1
                last_end = parameters.REGIONS[last_region][1]
            cleavages_visited += 1
            first_cleavage_in_region = cleavages_visited
        if parameters.FIGURE_ORIENTATION == 0:
            if start == end:
                label = str(start)
                x_0_line = start * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
                x_1_line = cleavages_visited * pixels_per_cleavage + utils.SEQUENCE_OFFSET
                y_3_line = y_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction
                
                plot_line_with_label_horizontal(fig,
                                 x_0_line, x_1_line,
                                 y_0_line, y_1_line, y_2_line, y_3_line,
                                 y_label,
                                 label, False, None, None)
            else:
                label = f'{start}-{end}'
                x_0_start_line = start * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
                x_0_end_line = end * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
                x_1_line = cleavages_visited * pixels_per_cleavage + utils.SEQUENCE_OFFSET
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
                y_0_line = utils.get_height() - start * utils.PIXELS_PER_PROTEIN - utils.SEQUENCE_OFFSET
                y_1_line = utils.get_height() - cleavages_visited * pixels_per_cleavage - utils.SEQUENCE_OFFSET
                x_3_line = x_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                x_label = x_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction

                plot_line_with_label_vertical(fig,
                                 x_0_line, x_1_line, x_2_line, x_3_line,
                                 y_0_line, y_1_line,
                                 x_label,
                                 label, False, None, None)
            else:
                label = f'{start}-{end}'
                y_0_start_line = utils.get_height() - start * utils.PIXELS_PER_PROTEIN - utils.SEQUENCE_OFFSET
                y_0_end_line = utils.get_height() - end * utils.PIXELS_PER_PROTEIN - utils.SEQUENCE_OFFSET
                y_1_line = utils.get_height() - cleavages_visited * pixels_per_cleavage - utils.SEQUENCE_OFFSET
                x_3_line = x_0_line + (label_plot_height - utils.get_label_length(label)) * group_direction
                x_label = x_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction

                plot_range_with_label_vertical(fig,
                                    x_0_line, x_1_line, x_2_line, x_3_line,
                                    y_0_start_line, y_0_end_line,
                                    y_1_line,
                                    x_label,
                                    label)
        cleavages_visited += 1
    # plot neuropathologies for last region
    if parameters.FIGURE_ORIENTATION == 0:
        x_0_neuropathologies = first_cleavage_in_region * pixels_per_cleavage + utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_cleavage_in_region-last_region:].columns)
        x_label = x_0_neuropathologies + (region_length * pixels_per_cleavage)//2 - dx//2
        y_label = y_0_neuropathologies + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction
        plot_neuropathologies_horizontal(fig, mean_values.iloc[:,first_cleavage_in_region-last_region:], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region, group_direction, False)

        create_custome_colorscale(fig, vertical_space_left, group_direction, x_0_neuropathologies, y_0_neuropathologies, region_length, pixels_per_cleavage, False)
    else:
        y_0_neuropathologies = utils.get_height() - first_cleavage_in_region * pixels_per_cleavage - utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_cleavage_in_region-last_region:].columns)
        y_label = y_0_neuropathologies - (region_length * pixels_per_cleavage)//2 + dy//2
        x_label = x_0_neuropathologies + len(mean_values.index)*dx + (5+utils.get_label_height()//2) * group_direction
        plot_neurophatologies_vertical(fig, mean_values.iloc[:,first_cleavage_in_region-last_region:], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region, group_direction, False)

        create_custome_colorscale(fig, horizontal_space_left, group_direction, x_0_neuropathologies, y_0_neuropathologies, region_length, pixels_per_cleavage, False)

def plot_ptms(fig: go.Figure, ptm_df: pd.DataFrame, pixels_per_ptm: int, label_plot_height: int, group: str, second_row: bool):
    group_direction = 1 if group == 'A' else -1
    mean_values, ptms = preprocess_neuropathologies(ptm_df, True)
    # For debugging purposes
    #new_row = pd.DataFrame([ptms], columns=mean_values.columns, index=['ptm_position'])
    #mean_values = pd.concat([new_row, mean_values])
    #mean_values.to_csv('plotting_data_ptms.csv', sep=',')

    label_length = utils.get_label_length(ptms[-1])
    # inverse index for group B
    if group == 'B':
        mean_values = mean_values.iloc[::-1]

    if parameters.FIGURE_ORIENTATION == 0:
        y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
        y_1_line = y_0_line + 10 * group_direction
        y_2_line = y_0_line + (label_plot_height - label_length - 10 - constants.DETAILS_PLOT_PTM_RECT_LENGTH - 10) * group_direction
        if second_row:
            y_2_line = y_0_line + (label_plot_height - 2*(label_length + 10) - constants.DETAILS_PLOT_PTM_RECT_LENGTH - 5) * group_direction
    else:
        x_0_line = utils.SEQUENCE_BOUNDARIES['x1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['x0']
        x_1_line = x_0_line + 10 * group_direction
        x_2_line = x_0_line + (label_plot_height - label_length - 10 - constants.DETAILS_PLOT_PTM_RECT_LENGTH - 10) * group_direction
        if second_row:
            x_2_line = x_0_line + (label_plot_height - 2*(label_length + 10) - constants.DETAILS_PLOT_PTM_RECT_LENGTH - 5) * group_direction

    ptms_visited = 0
    ptms_in_region = 0
    last_end = parameters.REGIONS[0][1]
    first_ptm_in_region = 0
    last_region = 0

    if parameters.FIGURE_ORIENTATION == 0:
        dx = pixels_per_ptm
        y_0_neuropathologies = y_0_line + (label_plot_height + 10) * group_direction
        vertical_space_left = utils.get_height() - y_0_neuropathologies if group == 'A' else y_0_neuropathologies
        # offset for region label
        dy_label = offset_region_label_from_angle()
        vertical_space_left -= dy_label*2 + 5
        dy = vertical_space_left//len(mean_values.index)*group_direction

        plot_neuropathology_labels_horizontal(fig, mean_values, y_0_neuropathologies, dy, dx)
    else:
        dy = pixels_per_ptm
        x_0_neuropathologies = x_0_line + (label_plot_height + 10) * group_direction
        horizontal_space_left = utils.get_width() - x_0_neuropathologies if group == 'A' else x_0_neuropathologies
        # offset for region label
        dx_label = offset_region_label_from_angle()
        horizontal_space_left -= dx_label*2 + 5
        dx = horizontal_space_left//len(mean_values.index)*group_direction

        plot_neuropathology_labels_vertical(fig, mean_values, x_0_neuropathologies, dx, dy)
    

    for i, ptm in enumerate(ptms):
        ptm_position = int(ptm[1:])
        if ptm_position > last_end:
            if parameters.FIGURE_ORIENTATION == 0:
                x_0_neuropathologies = first_ptm_in_region * pixels_per_ptm + utils.SEQUENCE_OFFSET
                x_divider = ptms_visited * pixels_per_ptm + utils.SEQUENCE_OFFSET                  
                x_label = x_0_neuropathologies + (x_divider-x_0_neuropathologies)//2 - dx//2
                y_label = y_0_neuropathologies + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction
                
                plot_neuropathologies_horizontal(fig, mean_values.iloc[:,i-ptms_in_region:i], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region, group_direction, True)

                fig.add_trace(go.Scatter(x=[x_divider,x_divider],
                            y=[y_0_neuropathologies, y_0_neuropathologies+len(mean_values.index)*dy],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
            else:
                y_0_neuropathologies = utils.get_height() - first_ptm_in_region * pixels_per_ptm - utils.SEQUENCE_OFFSET
                y_divider = utils.get_height() - ptms_visited * pixels_per_ptm - utils.SEQUENCE_OFFSET
                y_label = y_0_neuropathologies - (y_0_neuropathologies - y_divider)//2 + dy//2
                x_label = x_0_neuropathologies + len(mean_values.index)*dx + (5+utils.get_label_height()//2) * group_direction

                plot_neurophatologies_vertical(fig, mean_values.iloc[:,i-ptms_in_region:i], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region, group_direction, True)

                fig.add_trace(go.Scatter(x=[x_0_neuropathologies, x_0_neuropathologies+len(mean_values.index)*dx],
                            y=[y_divider, y_divider],
                            mode='lines',
                            line=dict(color="black", width=3), showlegend=False, hoverinfo='none'))
            
            while ptm_position > last_end:
                last_region += 1
                last_end = parameters.REGIONS[last_region][1]
            ptms_visited += 1
            first_ptm_in_region = ptms_visited
            ptms_in_region = 0
        ptms_in_region += 1
        if parameters.FIGURE_ORIENTATION == 0:
            x_0_line = ptm_position * utils.PIXELS_PER_PROTEIN + utils.SEQUENCE_OFFSET
            x_1_line = ptms_visited * pixels_per_ptm + utils.SEQUENCE_OFFSET
            y_3_line = y_2_line + 10 * group_direction
            if second_row and i % 2 == 1:
                x_1_line = ptms_visited * pixels_per_ptm + utils.SEQUENCE_OFFSET
                y_3_line = y_2_line + (label_length + 10 + 5) * group_direction
            y_label = y_3_line + (utils.get_label_length(ptm)+10) // 2 * group_direction
            text_color = parameters.MODIFICATIONS[str(ptm_df.iloc[0,i+2])][1]
            plot_line_with_label_horizontal(fig, x_0_line, x_1_line, y_0_line, y_1_line, y_2_line, y_3_line, y_label, ptm, True, text_color, str(ptm_df.iloc[0,i+2]))
            x_0_rect = x_1_line - dx//2
            fig.add_shape(type='rect',
                      x0 = x_0_rect,
                      x1 = x_0_rect + dx,
                      y0 = y_0_line + (label_plot_height-constants.DETAILS_PLOT_PTM_RECT_LENGTH)*group_direction,
                      y1 = y_0_line + label_plot_height*group_direction,
                      fillcolor=text_color,
                      line=dict(width=1, color='grey'),
                      showlegend=False,)
        else:
            y_0_line = utils.get_height() - ptm_position * utils.PIXELS_PER_PROTEIN - utils.SEQUENCE_OFFSET
            y_1_line = utils.get_height() - ptms_visited * pixels_per_ptm - utils.SEQUENCE_OFFSET
            x_3_line = x_2_line + 10 * group_direction
            if second_row and i % 2 == 1:
                y_1_line = utils.get_height() - ptms_visited * pixels_per_ptm - utils.SEQUENCE_OFFSET
                x_3_line = x_2_line + (label_length + 10 + 5) * group_direction
            x_label = x_3_line + (utils.get_label_length(ptm)+10) // 2 * group_direction
            text_color = parameters.MODIFICATIONS[str(ptm_df.iloc[0,i+2])][1]
            plot_line_with_label_vertical(fig, x_0_line, x_1_line, x_2_line, x_3_line, y_0_line, y_1_line, x_label, ptm, True, text_color, str(ptm_df.iloc[0,i+2]))
            y_0_rect = y_1_line - dy//2
            fig.add_shape(type='rect',
                      x0 = x_0_line + (label_plot_height-constants.DETAILS_PLOT_PTM_RECT_LENGTH)*group_direction,
                      x1 = x_0_line + label_plot_height*group_direction,
                      y0 = y_0_rect,
                      y1 = y_0_rect + dy,
                      fillcolor=text_color,
                      line=dict(width=1, color='grey'),
                      showlegend=False,)
        ptms_visited += 1
    # plot neuropathologies for last region
    if parameters.FIGURE_ORIENTATION == 0:
        x_0_neuropathologies = first_ptm_in_region * pixels_per_ptm + utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_ptm_in_region-last_region:].columns)
        x_label = x_0_neuropathologies + (region_length * pixels_per_ptm)//2 - dx//2
        y_label = y_0_neuropathologies + len(mean_values.index)*dy + (5+utils.get_label_height()//2) * group_direction
        plot_neuropathologies_horizontal(fig, mean_values.iloc[:,len(ptms)-ptms_in_region:], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region, group_direction, True)
        
        create_custome_colorscale(fig, vertical_space_left, group_direction, x_0_neuropathologies, y_0_neuropathologies, region_length, pixels_per_ptm, True)
        
    else:
        y_0_neuropathologies = utils.get_height() - first_ptm_in_region * pixels_per_ptm - utils.SEQUENCE_OFFSET
        region_length = len(mean_values.iloc[0:1,first_ptm_in_region-last_region:].columns)
        y_label = y_0_neuropathologies - (region_length * pixels_per_ptm)//2 + dy//2
        x_label = x_0_neuropathologies + len(mean_values.index)*dx + (5+utils.get_label_height()//2) * group_direction
        plot_neurophatologies_vertical(fig, mean_values.iloc[:,len(ptms)-ptms_in_region:], x_0_neuropathologies, y_0_neuropathologies, dx, dy, x_label, y_label, last_region, group_direction, True)

        create_custome_colorscale(fig, horizontal_space_left, group_direction, x_0_neuropathologies, y_0_neuropathologies, region_length, pixels_per_ptm, True)
        
def create_custome_colorscale(fig: go.Figure, vertical_space_left: int, group_direction: int, x_0_neuropathologies: int, y_0_neuropathologies: int, region_length: int, pixels_per_step: int, ptm: bool):
    if ptm:
        colorscale = [
            [0.0, parameters.PTM_SCALE_COLOR_LOW],
            [0.5, parameters.PTM_SCALE_COLOR_MID],
            [1.0, parameters.PTM_SCALE_COLOR_HIGH] 
        ]
        label = parameters.PTM_LEGEND_TITLE
    else:
        colorscale = [
            [0.0, parameters.CLEAVAGE_SCALE_COLOR_LOW],
            [0.5, parameters.CLEAVAGE_SCALE_COLOR_MID],
            [1.0, parameters.CLEAVAGE_SCALE_COLOR_HIGH] 
        ]
        label = parameters.CLEAVAGE_LEGEND_TITLE
    # Create a heatmap
    z = np.linspace(0, 1, 100).reshape(100, 1)

    if parameters.FIGURE_ORIENTATION == 0:
        dx = 15
        dy = 1
        scale_height = dy * 100 + 10 + utils.get_label_height() * label.count('<br>')
        y_offset = (vertical_space_left - scale_height) // 2 * group_direction
        x_bar = x_0_neuropathologies + region_length * pixels_per_step + 10
        y_bar = y_0_neuropathologies + y_offset
        if group_direction == -1:
            y_bar -= scale_height
    else:
        z = z.T
        dx = 1
        dy = 15
        y_scale = dy + 5
        x_offset = vertical_space_left // 2 * group_direction
        y_bar = y_0_neuropathologies - region_length * pixels_per_step - 5
        x_bar = x_0_neuropathologies + x_offset - dx * 50
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
        if parameters.FIGURE_ORIENTATION == 0:
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
                                family=parameters.FONT,
                                size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                                color='black',
                                ))
    if parameters.FIGURE_ORIENTATION == 0:
        x_legend_title = x_bar + 15
        y_legend_title = y_bar + scale_height
    else:
        x_legend_title = x_bar + dx * 50
        y_legend_title = y_scale - utils.get_label_height() * (label.count('<br>')+1)
    fig.add_annotation(x=x_legend_title,
                        y=y_legend_title,
                        text=label,
                        showarrow=False,
                        font=dict(
                            family=parameters.FONT,
                            size=parameters.SEQUENCE_PLOT_FONT_SIZE,
                            color='black',
                            ))

def filter_relevant_modification_sights(ptm_file: str, threshold: int):
    df = pd.read_csv(ptm_file)
    columns_to_keep = [col for col in df.columns if df[col].iloc[0] in parameters.MODIFICATIONS.keys() and
                       (df[col].iloc[1][:1] not in parameters.EXCLUDED_MODIFICATIONS.keys() or
                       df[col].iloc[0] not in parameters.EXCLUDED_MODIFICATIONS[df[col].iloc[1][:1]])]
    df_filtered = df[columns_to_keep]
    df_values = df_filtered.iloc[2:].astype(int)
    sums = df_values.sum()
    filtered_columns = sums[sums >= threshold].index
    # all filter options result in an empty dataframe, check relevant modifications and threshold to keep more columns
    assert len(filtered_columns) > 0
    filtered_df = df[filtered_columns]

    result_df = pd.concat([df.iloc[:, :2], filtered_df], axis=1)
    
    return result_df

def calculate_legend_space(ptm: bool):
    if parameters.FIGURE_ORIENTATION == 0:
        longest_label = ''
        if ptm:
            for string in parameters.PTM_LEGEND_TITLE.split('<br>'):
                if utils.get_label_length(string) > utils.get_label_length(longest_label):
                    longest_label = string
        else:
            for string in parameters.CLEAVAGE_LEGEND_TITLE.split('<br>'):
                if utils.get_label_length(string) > utils.get_label_length(longest_label):
                    longest_label = string
        if utils.get_label_length('100%') + 10 > utils.get_label_length(longest_label):
            return utils.get_label_length('100%') + 10
        return utils.get_label_length(longest_label)
    else:
        if ptm:
            title_height = utils.get_label_height() * (parameters.PTM_LEGEND_TITLE.count('<br')+1)
        else:
            title_height = utils.get_label_height() * (parameters.CLEAVAGE_LEGEND_TITLE.count('<br')+1)
        return utils.get_label_height() + title_height + 10
    
def create_details_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
    legend = None
    if not 'A' in parameters.INPUT_FILES.keys():
        if parameters.INPUT_FILES['B'][0] == 'PTM':
            legend = 'B'
        fig = sequence_plot.create_plot(input_file, 'A', legend)
    elif not 'B' in parameters.INPUT_FILES.keys():
        if parameters.INPUT_FILES['A'][0] == 'PTM':
            legend = 'A'
        fig = sequence_plot.create_plot(input_file, 'B', legend)
    else:
        if parameters.INPUT_FILES['A'][0] == 'PTM':
            legend = 'A'
        if parameters.INPUT_FILES['B'][0] == 'PTM':
            legend = 'B'
        fig = sequence_plot.create_plot(input_file, None, legend)
    cleavage_file_path = None
    ptm_file_path = None
    for group in parameters.INPUT_FILES.keys():
        match parameters.INPUT_FILES[group][0]:
            case 'Cleavage':
                cleavage_file_path = parameters.INPUT_FILES[group][1]
                cleavage_group = group
            case 'PTM':
                ptm_file_path = parameters.INPUT_FILES[group][1]
                ptm_group = group
    
    if parameters.FIGURE_ORIENTATION == 0:
        plot_space = utils.get_width()-utils.SEQUENCE_BOUNDARIES['x0']
    else:
        # first we calculate the missing space above the sequence and then subtract it from the total height
        plot_space = utils.get_height() - (utils.get_height()-utils.SEQUENCE_BOUNDARIES['y0'])

    label_plot_height = 150
    
    if cleavage_file_path:
        cleavage_df = pd.read_csv(cleavage_file_path)
        present_regions = get_present_regions_cleavage(cleavage_df)
        cleavages = cleavage_df.iloc[0:1,2:].values[0].tolist()
        cleavage_space = plot_space - calculate_legend_space(False)
        pixels_per_cleavage = cleavage_space // (len(cleavages) + present_regions.count(True)-1)
        assert(pixels_per_cleavage >= parameters.FONT_SIZE)

        plot_cleavages(fig, cleavage_df, pixels_per_cleavage, label_plot_height, cleavage_group)

    if ptm_file_path:
        ptm_df = filter_relevant_modification_sights(ptm_file_path, parameters.MODIFICATION_THRESHOLD)
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
        
        plot_ptms(fig, ptm_df, pixels_per_ptm, label_plot_height, ptm_group, second_row)
    
    utils.show_plot(fig, output_path)

create_details_plot(parameters.FASTA_INPUT_FILE, parameters.OUTPUT_FOLDER)