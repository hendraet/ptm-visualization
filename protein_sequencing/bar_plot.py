from collections import defaultdict
import importlib
import os
import plotly.graph_objects as go
import pandas as pd

from protein_sequencing import utils,sequence_plot

CONFIG = importlib.import_module('configs.default_config', 'configs')
PLOT_CONFIG = importlib.import_module('configs.default_bar', 'configs')

# TODO: linie mittig schwerpunk in plot

def get_bar_positions(modification_sights_all_A: dict[int, list[tuple[int, str, str]]], modification_sights_all_B: dict[int, list[tuple[int, str, str]]]):
   positions_A = []
   positions_B = []
   for aa_position in modification_sights_all_A.keys():
      for modification_sight in modification_sights_all_A[aa_position]:
         positions_A.append(aa_position)
   for aa_position in modification_sights_all_B.keys():
      for modification_sight in modification_sights_all_B[aa_position]:
         positions_B.append(aa_position)

   return positions_A, positions_B

def get_bar_plot_width(group_size_A: int, group_size_B: int):
   legend_height = utils.get_label_height() * (max(value.count('<br>') + 1 for value in PLOT_CONFIG.BAR_NEUROPATHOLOGIES.values()))
   if CONFIG.FIGURE_ORIENTATION == 0:
      legend_height += utils.get_label_length('100%')
      bar_width = (utils.get_width() - legend_height) // max(group_size_A, group_size_B)
      bar_plot_width = bar_width * max(group_size_A, group_size_B)
   else:
      legend_height += utils.get_label_height()
      bar_width = (utils.get_height() - legend_height) // max(group_size_A, group_size_B)
      bar_plot_width = bar_width * max(group_size_A, group_size_B)
   return bar_plot_width


def add_bar_plot(fig: go.Figure, group: str, modification_sights_all: dict[int, list[tuple[int, str, str, str]]], modification_sights_relevant: dict[int, list[tuple[int, str, str, str]]], df, group_positions: list, bar_plot_width: int, label_plot_height: int) -> go.Figure:
   group_direction = 1 if group == 'A' else -1
   bar_width = bar_plot_width // len(group_positions)
   assert bar_width >= CONFIG.FONT_SIZE, f"Too many bars to plot! Bar width: {bar_width} < FONT_SIZE: {CONFIG.FONT_SIZE}."
   bar_plot_margin = CONFIG.FONT_SIZE

   height_offset = 0
   modifications_visited = 0
   bar_percentages = {neuropathology: [] for neuropathology in PLOT_CONFIG.BAR_NEUROPATHOLOGIES.keys()}

   for aa_position in sorted(modification_sights_all.keys(), reverse=True):
      for modification_sight in modification_sights_all[aa_position]:
         if modification_sight not in modification_sights_relevant[aa_position]:
            for neuropathology in PLOT_CONFIG.BAR_NEUROPATHOLOGIES.keys():
               bar_percentages[neuropathology].append(0)
            modifications_visited += 1
            continue
         label, modification_type, _, isoform = modification_sight
         if CONFIG.FIGURE_ORIENTATION == 0:
            # x position for protein sequence
            position = utils.get_position_with_offset(aa_position, isoform)
            x_0_line = position * utils.PIXELS_PER_AA + utils.SEQUENCE_OFFSET
            x_0_line = utils.offset_line_for_exon(x_0_line, aa_position, CONFIG.FIGURE_ORIENTATION)
            # x position for bar plot
            x_1_line = utils.get_width() - (modifications_visited * bar_width + bar_width//2)
            y_0_line = utils.SEQUENCE_BOUNDARIES['y1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['y0']
            y_1_line = y_0_line + group_direction * height_offset
            y_3_line = y_0_line + label_plot_height * group_direction - (utils.get_label_length(label) + 10) * group_direction
            y_2_line = y_3_line - 10 * group_direction
            y_label = y_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction
            y_bar = y_3_line + (utils.get_label_length(label) + 5) * group_direction

            # plot line with label
            plot_line_with_label_horizontal(fig,
                                 x_0_line, x_1_line,
                                 y_0_line, y_1_line, y_2_line, y_3_line,
                                 y_label,
                                 CONFIG.MODIFICATIONS[modification_type][1],
                                 label, modification_type)
            
            space_above_sequence = utils.get_height()-y_0_line if group == 'A' else y_0_line
            space_per_neuropathalogy = (space_above_sequence - label_plot_height)  // (len(df["Neuropathology"].unique())-1)
            max_bar_height = space_per_neuropathalogy - 2*bar_plot_margin
            # plot bars
            for i, neuropathology in enumerate(PLOT_CONFIG.BAR_NEUROPATHOLOGIES.keys()):
               modification_column = [col for col in df.columns if col not in ['ID', 'Neuropathology'] and df[col].iloc[1] == label]
               bar_df = df[['ID', 'Neuropathology'] + modification_column]
               bar_df = bar_df[bar_df['Neuropathology'] == neuropathology]
               values = bar_df.iloc[:, 2].astype(int)
               percentage = values.mean()
               height = max_bar_height * percentage
               bar_percentages[neuropathology].append(percentage)
               x0 = x_1_line - bar_width // 2 * PLOT_CONFIG.BAR_WIDTH
               x1 = x_1_line + bar_width // 2 * PLOT_CONFIG.BAR_WIDTH
               y0 = y_bar + (i * space_per_neuropathalogy + 5) * group_direction
               y1 = y0 + height * group_direction
               if PLOT_CONFIG.INVERT_AXIS_GROUP_B and group == 'B':
                  y0 = y_bar + (i * space_per_neuropathalogy + 5 + max_bar_height) * group_direction
                  y1 = y0 - height * group_direction

               fig.add_shape(
                  type="rect",
                  x0=x0,
                  y0=y0,
                  x1=x1,
                  y1=y1,
                  line=dict(color="black", width=1),
                  fillcolor=CONFIG.MODIFICATIONS[modification_type][1]
               )

         else:
            position = utils.get_position_with_offset(aa_position, isoform)
            y_0_line = utils.get_height() - (position * utils.PIXELS_PER_AA + utils.SEQUENCE_OFFSET)
            y_0_line = utils.offset_line_for_exon(y_0_line, aa_position, CONFIG.FIGURE_ORIENTATION)
            y_1_line = modifications_visited * bar_width + bar_width//2
            x_0_line = utils.SEQUENCE_BOUNDARIES['x1'] if group == 'A' else utils.SEQUENCE_BOUNDARIES['x0']
            x_1_line = x_0_line + group_direction * height_offset
            x_3_line = x_0_line + label_plot_height * group_direction - (utils.get_label_length(label) + 10) * group_direction
            x_2_line = x_3_line - 10 * group_direction
            x_label = x_3_line + (utils.get_label_length(label) // 2 + 5) * group_direction
            x_bar = x_3_line + (utils.get_label_length(label) + 5) * group_direction

            # plot line with label
            plot_line_with_label_vertical(fig,
                                          x_0_line, x_1_line, x_2_line, x_3_line, x_label,
                                          y_0_line, y_1_line,
                                          CONFIG.MODIFICATIONS[modification_type][1],
                                          label, modification_type)
            
            space_above_sequence = utils.get_width()-x_0_line if group == 'A' else x_0_line
            space_per_neuropathalogy = (space_above_sequence - label_plot_height)  // (len(df["Neuropathology"].unique())-1)
            max_bar_height = space_per_neuropathalogy - 2*bar_plot_margin
            # plot bars
            for i, neuropathology in enumerate(PLOT_CONFIG.BAR_NEUROPATHOLOGIES.keys()):
               modification_column = [col for col in df.columns if col not in ['ID', 'Neuropathology'] and df[col].iloc[1] == label]
               bar_df = df[['ID', 'Neuropathology'] + modification_column]
               bar_df = bar_df[bar_df['Neuropathology'] == neuropathology]
               values = bar_df.iloc[:, 2].astype(int)
               percentage = values.mean()
               height = max_bar_height * percentage
               bar_percentages[neuropathology].append(percentage)
               y0 = y_1_line - bar_width // 2 * PLOT_CONFIG.BAR_WIDTH
               y1 = y_1_line + bar_width // 2 * PLOT_CONFIG.BAR_WIDTH
               x0 = x_bar + (i * space_per_neuropathalogy + 5) * group_direction
               x1 = x0 + height * group_direction
               if PLOT_CONFIG.INVERT_AXIS_GROUP_B and group == 'B':
                  x0 = x_bar + (i * space_per_neuropathalogy + 5 + max_bar_height) * group_direction
                  x1 = x0 - height * group_direction

               fig.add_shape(
                  type="rect",
                  x0=x0,
                  y0=y0,
                  x1=x1,
                  y1=y1,
                  line=dict(color="black", width=1),
                  fillcolor=CONFIG.MODIFICATIONS[modification_type][1]
               )
                                 
         modifications_visited += 1

      if modifications_visited > len(group_positions)/2:
         height_offset -= 2
      else:
         height_offset += 2
   
   max_lines = max(value.count('<br>')+1 for value in PLOT_CONFIG.BAR_NEUROPATHOLOGIES.values())
   if CONFIG.FIGURE_ORIENTATION == 0:
      for i, neuropathology in enumerate(PLOT_CONFIG.BAR_NEUROPATHOLOGIES.keys()):
         y_group = y_bar + (i * space_per_neuropathalogy + 5) * group_direction
         x_line_start = utils.get_width()-bar_plot_width
         fig.add_annotation(x=x_line_start-utils.get_label_length('100%')-max_lines*utils.get_label_height(),
                            y = y_group + space_per_neuropathalogy//2 * group_direction,
                            width=space_per_neuropathalogy,
                            text=PLOT_CONFIG.BAR_NEUROPATHOLOGIES[neuropathology],
                            textangle=-90,
                            font=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE, family=CONFIG.FONT, color="black"),
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
                                     line=dict(color='lightgray', width=1),
                                     showlegend=False,
                                     hoverinfo='none'))
            if j % 2 == 0:
               text = str(j*25) + '%'
               if PLOT_CONFIG.INVERT_AXIS_GROUP_B and group == 'B':
                  text = str(100-j*25) + '%'
               fig.add_annotation(x=x_line_start-utils.get_label_length(text)//2-5,
                                    y=y_trace,
                                    text=text,
                                    font=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE, family=CONFIG.FONT, color="black"),
                                    showarrow=False)
   else:
      for i, neuropathology in enumerate(PLOT_CONFIG.BAR_NEUROPATHOLOGIES.keys()):
         x_group = x_bar + (i * space_per_neuropathalogy + 5) * group_direction
         y_line_start = bar_plot_width
         fig.add_annotation(x=x_group + space_per_neuropathalogy//2 * group_direction,
                            y=y_line_start+utils.get_label_height()+max_lines*utils.get_label_height(),
                            text=PLOT_CONFIG.BAR_NEUROPATHOLOGIES[neuropathology],
                            font=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE, family=CONFIG.FONT, color="black"),
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
                                     line=dict(color='lightgray', width=1),
                                     showlegend=False,
                                     hoverinfo='none'))
            if j % 2 == 0:
               text = str(j*25) + '%'
               if PLOT_CONFIG.INVERT_AXIS_GROUP_B and group == 'B':
                  text = str(100-j*25) + '%'
               if j == 0:
                  x_pos = x_trace - utils.get_label_length(text)//2 if group == 'B' else x_trace + utils.get_label_length(text)//2
               elif j == 4:
                  x_pos = x_trace + utils.get_label_length(text)//2 if group == 'B' else x_trace - utils.get_label_length(text)//2
               else:
                  x_pos = x_trace
               fig.add_annotation(x=x_pos,
                                    y=y_line_start+5,
                                    text=text,
                                    font=dict(size=CONFIG.SEQUENCE_PLOT_FONT_SIZE, family=CONFIG.FONT, color="black"),
                                    showarrow=False)

   return fig

def plot_line_with_label_horizontal(fig, x_0, x_1, y_0, y_1, y_2, y_3, y_label, color, label, modification_type):
   fig.add_trace(go.Scatter(x=[x_0, x_0, x_1, x_1],
                            y=[y_0, y_1, y_2, y_3],
                            mode='lines',
                            line=dict(color=color, width=1), showlegend=False, hoverinfo='none'))
   fig.add_annotation(x=x_1, y=y_label,
                      text=label,
                      showarrow=False,
                      textangle=-90,
                      font=dict(
                         family=CONFIG.FONT,
                         size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                         color=color,
                         ))
   if f'{modification_type}({label[0]})@{label[1:]}' in CONFIG.PTMS_TO_HIGHLIGHT:
      fig.add_shape(type='rect',
                        x0 = x_1-utils.get_label_height()//2-1,
                        x1 = x_1+utils.get_label_height()//2+1,
                        y0 = y_label-utils.get_label_length(label)//2-3,
                        y1 = y_label+utils.get_label_length(label)//2+3,
                        line=dict(width=0),
                        fillcolor=CONFIG.PTM_HIGHLIGHT_LABEL_COLOR,
                        showlegend=False)
   return fig

def plot_line_with_label_vertical(fig, x_0, x_1, x_2, x_3, x_label, y_0, y_1, color, label, modification_type):
   fig.add_trace(go.Scatter(x=[x_0, x_1, x_2, x_3],
                            y=[y_0, y_0, y_1, y_1],
                            mode='lines',
                            line=dict(color=color, width=1), showlegend=False, hoverinfo='none'))
   fig.add_annotation(x=x_label, y=y_1,
                      text=label,
                      showarrow=False,
                      font=dict(
                         family=CONFIG.FONT,
                         size=CONFIG.SEQUENCE_PLOT_FONT_SIZE,
                         color=color,
                         ))
   if f'{modification_type}({label[0]})@{label[1:]}' in CONFIG.PTMS_TO_HIGHLIGHT:
      fig.add_shape(type='rect',
                            x0 = x_label-utils.get_label_length(label)//2-3,
                            x1 = x_label+utils.get_label_length(label)//2+3,
                            y0 = y_1-utils.get_label_height()//2-1,
                            y1 = y_1+utils.get_label_height()//2+1,
                            line=dict(width=0),
                            fillcolor=CONFIG.PTM_HIGHLIGHT_LABEL_COLOR,
                            showlegend=False)
   return fig

def filter_relevant_modification_sights(helper_file: str):
   # read csv into df
   df = pd.read_csv(helper_file)

   # only keep first two columns and columns that are in MODIFICATIONS
   columns_to_keep = list(CONFIG.INCLUDED_MODIFICATIONS.keys())
   df = df[[col for col in df.columns if df[col][0] in columns_to_keep or col in ['ID', 'Neuropathology']]]

   # create dict for all modification sights
   all_modification_sights = defaultdict(list)
   for column in df.columns:
      if column in ['ID', 'Neuropathology']:
         continue
      column_modification = df[column][0]
      column_label = df[column][1][0]
      if CONFIG.INCLUDED_MODIFICATIONS.get(column_modification):
         if column_label not in CONFIG.INCLUDED_MODIFICATIONS[column_modification]:
            continue
         if column_label == 'R' and column_modification == 'Deamidated':
            column_modification = 'Citrullination'
      isoform = df[column][2]
      all_modification_sights[int(df[column][1][1:])].append((df[column][1], column_modification, PLOT_CONFIG.MODIFICATIONS_GROUP[column_modification], isoform))

   # remove rows where Neuropathology is not in PLOT_CONFIG.BAR_Neurpathologies
   header_rows = df.iloc[:4, :]
   df = df[df['Neuropathology'].isin(PLOT_CONFIG.BAR_NEUROPATHOLOGIES.keys())]
   df = pd.concat([header_rows, df], ignore_index=True)

   # filter out columns that have no modifications observed for relevant neuropathologies
   data_to_filter = df.iloc[4:, 2:]
   header_columns = df.iloc[:, :2]
   for col in data_to_filter.columns:
      data_to_filter[col] = data_to_filter[col].astype(int)
   df = df[[col for col in data_to_filter if data_to_filter[col].sum(axis=0) > 0]]
   df = pd.concat([header_columns, df], axis=1)

   # create new dict for modification sights
   relevant_modification_sights = defaultdict(list)
   for column in df.columns:
      if column in ['ID', 'Neuropathology']:
         continue
      column_modification = df[column][0]
      column_label = df[column][1][0]
      if CONFIG.INCLUDED_MODIFICATIONS.get(column_modification):
         if column_label not in CONFIG.INCLUDED_MODIFICATIONS[column_modification]:
            continue
         if column_label == 'R' and column_modification == 'Deamidated':
            column_modification = 'Citrullination'
      isoform = df[column][2]
      relevant_modification_sights[int(df[column][1][1:])].append((df[column][1], column_modification, PLOT_CONFIG.MODIFICATIONS_GROUP[column_modification], isoform))
   return all_modification_sights, relevant_modification_sights, df


def create_bar_plot(input_file: str | os.PathLike, output_path: str | os.PathLike):
   sequence_plot.create_plot(input_file, None, 'A')
   all_positions, relevant_positions, df = filter_relevant_modification_sights(PLOT_CONFIG.BAR_INPUT_FILE)
   
   group_a_all, group_b_all = utils.separate_by_group(all_positions)
   group_a_relevant, group_b_relevant = utils.separate_by_group(relevant_positions)

   if len(group_a_relevant) == 0:
      fig = sequence_plot.create_plot(input_file, 'A', 'A')
   elif len(group_b_relevant) == 0:
      fig = sequence_plot.create_plot(input_file, 'B', 'B')
   else:
      fig = sequence_plot.create_plot(input_file, None, 'A')

   positions_A, positions_B = get_bar_positions(group_a_all, group_b_all)
   group_size_A = len(positions_A)
   group_size_B = len(positions_B)
   bar_plot_width = get_bar_plot_width(group_size_A, group_size_B)
   highest_position = positions_A[-1] if len(positions_A) > len(positions_B) else positions_B[-1]
   label_plot_height = (max(group_size_A, group_size_B) + utils.get_label_length(f'X{highest_position}') + 30)

   for (group_label, group_all, group_relevant, group_positions) in [('A', group_a_all, group_a_relevant, positions_A), ('B', group_b_all, group_b_relevant, positions_B)]:
      if len(group_positions) == 0:
         continue
      fig = add_bar_plot(fig, group_label, group_all, group_relevant, df, group_positions, bar_plot_width, label_plot_height)
   

   utils.show_plot(fig, output_path)