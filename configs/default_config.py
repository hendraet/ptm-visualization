# General Settings

# Sequence Settings
# First sequence is from (1, 44), second from (45, 73) and so on
# Region Name, Region End, Group, Region Abbreviation
# REGIONS = [
#     ('N-term', 44, 'A', 'N-term'),
#     ('N1', 73, 'B', 'N1'),
#     ('N2', 102, 'B', 'N2'),
#     ('2N4R-Tau', 150, 'A', 'Mid'),
#     ('Proline-rich region', 241, 'B', 'PRR'),
#     ('R1', 272, 'B', 'R1'),
#     ('R2', 303, 'B', 'R2'),
#     ('R3', 334, 'B', 'R3'),
#     ('R4', 371, 'B', 'R4'),
#     ('C-term', 441, 'A', 'C-term'),
# ]
REGIONS = [
    ("N-Term", 72, "A", "N"),
    ("1A", 104, "A", "1A"),
    ("", 115, "A", ""),
    ("1B", 214, "A", "1B"),
    ("", 230, "A", ""),
    ("2A", 252, "A", "2A"),
    ("", 256, "A", ""),
    ("2B", 377, "A", "2B"),
    ("", 390, "A", ""),
    ("α", 432, "B", "α"),
    ("ε", 431, "A", "ε"),
]

# Modification Settings
MODIFICATION_LEGEND_TITLE = 'PTMs'
MODIFICATIONS = {
    'Phospho': ('Phosphorylation', '#000000'),
    'Acetyl': ('Acetylation', '#93478F'),
    'Methyl': ('Methylation', '#C35728'),
    'GG': ('Ubiquitination', '#548056'),
    'Citrullination': ('Citrullination', '#FF17E3'),
    'Deamidated': ('Deamidation', '#34AEEB'),
}

INCLUDED_MODIFICATIONS = {'Phospho': ['S', 'T', 'Y'],
                          'Acetyl': ['K'],
                          'Methyl': ['K', 'R'],
                          'GG': ['K'],
                          'Citrullination': ['R'],
                          'Deamidated': ['N', 'Q', 'R'],}

# Input Output Settings
OUTPUT_FOLDER = 'output'

# Plot Settings
FIGURE_ORIENTATION = 0  # 0 for horizontal, 1 for vertical, note figure height and width are then automatically swapped

# just change width and height to change the size of the figure not the orientation
FIGURE_WIDTH = 1200
FIGURE_HEIGHT = 1000

FONT_SIZE = 12

PTMS_TO_HIGHLIGHT = ['Phospho(S)@61', 'Citrullination(R)@242', 'GG(K)@254', 'Acetyl(K)@267', 'Methyl(K)@311']
PTM_HIGHLIGHT_LABEL_COLOR = '#cfcfcf'


# Default Parameters
FONT = 'Arial'

# Margins for sequence Plot
# TODO remove margins and auto calculate based on legend
LEFT_MARGIN = 0.075
RIGHT_MARGIN = 0.025
TOP_MARGIN = 0.065
BOTTOM_MARGIN = 0.025

# Sequence Plot
SEQUENCE_PLOT_FONT_SIZE = FONT_SIZE
SEQUENCE_PLOT_HEIGHT = 50
# works best with an even number
EXONS_GAP = 10
MIN_EXON_LENGTH = 5

# Sequence Region Colors
SEQUENCE_REGION_COLORS = {
    'A': 'white',
    'B': 'lightgrey',
}
