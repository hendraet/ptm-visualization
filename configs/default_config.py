# Plot Settings
FIGURE_ORIENTATION = 0  # 0 for horizontal, 1 for vertical, note figure height and width are then automatically swapped
FIGURE_WIDTH = 1200
FIGURE_HEIGHT = 1000

FONT_SIZE = 12

# Sequence Settings
# First sequence is from (1, 44), second from (45, 73) and so on
# Region Name, Region End, Group, Region Abbreviation
REGIONS = [
    ('N-term', 44, 'A', 'N-term'),
    ('N1', 73, 'B', 'N1'),
    ('N2', 102, 'B', 'N2'),
    ('2N4R-Tau', 150, 'A', 'Mid'),
    ('Proline-rich region', 241, 'B', 'PRR'),
    ('R1', 272, 'B', 'R1'),
    ('R2', 303, 'B', 'R2'),
    ('R3', 334, 'B', 'R3'),
    ('R4', 371, 'B', 'R4'),
    ('C-term', 441, 'A', 'C-term'),
]

# Modification Settings
MODIFICATION_LEGEND_TITLE = 'PTMs'
MODIFICATIONS = {
    'Phospho': ('Phosphorylation', '#000000'),
    'Acetyl': ('Acetylation', '#93478F'),
    'Methyl': ('Methylation', '#C35728'),
    'GG': ('Ubiquitination', '#548056'),
    'Citrullination': ('Citrullination', '#FF17E3'),
}

PTMS_TO_HIGHLIGHT = ['Phospho(S)@61', 'Citrullination(R)@242', 'GG(K)@254', 'Acetyl(K)@267', 'Methyl(K)@311']
PTM_HIGHLIGHT_LABEL_COLOR = '#cfcfcf'

EXCLUDED_MODIFICATIONS = {'Q': None,
                          'X': None,
                          'S': ['GG'],}

# Input Output Settings
FASTA_INPUT_FILE = 'data/uniprot_data/tau_isoforms2N4R.fasta'
OUTPUT_FOLDER = 'output'


# Default Parameters
FONT = 'Arial'

# Margins for sequence Plot
# TODO remove margins and auto calculate based on legend
LEFT_MARGIN = 0.065
RIGHT_MARGIN = 0.025
TOP_MARGIN = 0.065
BOTTOM_MARGIN = 0.025

# Sequence Plot
SEQUENCE_PLOT_FONT_SIZE = FONT_SIZE
SEQUENCE_PLOT_HEIGHT = 50
# Sequence Region Colors
SEQUENCE_REGION_COLORS = {
    'A': 'white',
    'B': 'lightgrey',
}