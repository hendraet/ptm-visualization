# Plot Settings
FIGURE_ORIENTATION = 0  # 0 for horizontal, 1 for vertical, note figure height and width are then automatically swapped
FIGURE_WIDTH = 1500
FIGURE_HEIGHT = 1000

FONT_SIZE = 12

# Sequence Settings
REGIONS = [
    ('', 44, 'A'),
    ('N1', 73, 'B'),
    ('N2', 102, 'B'),
    ('2N4R-Tau', 150, 'A'),
    ('Proline-rich region', 241, 'B'),
    ('R1', 272, 'B'),
    ('R2', 303, 'B'),
    ('R3', 334, 'B'),
    ('R4', 371, 'B'),
    ('', 441, 'A'),
]

# Modification Settings
MODIFICATIONS = {
    'Phospho': ('Phosphorylation', '#000000', 'A'),
    'Acetyl': ('Acetylation', '#93478F', 'B'),
    'Methyl': ('Methylation', '#C35728', 'B'),
    'GG': ('Ubiquitination', '#7AB77C', 'B'),
    'Citrullination': ('Citrullination', '#FF17E3', 'B'),
}

EXCLUDED_MODIFICATIONS = {'Q': None,
                          'X': None,
                          'S': ['GG']}

# Input Output Settings
FASTA_INPUT_FILE = 'data/uniprot_data/tau_isoforms2N4R.fasta'
OUTPUT_FOLDER = 'output'


# Default Parameters
FONT = 'Arial'

# Margins
LEFT_MARGIN = 0.025
RIGHT_MARGIN = 0.025
TOP_MARGIN = 0.025
BOTTOM_MARGIN = 0.025

# Sequence Plot
SEQUENCE_PLOT_FONT_SIZE = FONT_SIZE
SEQUENCE_PLOT_HEIGHT = 50
# Sequence Region Colors
SEQUENCE_REGION_COLORS = {
    'A': 'white',
    'B': 'lightgrey',
}
# Sequence Minimum Line Length
SEQUENCE_MIN_LINE_LENGTH = 20
