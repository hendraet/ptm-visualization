FONT_SIZE = 12
FIGURE_WIDTH = 1500
FIGURE_HEIGHT = 400
FIGURE_ORIENTATION = 0  # 0 for horizontal, 1 for vertical

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

MODIFICATIONS = {
    'Acetyl': ('Acetylation', 'purple'),
    'Phospho': ('Phosphorylation', 'black'),
    'Methyl': ('Methylation', 'brown'),
    'GG': ('Ubiquitination', 'green'),
    'Deamidated': ('Citrullination', 'pink'),
}

# Default Parameters

# Margins
LEFT_MARGIN = 0.2
RIGHT_MARGIN = 0.1
TOP_MARGIN = 0.3
BOTTOM_MARGIN = 0.1

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
