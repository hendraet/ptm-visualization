# Plot Settings
FIGURE_ORIENTATION = 0  # 0 for horizontal, 1 for vertical, note figure height and width are then automatically swapped
FIGURE_WIDTH = 1200
FIGURE_HEIGHT = 1000

FONT_SIZE = 11

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
# TODO implement for all plots
MODIFICATION_LEGEND_TITLE = 'PTMs'
MODIFICATIONS = {
    'Phospho': ('Phosphorylation', '#000000'),
    'Acetyl': ('Acetylation', '#93478F'),
    'Methyl': ('Methylation', '#C35728'),
    'GG': ('Ubiquitination', '#548056'),
    'Citrullination': ('Citrullination', '#FF17E3'),
}

# TODO implement for all plots
PTMS_TO_HIGHLIGHT = ['Phospho(S)@68', 'Citrullination(R)@242', 'GG(K)@254', 'Acetyl(K)@267', 'Methyl(K)@311']
PTM_HIGHLIGHT_LABEL_COLOR = '#cfcfcf'
CLEAVAGES_TO_HIGHLIGHT = ['2-4', '15']
CLEAVAGE_HIGHLIGHT_LABEL_COLOR = '#ff0000'

EXCLUDED_MODIFICATIONS = {'Q': None,
                          'X': None,
                          'S': ['GG'],}

MODIFICATION_THRESHOLD = 10

INPUT_FILES = {
    'B': ('Cleavage', 'data/chris/cleavage_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_C_collmean_tarik.csv'),
    'A': ('PTM', 'data/chris/cleavage_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_tarik.csv'),
}

CLEAVAGE_LABEL_COLOR = '#333333'
CLEAVAGE_SCALE_COLOR_LOW = '#B35806'
CLEAVAGE_SCALE_COLOR_MID = '#F7F7F7'
CLEAVAGE_SCALE_COLOR_HIGH = '#542788'
CLEAVAGE_LEGEND_TITLE = 'Proteolytic<br>Cleavage<br>Patient<br>Frequency'
#CLEAVAGE_LEGEND_TITLE = 'Proteolytic Cleavage Patient Frequency'

PTM_SCALE_COLOR_LOW = '#B35806'
PTM_SCALE_COLOR_MID = '#F5F5F5'
PTM_SCALE_COLOR_HIGH = '#01665E'
PTM_LEGEND_TITLE = 'PTM Patient <br>Frequency'

NEUROPATHOLOGIES = {"CTR": (["CTR"], '#4DAF4A'),
                    "DLB": (["DLB"], '#8DD3C7'),
                    "PSP": (["PSP"], '#FF7F00'),
                    "PiD": (["PiD"], '#984EA3'),
                    "FTLD": (["FTLDTau"], '#17DFFF'),
                    "CBD": (["CBD"], '#FFD92F'),
                    "CTE": (["CTE"], '#1740B6'),
                    "AD": (["AD", "NPCAD"], '#E41A1C'),
                    "fAD": (["fAD"], '#9C0B0C'),}

REGION_LABEL_ANGLE_NEUROPATHOLOGIES = 30

# Input Output Settings
FASTA_INPUT_FILE = 'data/uniprot_data/tau_isoforms2N4R.fasta'
OUTPUT_FOLDER = 'output'


# Default Parameters
FONT = 'Arial'

# Margins for sequence Plot
# TODO remove margins and auto calculate based on legend
LEFT_MARGIN = 0.065
RIGHT_MARGIN = 0.025
TOP_MARGIN = 0.069
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
