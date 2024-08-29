# Details plot settings
MODIFICATION_THRESHOLD = 10

INPUT_FILES = {
    'B': ('Cleavage', 'data/chris/cleavage_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_C_collmean_tarik.csv'),
    'A': ('PTM', 'data/chris/cleavage_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_tarik.csv'),
}

CLEAVAGES_TO_HIGHLIGHT = ['2-4', '15']
CLEAVAGE_HIGHLIGHT_LABEL_COLOR = '#ff0000'

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
DETAILS_PLOT_PTM_RECT_LENGTH = 25
REGION_LABEL_ANGLE_NEUROPATHOLOGIES = 30