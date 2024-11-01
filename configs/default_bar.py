# Bar Plot Settings
# group from groups.csv and display name
BAR_NEUROPATHOLOGIES = {'FTLD-PiD': 'group 1',
                        'FTLD-Tau': 'group 2',
                        'CTRL': 'group 3',
                        }
BAR_WIDTH = 0.8
INVERT_AXIS_GROUP_B = True
BAR_INPUT_FILE = 'output/result_protein_pilot_mods.csv'
MODIFICATIONS_GROUP = {
    'Phospho': 'A',
    'Acetyl': 'A',
    'GG': 'A',
    'Citrullination': 'A',
    'Methyl': 'B',
}