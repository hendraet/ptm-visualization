# Overview Plot Settings
# MODIFICATIONS_GROUP = {
#     'Phospho': 'B',
#     'Acetyl': 'B',
#     'Methyl': 'B',
#     'GG': 'B',
#     'Citrullination': 'B',
#     'Deamidated': 'A',
# }
MODIFICATIONS_GROUP = {
    'Phospho': 'A',
    'Acetyl': 'A',
    'GG': 'A',
    'Citrullination': 'A',
    'Methyl': 'A',
    'Deamidated': 'B',
}
INPUT_FILE = 'tests/output/result_max_quant_mods.csv'

# Sequence Minimum Distance between label and sequence
SEQUENCE_MIN_LINE_LENGTH = 20
