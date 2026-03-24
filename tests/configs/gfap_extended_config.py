from configs.default_config import *

PTMS_TO_HIGHLIGHT = []  # TODO: properly use?
REGIONS = [
    ("N-Term", 72, "A", "N"),
    ("1Aaaa", 104, "A", "1A"),
    ("", 115, "A", ""),
    ("1Bbbb", 214, "A", "1B"),
    ("", 230, "A", ""),
    ("2Aaaa", 252, "A", "2A"),
    ("", 256, "A", ""),
    ("2Bbbb", 377, "A", "2B"),
    ("Pre-Exon", 390, "A", "P"),
    ("alpha", 432, "B", "α"),
    ("epsilon", 431, "A", "ε"),
]

INCLUDED_MODIFICATIONS = {
    "Phospho": ["S", "T", "Y"],
    "Acetyl": ["K"],
    "Methyl": ["K", "R"],
    "GG": ["K", "Q", "I", "R"],
    "Citrullination": ["R", "H", "G", "M", "E"],
    "Deamidated": ["N", "Q", "R"],
}
