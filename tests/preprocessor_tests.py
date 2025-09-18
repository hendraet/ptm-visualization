"""Test the Mascot preprocessor."""

import pandas as pd

from protein_sequencing.data_preprocessing.preprocessor import mascot, ms_fragger, protein_pilot, max_quant


def compare_files(file1, file2):
    """Compare two files."""
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    if set(df1.columns) != set(df2.columns):
        return False

    df2 = df2[df1.columns]
    for i in range(len(df1.columns) - 1):
        col1, col2 = df1.columns[i], df1.columns[i + 1]
        # Check if neighboring columns are swapped due to same index, but different modification
        if df1[col1].equals(df2[col2]) or df1[col2].equals(df2[col1]):
            df2 = df2.rename(columns={col1: "__temp__", col2: col1, "__temp__": col2})

    return df1.equals(df2)


def test_process_mascot_dir():
    """Test the process_mascot_file function."""

    # Process the Mascot files
    config = 'tests.configs.default_config'
    preprocessor_config = 'tests.configs.mascot_config'
    mascot(config, preprocessor_config)

    # Check the processed files
    compare_files("tests/output/result_mascot.csv", "tests/results/expected_result_mascot.csv")


def test_process_ms_fragger_file():
    """Test the process_ms_fragger_file function."""

    # Process the MS Fragger files
    config = 'tests.configs.default_config'
    preprocessor_config = 'tests.configs.ms_fragger_config'
    ms_fragger(config, preprocessor_config)

    # Check the processed files
    compare_files("tests/output/result_ms_fragger_mods.csv", "tests/results/expected_result_ms_fragger_mods.csv")
    compare_files("tests/output/result_ms_fragger_cleavages.csv",
                  "tests/results/expected_result_ms_fragger_cleavages.csv")


def test_process_protein_pilot_dir():
    """Test the process_protein_pilot_file function."""

    # Process the Protein Pilot files
    config = 'tests.configs.default_config'
    preprocessor_config = 'tests.configs.protein_pilot_config'
    protein_pilot(config, preprocessor_config)

    # Check the processed files
    compare_files("tests/output/result_protein_pilot_mods.csv", "tests/results/expected_result_protein_pilot_mods.csv")
    compare_files("tests/output/result_protein_pilot_cleavages.csv",
                  "tests/results/expected_result_protein_pilot_cleavages.csv")


def test_max_quant_file():
    """Test the process_max_quant_file function."""

    # Process the MaxQuant files
    config = 'tests.configs.default_config'
    preprocessor_config = 'tests.configs.max_quant_config'
    max_quant(config, preprocessor_config)

    # Check the processed files
    compare_files("tests/output/result_max_quant_mods.csv", "tests/results/expected_result_max_quant_mods.csv")
    compare_files("tests/output/result_max_quant_cleavages.csv",
                  "tests/results/expected_result_max_quant_cleavages.csv")
