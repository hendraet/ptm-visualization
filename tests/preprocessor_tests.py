"""Test the Mascot preprocessor."""

import filecmp
import os

def test_process_mascot_dir():
    """Test the process_mascot_file function."""

    # Process the Mascot files
    os.system("python3 protein_sequencing/data_preprocessing/preprocessor.py -p ma -pc tests.configs.preprocessor_config -c tests.configs.default_config")
    result = filecmp.cmp("tests/output/result_mascot.csv", "tests/results/expected_result_mascot.csv", shallow=False)
    assert result
