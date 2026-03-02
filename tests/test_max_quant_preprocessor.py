from pathlib import Path

import importlib
import pandas as pd
import pytest

from protein_sequencing.data_preprocessing.max_quant_preprocessor import (
    MaxQuantPreprocessor,
)

# TODO: also need files for SATB1 and Tau

GFAP_DATA_DIR = Path("tests/test_data/GFAP_P14136")
GFAP_RESULTS_DIR = Path("tests/results/GFAP_P14136")


def get_gfap_evidence_df():
    return pd.read_csv(GFAP_DATA_DIR / "evidence_df.csv", sep="\t")


def get_gfap_metadata_df():
    return pd.read_csv(GFAP_DATA_DIR / "metadata_df.csv")


@pytest.mark.parametrize(
    "evidence_df, metadata_df, metadata_column, expected_mods_file, expected_cleavages_file",
    [
        (
            get_gfap_evidence_df(),
            None,
            None,
            GFAP_RESULTS_DIR / "expected_modifications.csv",
            GFAP_RESULTS_DIR / "expected_cleavages.csv",
        ),
        (
            get_gfap_evidence_df(),
            get_gfap_metadata_df(),
            "Group",
            GFAP_RESULTS_DIR / "expected_modifications_with_groups.csv",
            GFAP_RESULTS_DIR / "expected_cleavages_with_groups.csv",
        ),
    ],
)
def test_max_quant_file(
    evidence_df,
    metadata_df,
    metadata_column,
    expected_mods_file,
    expected_cleavages_file,
):
    config = "tests.configs.gfap_config"
    preprocessor_config = "tests.configs.max_quant_config"

    MaxQuantPreprocessor(
        importlib.import_module(config, "configs"),
        importlib.import_module(preprocessor_config, "configs"),
        evidence_df=evidence_df,
        metadata_df=metadata_df,
        metadata_column=metadata_column,
    )

    resulting_mods = pd.read_csv("tests/output/result_max_quant_mods.csv")
    resulting_mods = resulting_mods[sorted(resulting_mods.columns)]
    expected_mods = pd.read_csv(expected_mods_file)
    expected_mods = expected_mods[sorted(expected_mods.columns)]
    pd.testing.assert_frame_equal(
        resulting_mods,
        expected_mods,
        check_dtype=False,
    )

    resulting_cleavages = pd.read_csv("tests/output/result_max_quant_cleavages.csv")
    resulting_cleavages = resulting_cleavages[sorted(resulting_cleavages.columns)]
    expected_cleavages = pd.read_csv(expected_cleavages_file)
    expected_cleavages = expected_cleavages[sorted(expected_cleavages.columns)]
    pd.testing.assert_frame_equal(
        resulting_cleavages,
        expected_cleavages,
        check_dtype=False,
    )

# TODO: needs more tests to check for edge cases
