from pathlib import Path

import importlib
import pandas as pd
import pytest

from protein_sequencing.data_preprocessing.max_quant_preprocessor import (
    MaxQuantPreprocessor,
)
from tests.paths import (
    GFAP_DATA_DIR,
    TAU_DATA_DIR,
    SATB1_DATA_DIR,
    GFAP_EXTENDED_DATA_DIR,
)


def get_evidence_df(data_dir: Path):
    return pd.read_csv(data_dir / "evidence_df.csv", sep="\t")


def get_metadata_df(data_dir):
    return pd.read_csv(data_dir / "metadata_df.csv")


@pytest.fixture
def gfap_evidence_df():
    return get_evidence_df(GFAP_DATA_DIR)


@pytest.fixture
def gfap_extended_evidence_df():
    return get_evidence_df(GFAP_EXTENDED_DATA_DIR)


@pytest.fixture
def tau_evidence_df():
    return get_evidence_df(TAU_DATA_DIR)


@pytest.fixture
def gfap_metadata_df():
    return get_metadata_df(GFAP_DATA_DIR)


@pytest.fixture
def gfap_config():
    return importlib.import_module("tests.configs.gfap_config", "configs")


@pytest.fixture
def gfap_extended_config():
    return importlib.import_module("tests.configs.gfap_extended_config", "configs")


@pytest.fixture
def gfap_preprocessor_config():
    return importlib.import_module(
        "tests.configs.gfap_mq_preprocessor_config", "configs"
    )


@pytest.fixture
def tau_config():
    return importlib.import_module("tests.configs.tau_config", "configs")


@pytest.fixture
def tau_preprocessor_config():
    return importlib.import_module(
        "tests.configs.tau_mq_preprocessor_config", "configs"
    )


@pytest.fixture
def satb1_config():
    return importlib.import_module("tests.configs.satb1_config", "configs")


@pytest.fixture
def satb1_preprocessor_config():
    return importlib.import_module(
        "tests.configs.satb1_mq_preprocessor_config", "configs"
    )


@pytest.mark.parametrize(
    "evidence_df, metadata_df, metadata_column, config_fixture, preprocessor_config_fixture, expected_mods_file, expected_cleavages_file",
    [
        (
            get_evidence_df(GFAP_DATA_DIR),
            None,
            None,
            "gfap_config",
            "gfap_preprocessor_config",
            GFAP_DATA_DIR / "expected_modifications.csv",
            GFAP_DATA_DIR / "expected_cleavages.csv",
        ),
        (
            get_evidence_df(GFAP_DATA_DIR),
            get_metadata_df(GFAP_DATA_DIR),
            "Group",
            "gfap_config",
            "gfap_preprocessor_config",
            GFAP_DATA_DIR / "expected_modifications_with_groups.csv",
            GFAP_DATA_DIR / "expected_cleavages_with_groups.csv",
        ),
        (
            get_evidence_df(TAU_DATA_DIR),
            None,
            None,
            "tau_config",
            "tau_preprocessor_config",
            TAU_DATA_DIR / "expected_modifications.csv",
            TAU_DATA_DIR / "expected_cleavages.csv",
        ),
        (
            get_evidence_df(TAU_DATA_DIR),
            get_metadata_df(TAU_DATA_DIR),
            "Group",
            "tau_config",
            "tau_preprocessor_config",
            TAU_DATA_DIR / "expected_modifications_with_groups.csv",
            TAU_DATA_DIR / "expected_cleavages_with_groups.csv",
        ),
        (
            get_evidence_df(SATB1_DATA_DIR),
            None,
            None,
            "satb1_config",
            "satb1_preprocessor_config",
            SATB1_DATA_DIR / "expected_modifications.csv",
            SATB1_DATA_DIR / "expected_cleavages.csv",
        ),
        (
            get_evidence_df(SATB1_DATA_DIR),
            get_metadata_df(SATB1_DATA_DIR),
            "Group",
            "satb1_config",
            "satb1_preprocessor_config",
            SATB1_DATA_DIR / "expected_modifications_with_groups.csv",
            SATB1_DATA_DIR / "expected_cleavages_with_groups.csv",
        ),
    ],
)
def test_max_quant_preprocessing(
    evidence_df,
    metadata_df,
    metadata_column,
    config_fixture,
    preprocessor_config_fixture,
    expected_mods_file,
    expected_cleavages_file,
    request,
):
    config = request.getfixturevalue(config_fixture)
    preprocessor_config = request.getfixturevalue(preprocessor_config_fixture)

    MaxQuantPreprocessor(
        config,
        preprocessor_config,
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


def test_max_quant_preprocessing_ptms_key_locations(
        gfap_extended_evidence_df, gfap_metadata_df, gfap_extended_config, gfap_preprocessor_config,
):
    # The extended GFAP evidence has "mocked" PTMs at key locations: beginning of sequence, before exons, at exon
    # starts and ends (which are also sequence ends)
    MaxQuantPreprocessor(
        gfap_extended_config,
        gfap_preprocessor_config,
        evidence_df=gfap_extended_evidence_df,
        metadata_df=gfap_metadata_df,
        metadata_column="Group",
    )

    resulting_mods = pd.read_csv("tests/output/result_max_quant_mods.csv")
    resulting_mods = resulting_mods[sorted(resulting_mods.columns)]
    expected_mods = pd.read_csv(GFAP_EXTENDED_DATA_DIR / "expected_modifications_with_groups.csv")
    expected_mods = expected_mods[sorted(expected_mods.columns)]
    pd.testing.assert_frame_equal(
        resulting_mods,
        expected_mods,
        check_dtype=False,
    )

    resulting_cleavages = pd.read_csv("tests/output/result_max_quant_cleavages.csv")
    resulting_cleavages = resulting_cleavages[sorted(resulting_cleavages.columns)]
    expected_cleavages = pd.read_csv(GFAP_EXTENDED_DATA_DIR / "expected_cleavages_with_groups.csv")
    expected_cleavages = expected_cleavages[sorted(expected_cleavages.columns)]
    pd.testing.assert_frame_equal(
        resulting_cleavages,
        expected_cleavages,
        check_dtype=False,
    )


# TODO: would be nice to have a test for the cassette exon, e.g. Tau isoform 7 vs 8


def test_fasta_non_matching_isoform_ids(
    gfap_evidence_df, gfap_config, gfap_preprocessor_config
):
    gfap_preprocessor_config.FASTA_FILE = (
        GFAP_DATA_DIR / "non_matching_isoform_ids.fasta"
    )
    with pytest.raises(
        ValueError,
        match=r"There seem to be isoforms of different proteins in the fasta file.*",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=gfap_evidence_df,
            metadata_df=None,
            metadata_column=None,
        )


def test_fasta_proteins_not_in_evidence_df(
    gfap_evidence_df, gfap_config, gfap_preprocessor_config
):
    gfap_preprocessor_config.FASTA_FILE = GFAP_DATA_DIR / "wrong.fasta"
    with pytest.raises(
        ValueError,
        match="No matching isoform IDs found between the uploaded evidence file and the "
        "fasta file.",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=gfap_evidence_df,
            metadata_df=None,
            metadata_column=None,
        )


def test_malformed_fasta(gfap_evidence_df, gfap_config, gfap_preprocessor_config):
    gfap_preprocessor_config.FASTA_FILE = GFAP_DATA_DIR / "malformed.fasta"
    with pytest.raises(
        ValueError,
        match=r"Error parsing fasta file. Please check the format of the fasta file. "
        r"Uniprot style is recommended.",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=gfap_evidence_df,
            metadata_df=None,
            metadata_column=None,
        )


def test_fasta_shortened_sequence(
    gfap_evidence_df, gfap_config, gfap_preprocessor_config
):
    gfap_preprocessor_config.FASTA_FILE = GFAP_DATA_DIR / "short.fasta"
    with pytest.raises(
        ValueError,
        match=r"The longest original sequence has a length of [0-9]+, but the regions file only "
        r"contains regions up to position [0-9]+. Please adjust the regions in the "
        "corresponding CSV file to match the sequence length.",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=gfap_evidence_df,
            metadata_df=None,
            metadata_column=None,
        )


def test_regions_not_matching_protein(
    gfap_evidence_df, gfap_config, gfap_preprocessor_config
):
    # Removes first exon
    gfap_config.REGIONS = gfap_config.REGIONS[:-2] + [gfap_config.REGIONS[-1]]
    with pytest.raises(
        ValueError,
        match=r"The longest original sequence has a length of .*, but the regions file only "
        r"contains regions up to position .*. Please adjust the regions in the corresponding "
        r"CSV file to match the sequence length.",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=gfap_evidence_df,
            metadata_df=None,
            metadata_column=None,
        )


def test_regions_too_short_but_matching_region_end(
    tau_evidence_df, gfap_config, tau_preprocessor_config
):
    gfap_config.REGIONS = [("Region 1", 44, "A", "1"), ("Region 1", 733, "B", "M")]

    with pytest.raises(
        ValueError,
        match=r"The longest original sequence has a length of .*, but the regions file only "
        r"contains regions up to position .*. Please adjust the regions in the corresponding "
        r"CSV file to match the sequence length.",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            tau_preprocessor_config,
            evidence_df=tau_evidence_df,
            metadata_df=None,
            metadata_column=None,
        )


def test_metadata_file_names_differing(
    gfap_evidence_df, gfap_metadata_df, gfap_config, gfap_preprocessor_config
):
    gfap_metadata_df["Sample"] = gfap_metadata_df["Sample"].apply(
        lambda x: f"Different_{x}"
    )
    with pytest.raises(
        ValueError,
        match=r"The following samples from the evidence file are missing in the metadata file: .*",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=gfap_evidence_df,
            metadata_df=gfap_metadata_df,
            metadata_column="Group",
        )


def test_no_groups_in_metadata(
    gfap_evidence_df, gfap_metadata_df, gfap_config, gfap_preprocessor_config
):
    metadata_df = pd.DataFrame(columns=gfap_metadata_df.columns)
    with pytest.raises(
        ValueError,
        match=r"The following samples from the evidence file are missing in the metadata file: .*",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=gfap_evidence_df,
            metadata_df=metadata_df,
            metadata_column="Group",
        )


def test_group_col_not_in_metadata_df(
    gfap_evidence_df, gfap_metadata_df, gfap_config, gfap_preprocessor_config
):
    with pytest.raises(
        ValueError,
        match=r"Metadata column '.*' not found in metadata DataFrame columns: .*",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=gfap_evidence_df,
            metadata_df=gfap_metadata_df,
            metadata_column="NonExistingColumn",
        )


def test_evidence_samples_not_matching_metadata(
    gfap_evidence_df, gfap_metadata_df, gfap_config, gfap_preprocessor_config
):
    new_evidence_df = gfap_evidence_df
    new_evidence_df["Sample"] = new_evidence_df["Sample"].apply(
        lambda x: f"Different_{x}" if "AD" in x or "CTR" in x else x
    )

    with pytest.raises(
        ValueError,
        match=f"The following samples from the evidence file are missing in the metadata file: .*",
    ):
        MaxQuantPreprocessor(
            gfap_config,
            gfap_preprocessor_config,
            evidence_df=new_evidence_df,
            metadata_df=gfap_metadata_df,
            metadata_column="Group",
        )
