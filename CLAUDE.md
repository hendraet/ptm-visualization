# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

All commands run from the repo root.

- Install deps: `poetry install`
- Make Clustal Omega executable (one-time, Linux): `chmod +x protein_sequencing/clustal-omega/clustalo-1.2.4-Ubuntu-x86_64`
- Run tests: `poetry run python -m pytest tests/test_max_quant_preprocessor.py`
- Run a single test: `poetry run python -m pytest tests/test_max_quant_preprocessor.py::test_malformed_fasta`
- Generate a plot: `poetry run python -m protein_sequencing.plots -p {bar|details|overview} -f PATH/TO/FASTA [-o OUTPUT]`
- Run a preprocessor: `poetry run python -m protein_sequencing.data_preprocessing.preprocessor -p {ma|pp|mq|ms}`

Use `python -m`, not `python protein_sequencing/plots.py` — the modules use absolute `protein_sequencing.*` imports and only resolve with the repo root on `sys.path`. Tests likewise import `tests.*` and `protein_sequencing.*` as packages, so they must be invoked from the root.

## Architecture

This is a two-stage pipeline that visualizes post-translational modifications (PTMs) and cleavages on protein sequences.

**Stage 1 — Preprocessors** (`protein_sequencing/data_preprocessing/`) parse mass-spec search-engine output and emit CSV files (`<output>/result_<tool>_mods.csv` and `result_<tool>_cleavages.csv`). One preprocessor class per tool — `MaxQuantPreprocessor`, `MascotPreprocessor`, `ProteinPilotPreprocessor`, `MSFraggerPreprocessor` — all dispatched via `data_preprocessing/preprocessor.py`. Shared logic lives in `preprocessor_helper.py`. MaxQuant is the actively maintained/tested path; the others are older.

**Stage 2 — Plotters** (`protein_sequencing/`) read those CSVs and render Plotly figures. `Plotter` (`plotter.py`) is the base class handling sequence/region/exon layout and pixel math; `BarPlotter`, `DetailsPlotter`, `OverviewPlotter` subclass it. Dispatched via `plots.py`.

**Constructors do the work.** Both preprocessor and plotter classes run the entire pipeline as a side effect of `__init__` — instantiating the class writes the output files. There is no separate `.run()`.

### Configuration

Behavior is driven entirely by Python config modules under `configs/`, loaded by name with `importlib.import_module(...)`. Two config kinds are passed to every entry point:

- A **general config** (`default_config.py`) — `REGIONS`, `MODIFICATIONS`, `INCLUDED_MODIFICATIONS`, figure dimensions, `OUTPUT_FOLDER`, `MIN_EXON_LENGTH`.
- A **stage-specific config** — `preprocessor_config.py` (FASTA path, thresholds) or `default_{bar,details,overview}.py` (plot tuning). `max_quant_*.py` are MaxQuant-specific variants.

`REGIONS` is a list of `(name, end_position, group, abbreviation)` tuples that must collectively span the full protein length, or the preprocessor raises a length-mismatch `ValueError`.

### Sequence alignment and exons

`uniprot_align.py` aligns FASTA isoforms with the bundled Clustal Omega binary (platform-specific binaries under `clustal-omega/`; single-isoform files skip alignment). `exon_helper.py` detects alternative/cassette exons. Exon handling is the main source of complexity: it threads through `Plotter._create_plot` (region/pixel layout, `region_plot_type` 0–5) and the preprocessors' offset math (`calculate_exon_offset`, `get_isoform_for_offset`). When changing exon logic, check both layers.

### Metadata / groups

Preprocessors optionally take a `metadata_df` + `metadata_column` to split modifications into experimental groups (e.g. AD vs. control), producing per-group columns in the output CSVs. Sample names in the evidence file must all be present in the metadata file.

## Testing

`tests/test_max_quant_preprocessor.py` is the current suite — it runs `MaxQuantPreprocessor` against fixtures in `tests/test_data/` (GFAP, TAU, SATB1 proteins) and compares the CSVs in `tests/output/` against checked-in `expected_*.csv` files. Test configs live in `tests/configs/`. `tests/preprocessor_tests.py` is an older suite for the other preprocessors.

To regenerate expected fixtures after an intentional output change, compare against `tests/output/result_max_quant_*.csv` produced by a test run before overwriting the `expected_*` files.