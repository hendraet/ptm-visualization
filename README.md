# ptm-visualization
## Setup
1. Check out repository and cd to it
1. Set up a virtualenv for the project with Python >=3.11 and activate it (e.g., `python3 -m venv venv` and `source venv/bin/activate`)
1. Install poetry (if not already installed): `curl -sSL https://install.python-poetry.org/ | python -`
1. Install dependencies with `poetry install`
1. Make the clustal-omega (`http://www.clustal.org/omega/`) file executable with `chmod +x clustal-omega\clustalo-1.2.4-Ubuntu-x86_64`
## Run
1. To run a preprocessor, you must execute the corresponding script by running, e.g., `python3 protein_pilot_preprocessor.py`. Be sure to supply a FASTA file, group.csv and the `preprocessor_config.py`.
1. To run the plotting script, run with `python3 plots.py -p PLOT_TYPE -f PATH/TO/FASTA`. The plot type can be `overview,` `bar`, or `details`. Be sure to alter the settings to your needs in the configuration files.