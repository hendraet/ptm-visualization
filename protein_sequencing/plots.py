"""
Main module to create different plots.
Execute with python3 plots.py -p <plot_type> -f <fasta_file> -o <output_folder>
Optional arguments:
    -pc, --plot-config: Path to plot specific configuration file.
    -c, --config: Path to configuration file.
"""

import argparse
import importlib
from protein_sequencing import utils, sequence_plot
from protein_sequencing.bar_plot import BarPlotter
from protein_sequencing.details_plot import DetailsPlotter
from protein_sequencing.overview_plot import OverviewPlotter


def generate_bar_plot(config, plot_config, fasta, output):
    """Generate bar plot."""
    BarPlotter(importlib.import_module(config, 'configs'),
               importlib.import_module(plot_config, 'configs'),
               fasta, output)


def generate_details_plot(config, plot_config, fasta, output):
    """Generate details plot."""
    DetailsPlotter(importlib.import_module(config, 'configs'),
                   importlib.import_module(plot_config, 'configs'),
                   fasta, output)


def generate_overview_plot(config, plot_config, fasta, output):
    """Generate overview plot."""
    OverviewPlotter(importlib.import_module(config, 'configs'),
                    importlib.import_module(plot_config, 'configs'),
                    fasta, output)


DEFAULT_CONFIGS = {
    'bar': 'configs.default_bar',
    'details': 'configs.default_details',
    'overview': 'configs.default_overview',
    'config': 'configs.default_config',
}


def main():
    """Main function to generate protein sequencing plots."""
    parser = argparse.ArgumentParser(description='Generate protein sequencing plots.')

    parser.add_argument(
        '-p', '--plot',
        required=True,
        choices=DEFAULT_CONFIGS.keys(),  # TODO: so 'config' would apparently also be a valid plot type
        help='Type of plot to generate (bar, details, overview).'
    )
    parser.add_argument(
        '-pc', '--plot-config',
        required=False,
        help='Path to plot specific configuration file. Default=configs.default_<plot_type>'
    )
    parser.add_argument('-c',
                        '--config',
                        required=False,
                        default=DEFAULT_CONFIGS['config'],
                        help='Path to configuration file. Default=configs.default_config')
    parser.add_argument('-f',
                        '--fasta',
                        required=True,
                        help='Path to Fasta file (e.g., data/uniprot_data/tau_isoforms2N4R.fasta)')
    parser.add_argument('-o', '--output',
                        required=False,
                        default='output',
                        help='Path to output folder, default=output')
    args = parser.parse_args()

    if args.plot_config:
        plot_config = args.plot_config
    else:
        plot_config = DEFAULT_CONFIGS[args.plot]

    sequence_plot.CONFIG = importlib.import_module(args.config, 'configs')
    utils.CONFIG = importlib.import_module(args.config, 'configs')

    if args.plot == 'bar':
        generate_bar_plot(args.config, plot_config, args.fasta, args.output)
    elif args.plot == 'details':
        generate_details_plot(args.config, plot_config, args.fasta, args.output)
    elif args.plot == 'overview':
        generate_overview_plot(args.config, plot_config, args.fasta, args.output)
    else:
        print(f"Unknown plot type: {args.plot}. Please choose from 'bar', 'details', 'overview'.")


if __name__ == '__main__':
    main()
