import argparse
from protein_sequencing import bar_plot, details_plot, overview_plot

# Define functions for each plot type
def generate_bar_plot(config_path, fasta, output):
    bar_plot.create_bar_plot(fasta, output, config_path)

def generate_details_plot(config_path, fasta, output):
    details_plot.create_details_plot(fasta, output)

def generate_overview_plot(config_path, fasta, output):
    overview_plot.create_overview_plot(fasta, output)


DEFAULT_CONFIGS = {
    'bar': './configs/default_bar.yml',
    'details': './configs/default_details.yml',
    'overview': './configs/default_overview.yml'
}

def main():
    parser = argparse.ArgumentParser(description='Generate protein sequencing plots.')
    
    parser.add_argument(
        '-p', '--plot', 
        required=True,
        choices=DEFAULT_CONFIGS.keys(),
        help='Type of plot to generate (bar, details, overview).'
    )
    parser.add_argument(
        '-c', '--config', 
        help='Path to configuration file. If not provided, a default config file will be used based on the plot type.'
    )
    parser.add_argument('-f',
                         '--fasta',
                         required=False,
                         default='data/uniprot_data/tau_isoforms2N4R.fasta',
                         help='Path to Fasta file (e.g., data/uniprot_data/tau_isoforms2N4R.fasta)')
    parser.add_argument('-o', '--output', required=False, default='output', help='Path to output folder, default=output')
    args = parser.parse_args()

    if args.config:
        config_path = args.config
    else:
        config_path = DEFAULT_CONFIGS[args.plot]

    if args.plot == 'bar':
        generate_bar_plot(config_path, args.fasta, args.output)
    elif args.plot == 'details':
        generate_details_plot(config_path, args.fasta, args.output)
    elif args.plot == 'overview':
        generate_overview_plot(config_path, args.fasta, args.output)
    else:
        print(f"Unknown plot type: {args.plot}. Please choose from 'bar', 'details', 'overview'.")

if __name__ == '__main__':
    main()
