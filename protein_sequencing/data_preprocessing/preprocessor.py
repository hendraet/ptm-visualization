"""Common interface to execute preprocessors."""
import argparse
import importlib
from protein_sequencing.data_preprocessing.protein_pilot_preprocessor import ProteinPilotPreprocessor
from protein_sequencing.data_preprocessing.mascot_preprocessor import MascotPreprocessor
from protein_sequencing.data_preprocessing.ms_fragger_preprocessor import MSFraggerPreprocessor
from protein_sequencing.data_preprocessing.max_quant_preprocessor import MaxQuantPreprocessor

def mascot(config, pre_config):
    """Mascot preprocessor."""
    MascotPreprocessor(importlib.import_module(config, 'configs'), importlib.import_module(pre_config, 'configs'))

def protein_pilot(config, pre_config):
    """Protein Pilot preprocessor."""
    ProteinPilotPreprocessor(importlib.import_module(config, 'configs'), importlib.import_module(pre_config, 'configs'))

def ms_fragger(config, pre_config):
    """MS Fragger preprocessor."""
    MSFraggerPreprocessor(importlib.import_module(config, 'configs'), importlib.import_module(pre_config, 'configs'))

def max_quant(config, pre_config):
    """MaxQuant preprocessor."""
    MaxQuantPreprocessor(importlib.import_module(config, 'configs'), importlib.import_module(pre_config, 'configs'))

DEFAULT_CONFIGS = {
    'config': 'configs.default_config',
    'preprocessor': 'configs.preprocessor_config',
}

def main():
    """Main function to generate protein sequencing plots."""
    parser = argparse.ArgumentParser(description='Generate protein sequencing plots.')

    parser.add_argument(
        '-p', '--preprocessor', 
        required=True,
        choices=['ma', 'pp', 'mq', 'ms'],
        help='Type of preprocessor (ma (Mascot), pp (ProteinPilot), mq (MaxQuant), ms (MS Fragger)).'
    )
    parser.add_argument(
        '-pc', '--preprocessor-config',
        required=False,
        default=DEFAULT_CONFIGS['preprocessor'],
        help='Path to plot specific configuration file. Default=configs.default_preprocessor'
    )
    parser.add_argument('-c',
                        '--config',
                        required=False,
                        default=DEFAULT_CONFIGS['config'],
                        help='Path to configuration file. Default=configs.default_config')
    args = parser.parse_args()


    if args.preprocessor == 'ma':
        mascot(args.config, args.preprocessor_config)
    elif args.preprocessor == 'pp':
        protein_pilot(args.config, args.preprocessor_config)
    elif args.preprocessor == 'mq':
        max_quant(args.config, args.preprocessor_config)
    elif args.preprocessor == 'ms':
        ms_fragger(args.config, args.preprocessor_config)
    else:
        print(f"Unknown preprocessor type: {args.preprocessor}. Currently supported preprocessors are: ma (Mascot), pp (ProteinPilot), mq (MaxQuant), ms (MS Fragger).")

if __name__ == '__main__':
    main()
