################################################################################
# Runs a specific preprocessing task to properly merge genomic datasets 
# in .vcf format
################################################################################

import argparse
import sys

from utils.logger import get_logger, parser_msg
from utils.misc import check_arguments, check_configurations
from modules.partition_and_rename_chr import partition_by_chromosome

# Define parser to read from command line
parser = argparse.ArgumentParser(usage=parser_msg())

# Define subparser to choose a command
subparsers = parser.add_subparsers(help='Choose a command', dest='command')

# Define subparser for partition command
partition_parser = subparsers.add_parser('partition', help='To partition the data in a separate file per chromosome.')
partition_parser.add_argument('-f', '--file', type=str, required=True, help='Path to input .vcf file.')
partition_parser.add_argument('-o', '--output-folder', type=str, required=True, help='Path to output folder.')
partition_parser.add_argument('-r', '--rename-chr', action='store_true', help='Rename (or not) chromosome notation.')
partition_parser.add_argument('-d', '--debug', type=str, required=False, help='Path to file storing info/debug messages.')

# Parse the arguments
args = parser.parse_args()

# Create logger
logger = get_logger(__name__, args.debug)

if args.command == 'partition':
    
    # Check the input file exists and its format is supported
    check_arguments([args.file])
    
    # Partition the file in a separate file per chromosome
    partition_by_chromosome(args.file, args.output_folder, args.rename_chr, logger)
    

    
