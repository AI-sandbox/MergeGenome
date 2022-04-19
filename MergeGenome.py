################################################################################
# Runs a specific preprocessing task to properly merge genomic datasets 
# in .vcf format
################################################################################

import argparse
import sys

from utils.logger import get_logger, parser_msg
from utils.misc import check_arguments, check_configurations
from modules.partition import partition_by_chromosome
from modules.rename import rename_chromosome
from modules.clean import clean_genomic_data

# Define parser to read from command line
parser = argparse.ArgumentParser(usage=parser_msg())

# Define subparser to choose a command
subparsers = parser.add_subparsers(help='Choose a command', dest='command')

# Define subparser for partition command
partition_parser = subparsers.add_parser('partition', help='To partition the data in a separate file per chromosome.')
partition_parser.add_argument('-f', '--file', required=True, help='Path to input .vcf file.')
partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder to store separate .vcf files.')
partition_parser.add_argument('-r', '--rename-chr', action='store_true', help='Rename chromosome notation.')
partition_parser.add_argument('-m', '--rename-map', required=False, help='Dictionary with mapping from actual to new chromosome notation.')
partition_parser.add_argument('-d', '--debug', required=False, help='Path to file to store info/debug messages.')

# Define subparser for rename command
partition_parser = subparsers.add_parser('rename', help='.')
partition_parser.add_argument('-f', '--file', required=True, help='Path to .vcf file with data for a particular chromosome.')
partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder to store the modified .vcf file.')
partition_parser.add_argument('-m', '--rename-map', required=False, help='Dictionary with mapping from actual to new chromosome notation.')
partition_parser.add_argument('-d', '--debug', required=False, help='Path to file to store info/debug messages.')

# Define subparser for clean command
partition_parser = subparsers.add_parser('clean', help='.')
partition_parser.add_argument('-r', '--reference', required=True, nargs="*", help='Paths to reference .vcf files with data for each chromosome.')
partition_parser.add_argument('-q', '--query', required=True, nargs="*", help='Paths to query .vcf files with data for each chromosome.')
partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder to store the modified .vcf files.')
partition_parser.add_argument('-s', '--remove-sample-ID', required=False, nargs="*", help='Sample IDs or substring of sample IDs to remove.')
partition_parser.add_argument('-a', '--remove-ambiguous-snps', action='store_true', help='Remove (or not) ambiguous SNPs.')
partition_parser.add_argument('-f', '--correct-snp-flips', action='store_true', help='Correct (or not) SNP in the query with respect to the reference.')
partition_parser.add_argument('-m', '--remove-mismatching-snps', action='store_true', help='Remove (or not) mismatching SNPs.')
partition_parser.add_argument('-v', '--rename-map-query', required=False, help='Dictionary with mapping from old to new missing notation for the '\
                              'query.')
partition_parser.add_argument('-w', '--rename-map-reference', required=False, help='Dictionary with mapping from old to new missing notation for the'\
                              'reference.')
partition_parser.add_argument('-d', '--debug', required=False, help='Path to file to store info/debug messages.')

# Parse the arguments
args = parser.parse_args()

# Create logger
logger = get_logger(__name__, args.debug)

if args.command == 'partition':
    
    # Check the input arguments are correct
    check_arguments([args.file], args.rename_chr, args.rename_map)
    
    # Partition the file in a separate file per chromosome
    partition_by_chromosome(args.file, args.output_folder, args.rename_chr, args.rename_map, logger)
    

if args.command == 'rename':
    
    # Check the input arguments are correct
    check_arguments([args.file])
    
    # Rename chromosomme notation
    rename_chromosome(args.file, args.output_folder, args.rename_map, logger)
    
if args.command == 'clean':
    
    # Check the input arguments are correct
    check_arguments(args.reference+args.query)
    
    # Clean genomic data
    clean_genomic_data(args.reference, args.query, args.output_folder, args.remove_sample_ID, args.remove_ambiguous_snps, 
                       args.correct_snp_flips, args.remove_mismatching_snps, args.rename_map_query, 
                       args.rename_map_reference, logger)
    
    