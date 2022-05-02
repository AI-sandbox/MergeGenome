################################################################################
# Runs a specific preprocessing task to properly merge genomic datasets 
# in VCF format
################################################################################

import argparse
import sys

from utils.logger import get_logger, parser_msg
from utils.misc import check_arguments, check_configurations
from modules.partition import partition_by_chromosome
from modules.rename import rename_chromosome
from modules.clean import clean_genomic_data
from modules.subset import subset_common_markers
from scripts.store_allele_data import store_allele_data_in_npy_or_h5

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
partition_parser = subparsers.add_parser('rename', help='To rename chromosome notation.')
partition_parser.add_argument('-f', '--file', required=True, help='Path to .vcf file with data for a particular chromosome.')
partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder to store the modified .vcf file.')
partition_parser.add_argument('-m', '--rename-map', required=False, help='Dictionary with mapping from actual to new chromosome notation.')
partition_parser.add_argument('-d', '--debug', required=False, help='Path to file to store info/debug messages.')

# Define subparser for clean command
partition_parser = subparsers.add_parser('clean', help='To clean genomic sequences.')
partition_parser.add_argument('-q', '--query', required=True, nargs="*", help='Paths to query .vcf files with data for each chromosome.')
partition_parser.add_argument('-r', '--reference', required=False, nargs="*", help='Paths to reference .vcf files with data for each chromosome.')
partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder to store the modified .vcf files.')
partition_parser.add_argument('-s', '--remove-sample-ID-query', required=False, nargs="*", help='Sample IDs or substring of sample IDs to remove from the query.')
partition_parser.add_argument('-t', '--remove-sample-ID-reference', required=False, nargs="*", help='Sample IDs or substring of sample IDs to remove from the reference.')
partition_parser.add_argument('-a', '--remove-ambiguous-snps-query', action='store_true', help='Remove (or not) ambiguous SNPs from the query.')
partition_parser.add_argument('-b', '--remove-ambiguous-snps-reference', action='store_true', help='Remove (or not) ambiguous SNPs from the reference.')
partition_parser.add_argument('-f', '--correct-snp-flips', action='store_true', help='Correct (or not) SNP in the query with respect to the reference.')
partition_parser.add_argument('-m', '--remove-mismatching-snps', action='store_true', help='Remove (or not) mismatching SNPs.')
partition_parser.add_argument('-v', '--rename-map-query', required=False, help='Dictionary with mapping from old to new missing notation for the query.')
partition_parser.add_argument('-w', '--rename-map-reference', required=False, help='Dictionary with mapping from old to new missing notation for the reference.')
partition_parser.add_argument('-d', '--debug', required=False, help='Path to file to store info/debug messages.')

# Define subparser for subset command
partition_parser = subparsers.add_parser('subset', help='To subset the data to the common markers.')
partition_parser.add_argument('-q', '--query', required=True, nargs="*", help='Paths to query .vcf files with data for each chromosome.')
partition_parser.add_argument('-r', '--reference', required=True, nargs="*", help='Paths to reference .vcf files with data for each chromosome.')
partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder to store the modified .vcf files.')
partition_parser.add_argument('-d', '--debug', required=False, help='Path to file to store info/debug messages.')

# Define subparser for store-allele command
partition_parser = subparsers.add_parser('store-allele', help='To store formatted allele data in .npy or .h5 format.')
partition_parser.add_argument('-q', '--query', required=True, help='Path to input .vcf file.')
partition_parser.add_argument('-s', '--data-format', required=True, choices=['separated', 'averaged.'], 
                              help='Separate or average maternal and paternal strands.')
partition_parser.add_argument('-f', '--file-format', required=True, choices=['.npy', '.h5'], help='Format of the output file.')
partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder to store the formatted allele data.')
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
    check_arguments(args.query)
    if args.reference is not None:
        check_arguments(args.reference)
    
    # Clean genomic data
    clean_genomic_data(args.query, args.reference, args.output_folder, args.remove_sample_ID_query, args.remove_sample_ID_reference,
                       args.remove_ambiguous_snps_query, args.remove_ambiguous_snps_reference,
                       args.correct_snp_flips, args.remove_mismatching_snps, args.rename_map_query, 
                       args.rename_map_reference, logger)
    
if args.command == 'subset':
    
    # Check the input arguments are correct
    check_arguments(args.query + args.reference)
    
    # Subset genomic data
    subset_common_markers(args.reference, args.query, args.output_folder, logger)
    
if args.command == 'store-allele':
    
    # Check the input arguments are correct
    check_arguments([args.query])
    
    # Store averaged or separated maternal and paternal strands in .npy or .h5
    store_allele_data_in_npy_or_h5(args.query, args.output_folder, args.data_format, args.file_format, logger)
    
    