################################################################################
# Utilities to define parser and check the validity of arguments
# from the Command Line Interface (CLI).
################################################################################

import os
import argparse
from typing import Dict, List

from utils.vcf_utils import obtain_chromosomes

def parser_msg() -> str:
    """
    Specifies the custom usage of the parser.
    
    Returns
        (str)
    
    """
    return '''MergeGenome.py
        partition
            [-f, --file: path to input .vcf file)
            [-o, --output-folder: path to output folder]
            [-r, --rename-chr: rename (or not) chromosome notation]
            [-d, --debug: path to file storing info/debug messages]
        '''

def define_parser() -> argparse.ArgumentParser:
    """
    Defines the parser to read the arguments from the Command 
    Line Interface (CLI). Each command has its own parser arguments.
    
    Args:
        (None)
    
    Returns:
        parser (argparse.ArgumentParser): object to parse arguments 
        from CLI.
        
    """
    
    # Define parser to read arguments from CLI
    parser = argparse.ArgumentParser(usage=parser_msg())

    # Define subparser for the commands
    subparsers = parser.add_subparsers(help='Choose a command', dest='command')
    
    # Define subparser for 'partition' command
    partition_parser = subparsers.add_parser('partition', help='To partition .vcf data in a separate .vcf file per chromosome.')
    partition_parser.add_argument('-q', '--query', required=True, 
                                  help='Path to input .vcf file with data for multiple chromosomes.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-r', '--rename-chr', action='store_true', help='To rename chromosome notation.')
    partition_parser.add_argument('-m', '--rename-map', required=False, help='Mapping from actual to new chromosome notation.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')
    
    # Define subparser for 'rename' command
    partition_parser = subparsers.add_parser('rename', help='To rename chromosome notation (variants/CHROM field).')
    partition_parser.add_argument('-q', '--query', required=True, 
                                  help='Path to input .vcf file with data for a single or multiple chromosomes.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-m', '--rename-map', required=False, 
                                  help='Mapping from actual to new chromosome notation.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')
    
    # Define subparser for 'clean' command
    partition_parser = subparsers.add_parser('clean', help='To clean genomic sequences.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for a single chromosome each.')
    partition_parser.add_argument('-r', '--reference', required=False, nargs="*", 
                                  help='Paths to reference .vcf files with data for a single chromosome each.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-s', '--remove-sample-ID-query', required=False, nargs="*", 
                                  help='Sample IDs or substring of sample IDs to remove from the query.')
    partition_parser.add_argument('-t', '--remove-sample-ID-reference', required=False, nargs="*", 
                                  help='Sample IDs or substring of sample IDs to remove from the reference.')
    partition_parser.add_argument('-a', '--remove-ambiguous-snps-query', action='store_true', 
                                  help='To remove (or not) ambiguous SNPs from the query.')
    partition_parser.add_argument('-b', '--remove-ambiguous-snps-reference', action='store_true', 
                                  help='To remove (or not) ambiguous SNPs from the reference.')
    partition_parser.add_argument('-f', '--correct-snp-flips', action='store_true', 
                                  help='To correct (or not) SNP flips in the query with respect to the reference.')
    partition_parser.add_argument('-m', '--remove-mismatching-snps', action='store_true', help='To remove (or not) mismatching SNPs.')
    partition_parser.add_argument('-v', '--rename-map-query', required=False, 
                                  help='Mapping from old to new missing values notation for the query.')
    partition_parser.add_argument('-w', '--rename-map-reference', required=False, 
                                  help='Mapping from old to new missing values notation for the reference.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')

    # Define subparser for 'subset' command
    partition_parser = subparsers.add_parser('subset', help='To subset the data to the common markers.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for each chromosome.')
    partition_parser.add_argument('-r', '--reference', required=True, nargs="*", 
                                  help='Paths to reference .vcf files with data for each chromosome.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')

    # Define subparser for 'plot-snp-means' command
    partition_parser = subparsers.add_parser('plot-snp-means', 
                                             help='To plot the SNP means for the common markers between the query and the reference.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for each chromosome.')
    partition_parser.add_argument('-r', '--reference', required=True, nargs="*", 
                                  help='Paths to reference .vcf files with data for each chromosome.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-x', '--x-axis-name', required=False, default="query", 
                                  help='Name given to the query dataset that will appear in the x-axis.')
    partition_parser.add_argument('-y', '--y-axis-name', required=False, default="reference", 
                                  help='Name given to the reference dataset that will apprear in the y-axis.')
    partition_parser.add_argument('-f', '--fontsize', required=False, default=25, help='Fontsize in the plot.')
    partition_parser.add_argument('-w', '--figure-width', required=False, default=26, help='Figure width of the plot.')
    partition_parser.add_argument('-i', '--figure-height', required=False, default=15, help='Figure height of the plot.')
    partition_parser.add_argument('-s', '--size-points', required=False, default=0.1, help='Size of the points in the plot.')
    partition_parser.add_argument('-c', '--color-points', required=False, default='#306998', help='Color of the points in the plot.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')

    # Define subparser for 'plot-pca' command
    partition_parser = subparsers.add_parser('plot-pca', help='To plot the PCA of the SNPs data.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for each chromosome.')
    partition_parser.add_argument('-r', '--reference', required=True, nargs="*", 
                                  help='Paths to reference .vcf files with data for each chromosome.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-t', '--train-both', action='store_true', help='Train on both the query and the reference.')
    partition_parser.add_argument('-f', '--fontsize', required=False, default=25, help='Fontsize in the plot.')
    partition_parser.add_argument('-w', '--figure-width', required=False, default=26, help='Figure width of the plot.')
    partition_parser.add_argument('-i', '--figure-height', required=False, default=15, help='Figure height of the plot.')
    partition_parser.add_argument('-s', '--size-points', required=False, default=0.1, help='Size of the points in the plot.')
    partition_parser.add_argument('-cq', '--color-points-query', required=False, default='#259988', 
                                  help='Color of query points in the plot.')
    partition_parser.add_argument('-cr', '--color-points-reference', required=False, default='#EBD0A1', 
                                  help='Color of reference points in the plot.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')

    # Define subparser for 'store-allele' command
    partition_parser = subparsers.add_parser('store-allele', help='To store formatted allele data in .npy or .h5 format.')
    partition_parser.add_argument('-q', '--query', required=True, help='Path to input .vcf file.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-s', '--data-format', required=True, choices=['separated', 'averaged'], 
                                  help='Separate or average maternal and paternal strands.')
    partition_parser.add_argument('-f', '--file-format', required=True, choices=['.npy', '.h5'], 
                                  help='Format of the output file.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')

    # Define subparser for 'store-indexes-common' command
    partition_parser = subparsers.add_parser('store-indexes-common', help='Store indexes of common markers.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for each chromosome.')
    partition_parser.add_argument('-r', '--reference', required=False, nargs="*", 
                                  help='Paths to reference .vcf files with data for each chromosome.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-f', '--file-format', required=True, choices=['.npy', '.h5'], help='Format of the output file.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')

    # Define subparser for 'remove-snps-different-means' command
    partition_parser = subparsers.add_parser('remove-snps-different-means', help='Store indexes of common markers.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for each chromosome.')
    partition_parser.add_argument('-r', '--reference', required=False, nargs="*", 
                                  help='Paths to reference .vcf files with data for each chromosome.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-t', '--threshold', required=False, default=0.1, 
                                  help='All common SNPs with a mean absolute difference higher than the threshold will be removed.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')
    
    return parser


def check_paths(paths: List[str]) -> None:
    """
    Checks that all the file paths and formats are correct.
    
    Args:
        paths (List[str]): paths to files.
    
    Returns:
        (None)
    """

    for path in paths:
        # For each path...
        # Check the path exists
        assert os.path.isfile(path), f'{path} not found. Is it a file?'
        
        # Check the path is a .vcf file
        assert path.endswith('.vcf'), f'{path} is not in a supported format. Convert it to .vcf file.'
    

def check_arguments(args: argparse.Namespace) -> None:
    """
    Checks that the file paths and formats are correct. 
    Also checks that all required arguments provided 
    are appropriate for each command in particular.
    
    Args:
        args (argparse.Namespace): parser with all arguments
        information.
    
    Returns:
        (None)
    """
    
    if args.command == 'partition':
        # Check the query path exists and is a .vcf file
        check_paths([args.query])
    
        if args.rename_map is not None: 
            if args.rename_chr == False:
                # Ensure if rename_map is provided rename_chr = True
                raise argparse.ArgumentTypeError(f'--rename_map is only supported when --rename_chr=True.')
            # Convert args.rename_map str to dict
            args.rename_map = eval(args.rename_map)
        
    if args.command == 'rename':
        # Check the query path exists and is a .vcf file
        check_paths([args.query])
        
        if args.rename_map is not None:
            # Convert args.rename_map str to dict
            args.rename_map = eval(args.rename_map)
        
    elif args.command == 'clean':
        # Check all the query paths exist and are a .vcf file
        check_paths(args.query)
        
        if args.reference == None:
            # Put reference in list to be compatible with other scripts
            args.reference = [None]
            
            # Ensure all tasks specified only work on the query
            if remove_sample_ID_reference is not None:
                raise argparse.ArgumentTypeError(f'--remove_sample_ID_reference is only supported when a '\
                                                 'reference is provided.')
            if remove_ambiguous_snps_reference is not None:
                raise argparse.ArgumentTypeError(f'--remove_ambiguous_snps_reference is only supported when a '\
                                                 'reference is provided.')
            if correct_snp_flips is not None:
                raise argparse.ArgumentTypeError(f'--correct_snp_flips is only supported when a reference is '\
                                                 'provided.')
            if remove_mismatching_snps is not None:
                raise argparse.ArgumentTypeError(f'--remove_mismatching_snps is only supported when a '\
                                                 'reference is provided.')
            if rename_map_reference is not None:
                raise argparse.ArgumentTypeError(f'--rename_map_reference is only supported when a '\
                                                 'reference is provided.')
        else:
            # Check all the reference paths exist and are a .vcf file
            check_paths(args.reference)
        if args.rename_map_query is not None:
            # Convert args.rename_map_query str to dict
            args.rename_map_query = eval(args.rename_map_query)
        if args.rename_map_reference is not None:
            # Convert args.rename_map_reference str to dict
            args.rename_map_reference = eval(args.rename_map_reference)
            
    
def check_chromosome(query: Dict, reference: Dict):
    """
    Ensures the query and the reference contain data for a single
    chromosome and that such chromosome is the same in both datasets.
    Returns the chromosome in particular.
    
    Args:
        query (Dict): query data.
        reference (Dict): reference data.
    
    Returns:
        chrom (str): chromosome in the query and the reference.
    
    """
    
    # Obtain list with chromosomes in the query
    chroms_query = obtain_chromosomes(query)
    
    if reference is not None:
        # If the reference is not None...
        # Obtain list with chromosomes in the reference
        chroms_reference = obtain_chromosomes(reference)

        # Ensure the query and the reference contain data for a single chromosome
        assert len(chroms_query) == 1, 'The query contains data for more than one chromosome.'\
        'Use partition command to obtain a separate VCF file per chromosome.'
        assert len(chroms_reference) == 1, 'The reference contains data for more than one chromosome.'\
        'Use partition command to obtain a separate VCF file per chromosome.'

        # Ensure the chromosome is the saame
        assert chroms_query[0] == chroms_reference[0], f'The query contains data for chromosome {chroms_query[0]}.'\
        f'while the query contains data for chromosome {chroms_reference[0]}, when they must be the same.'\
        'Check the input files.'
    
    return chroms_query[0]
    
    
    