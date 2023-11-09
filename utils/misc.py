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
    partition_parser.add_argument('-m', '--remove-mismatching-snps', action='store_true', 
                                  help='To remove (or not) mismatching SNPs.')
    partition_parser.add_argument('-v', '--rename-map-query', required=False, 
                                  help='Mapping from old to new missing values notation for the query.')
    partition_parser.add_argument('-w', '--rename-map-reference', required=False, 
                                  help='Mapping from old to new missing values notation for the reference.')
    partition_parser.add_argument('-d', '--debug', required=False, 
                                  help='Path to .log/.txt file to store info/debug messages.')

    # Define subparser for 'subset' command
    partition_parser = subparsers.add_parser('subset', help='To subset the SNPs to the common markers.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for a single chromosome each.')
    partition_parser.add_argument('-r', '--reference', required=True, nargs="*", 
                                  help='Paths to reference .vcf files with data for a single chromosome each.')
    partition_parser.add_argument('-o', '--output-folder', required=True, 
                                  help='Path to output folder.')
    partition_parser.add_argument('-d', '--debug', required=False, 
                                  help='Path to .log/.txt file to store info/debug messages.')
    
    # Define subparser for 'plot-snp-means' command
    partition_parser = subparsers.add_parser('plot-snp-means', 
                                             help='To plot the SNP means of the common markers between the query and the reference.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for a single or multiple chromosomes each.')
    partition_parser.add_argument('-r', '--reference', required=True, nargs="*", 
                                  help='Paths to reference .vcf files with data for a single or multiple chromosomes each.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-idx', '--indexes', required=False, help='Path to SNP indexes with different color.')
    partition_parser.add_argument('-x', '--x-axis-name', required=False, default="reference", 
                                  help='Name given to reference dataset (x-axis).')
    partition_parser.add_argument('-y', '--y-axis-name', required=False, default="query", 
                                  help='Name given to query dataset (y-axis).')
    partition_parser.add_argument('-f', '--fontsize', required=False, type=float, default=35, help='Fontsize of all text in plot.')
    partition_parser.add_argument('-w', '--figure-width', required=False, type=float, default=20, help='Figure width of plot.')
    partition_parser.add_argument('-i', '--figure-height', required=False, type=float, default=20, help='Figure height of plot.')
    partition_parser.add_argument('-s', '--size-points', required=False, type=float, default=0.7, help='Size of points in plot.')
    partition_parser.add_argument('-a', '--alpha', required=False, default=0.7, type=float, 
                                  help='Transparency of points in plot.')
    partition_parser.add_argument('-c', '--color-points', required=False, default='#FA8460', 
                                  help='Color of all points not in indexes.')
    partition_parser.add_argument('-ci', '--color-points-indexes', required=False, 
                                  default='#7793F5', help='Color of all points in indexes.')
    partition_parser.add_argument('-d', '--debug', required=False, help='Path to .log/.txt file to store info/debug messages.')
    partition_parser.add_argument('-l', '--legend-points', required=False, default="SNPs not in indexes", 
                                  help='Name given to all points not in indexes.')
    partition_parser.add_argument('-li', '--legend-points-indexes', required=False, default="SNPs in indexes", 
                                  help='Name given to all points in indexes.')
    
    # Define subparser for 'plot-pca' command
    partition_parser = subparsers.add_parser('plot-pca', help='To plot the PCA of the provided data.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for a single or multiple chromosomes each.')
    partition_parser.add_argument('-r', '--reference', required=False, nargs="*", 
                                  help='Paths to reference .vcf files with data for a single or multiple chromosomes each.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-t', '--train-query', action='store_true', 
                                  help='To train the PCA the PCA only on the query instead of on both datasets.')
    partition_parser.add_argument('-f', '--fontsize', required=False, default=35, type=float, 
                                  help='Fontsize of all text in plot.')
    partition_parser.add_argument('-w', '--figure-width', required=False, default=20, type=float, 
                                  help='Figure width of plot.')
    partition_parser.add_argument('-i', '--figure-height', required=False, default=20, type=float, 
                                  help='Figure height of plot.')
    partition_parser.add_argument('-s', '--size-points', required=False, default=20, type=float, 
                                  help='Size of points in plot.')
    partition_parser.add_argument('-a', '--alpha', required=False, default=0.7, type=float, 
                                  help='Transparency of points in plot.')
    partition_parser.add_argument('-cq', '--color-points-query', required=False, default='#EBD0A1', 
                                  help='Color of query points in plot.')
    partition_parser.add_argument('-cr', '--color-points-reference', required=False, default='#259988', 
                                  help='Color of reference points in plot.')
    partition_parser.add_argument('-d', '--debug', required=False, 
                                  help='Path to .log/.txt file to store info/debug messages.')
    partition_parser.add_argument('-lq', '--legend-query', required=False, default="Query", 
                                  help='Name given to all points in query.')
    partition_parser.add_argument('-lr', '--legend-reference', required=False, default="Reference", 
                              help='Name given to all points not in reference.')

    # Define subparser for 'remove-snps-different-means' command
    partition_parser = subparsers.add_parser('remove-snps-different-means', 
                                             help='To remove common markers with a mean absolute difference higher than a threshold.')
    partition_parser.add_argument('-q', '--query', required=True, nargs="*", 
                                  help='Paths to query .vcf files with data for a single chromosome each.')
    partition_parser.add_argument('-r', '--reference', required=False, nargs="*", 
                                  help='paths to reference .vcf files with data for a single chromosome each.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-t', '--threshold', required=False, type=float, default=0.1, 
                                  help='All common SNPs with a mean absolute difference higher than the threshold will be removed.')
    partition_parser.add_argument('-d', '--debug', required=False, 
                                  help='Path to .log/.txt file to store info/debug messages.')

    # Define subparser for 'store-common-indexes' command
    partition_parser = subparsers.add_parser('store-common-indexes', help='Store indexes of common markers.')
    partition_parser.add_argument('-q', '--query', required=True,
                                  help='Path to query .vcf file with data for a single or multiple chromosomes.')
    partition_parser.add_argument('-r', '--reference', required=True,
                                  help='Paths to reference .vcf file with data for a single or multiple chromosomes.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-f', '--file-format', required=False, choices=['.npy', '.h5'], default='.npy', 
                                  help='Format of the output file.')
    partition_parser.add_argument('-d', '--debug', required=False, 
                                  help='Path to .log/.txt file to store info/debug messages.')
    
    # Define subparser for 'store-allele' command
    partition_parser = subparsers.add_parser('store-npy', help='To store formatted allele data in .npy or .h5 format.')
    partition_parser.add_argument('-q', '--query', required=True, help='Path to query .vcf file.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-s', '--data-format', required=False, choices=['separated', 'averaged'], 
                                  default='separated', help='Separate or average maternal and paternal strands.')
    partition_parser.add_argument('-f', '--file-format', required=False, choices=['.npy', '.h5'], default='.npy', 
                                  help='Format of the output file.')
    partition_parser.add_argument('-d', '--debug', required=False, 
                                  help='Path to .log/.txt file to store info/debug messages.')
    
    # Define subparser for 'store-vcf' command
    partition_parser = subparsers.add_parser('store-vcf', help='To store separated allele data in .vcf.')
    partition_parser.add_argument('-q', '--vcf-query', required=True, help='Path to input .vcf file.')
    partition_parser.add_argument('-n', '--npy-query', required=True, help='Path to input .npy file.')
    partition_parser.add_argument('-o', '--output-folder', required=True, help='Path to output folder.')
    partition_parser.add_argument('-b', '--binary-indexes', required=False, 
                                  help='Path to binary indexes, where 1 means the SNP is correct and 0 that it is to be removed.')
    partition_parser.add_argument('-d', '--debug', required=False, 
                                  help='Path to .log/.txt file to store info/debug messages.')
    
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
    
    elif args.command == 'subset':
        # Check all the query and reference paths exist and are a .vcf file
        check_paths(args.query+args.reference)
        
        # Check amount of query and reference paths is the same
        assert len(args.query) == len(args.reference), f'The amount of query and reference paths does not coincide '\
        f'{len(args.query)} != {len(args.reference)}.'
        
    elif args.command in 'plot-snp-means':
        # Check all the query and reference paths exist and are a .vcf file
        check_paths(args.query+args.reference)
        
        # Check amount of query and reference paths is the same
        assert len(args.query) == len(args.reference), f'The amount of query and reference paths does not coincide '\
        f'{len(args.query)} != {len(args.reference)}.'
                
    elif args.command == 'plot-snp-means':
        # Check all the query paths exist and are a .vcf file
        check_paths(args.query)
        
        if args.reference == None:
            # Put reference in list to be compatible with other scripts
            args.reference = [None]
            
            # Ensure all tasks specified only work on the query
            if train_both is not None:
                raise argparse.ArgumentTypeError(f'--train_both is only supported when a '\
                                                 'reference is provided.')
        else:
            # Check all the reference paths exist and are a .vcf file
            check_paths(args.reference)
    
    elif args.command == 'remove-snps-different-means':
        # Check all the query and reference paths exist and are a .vcf file
        check_paths(args.query+args.reference)
        
        # Check amount of query and reference paths is the same
        assert len(args.query) == len(args.reference), f'The amount of query and reference paths does not coincide '\
        f'{len(args.query)} != {len(args.reference)}.'
        
        if args.threshold is not None:
        
            if args.threshold > 1.0 or args.threshold < 0.0:
                raise argparse.ArgumentTypeError(f'--threshold is out of the [0.0, 1.0] range.')
    
    elif args.command == 'store-common-indexes':
        # Check the query and reference paths exist
        check_paths([args.query, args.reference])
        
    elif args.command == 'store-npy':
        # Check the query path exists
        check_paths([args.query])
    
    elif args.command == 'store-vcf':
        # Check all the paths exist
        check_paths([args.vcf_query, args.npy_query, args.binary_indexes])
    

def check_chromosome(query: Dict, reference: Dict, single: bool = True):
    """
    By default, ensures the query and the reference contain data for a single
    chromosome and that such chromosome is the same in both datasets.
    Returns the chromosome/s in particular. If single = False, ensures the
    query and the reference contain the same ordered set of chromosomes. 
    
    
    Args:
        query (Dict): query data.
        reference (Dict): reference data.
        single (bool): to ensure data comes from a single chromosome.
    
    Returns:
        chrom (str): chromosome in the query and the reference.
    
    """
    
    # Obtain list with chromosomes in the query
    chroms_query = obtain_chromosomes(query)
    
    if reference is not None:
        # If the reference is not None...
        # Obtain list with chromosomes in the reference
        chroms_reference = obtain_chromosomes(reference)
        
        if single:
            # Ensure the query and the reference contain data for a single chromosome
            assert len(chroms_query) == 1, 'The query contains data for more than one chromosome. '\
            'Use partition command to obtain a separate VCF file per chromosome.'
            assert len(chroms_reference) == 1, 'The reference contains data for more than one chromosome. '\
            'Use partition command to obtain a separate VCF file per chromosome.'

        # Ensure the chromosomes is/are the same (in the same order)
        for chrom_query, chrom_reference in zip(chroms_query, chroms_reference):
            assert chrom_query == chrom_reference, f'The query contains data for chromosome {chroms_query[0]} '\
            f'but the reference contains data for chromosome {chroms_reference[0]}. Terminate.'
    
    if len(chroms_query) == 1: return chroms_query[0]
    else: return chroms_query


def define_plot_configuration(args: argparse.Namespace) -> Dict:
    """
    Defines plot hyperparams for plot-snp-means or plot-pca
    commands.
    
    Args:
        args (argparse.Namespace): parser with all arguments
        information.
    
    Returns:
        plot_dict (Dict): dictionary with plot configurations.
    
    """
    
    if args.command == 'plot-snp-means':
        # Define plot configurations for
        # 'plot-pca' commanad...
        plot_dict = {"x_axis_name": args.x_axis_name, 
                     "y_axis_name": args.y_axis_name, 
                     "fontsize" : args.fontsize, 
                     "fig_width": args.figure_width,
                     "fig_height": args.figure_height, 
                     "s": args.size_points,
                     "alpha" : args.alpha,
                     "color_points": args.color_points, 
                     "color_points_indexes": args.color_points_indexes,
                     "legend_points": args.legend_points, 
                     "legend_points_indexes": args.legend_points_indexes}
    
    elif args.command == 'plot-pca':
        # Define plot configurations for
        # 'plot-pca' commanad...
        plot_dict = {"fontsize": args.fontsize, 
                     "fig_width": args.figure_width,
                     "fig_height": args.figure_height, 
                     "s": args.size_points,
                     "alpha" : args.alpha,
                     "color_query": args.color_points_query,
                     "color_reference": args.color_points_reference,
                     "legend_query": args.legend_query,
                     "legend_reference": args.legend_reference}
    
    return plot_dict