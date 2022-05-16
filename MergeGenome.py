################################################################################
# Runs a specific command from MergeGenome to:
# * Preprocess .vcf files to properly merge genomic datasets.
# * Evalute the quality of merging.
################################################################################

import sys
import argparse

from utils.logger import get_logger
from utils.misc import define_parser, check_arguments, define_plot_configuration

from modules.partition import partition_by_chromosome
from modules.rename import rename_chromosome
from modules.clean import clean_genomic_data
from modules.subset import subset_common_markers

from evaluation.evaluation_plots import plot_snp_means, plot_pca

from scripts.store_allele_data import store_allele_data_in_npy_or_h5, store_allele_data_in_vcf
from scripts.search_indexes_common_snps import store_indexes_common_markers
from scripts.remove_snps_different_means import remove_snps_with_different_means

# Define parser
parser = define_parser()

# Parse the arguments
args = parser.parse_args()

# Check the input arguments are correct
check_arguments(args)

# Create logger
logger = get_logger(__name__, args.debug)

if args.command == 'partition':
    # Partition .vcf data in a separate .vcf file per chromosome and,
    # if rename-chr is set to True, rename chromosome nomenclature
    partition_by_chromosome(args.query, args.output_folder, args.rename_chr, args.rename_map, logger)
    
elif args.command == 'rename':    
    # Rename nomenclature in variants/CHROM field
    rename_chromosome(args.query, args.output_folder, args.rename_map, logger)
    
elif args.command == 'clean':
    # Clean genomic data
    clean_genomic_data(args.query, args.reference, args.output_folder, args.remove_sample_ID_query, 
                       args.remove_sample_ID_reference, args.remove_ambiguous_snps_query, 
                       args.remove_ambiguous_snps_reference, args.correct_snp_flips, 
                       args.remove_mismatching_snps, args.rename_map_query, 
                       args.rename_map_reference, logger)

elif args.command == 'subset':
    # Subset genomic data
    subset_common_markers(args.query, args.reference, args.output_folder, logger)
    
elif args.command == 'plot-snp-means':
    # Define plot configuration for this command
    plot_dict = define_plot_configuration(args)
    
    # Plot the SNP means for the query and the reference
    plot_snp_means(args.query, args.reference, plot_dict, args.output_folder, logger)
    
elif args.command == 'plot-pca':
    # Check the input arguments are correct
    #check_arguments(args.query)
    #if args.reference is not None:
    #    check_arguments(args.reference)
    
    # Define plot configuration for this command
    plot_dict = define_plot_configuration(args)
    
    # Plot the SNP means for the query and the reference
    plot_pca(args.query, args.reference, args.train_both, plot_dict, args.output_folder, logger)
    
elif args.command == 'store-npy':
    
    # Check the input arguments are correct
    #check_arguments([args.query])
    
    # Store averaged or separated maternal and paternal strands in .npy or .h5
    store_allele_data_in_npy_or_h5(args.query, args.output_folder, args.data_format, args.file_format, logger)
    
elif args.command == 'store-common-indexes':
    
    # Check the input arguments are correct
    #check_arguments([args.query]+[args.reference])
    
    # Store indexes of common markers between the query and the reference in .npy or .h5
    store_indexes_common_markers(args.query, args.reference, args.output_folder, args.file_format, logger)
    
elif args.command == 'remove-snps-different-means':
    
    # Check the input arguments are correct
    #check_arguments(args.query+args.reference)
    
    # Store indexes of common markers between the query and the reference in .npy or .h5
    remove_snps_with_different_means(args.query, args.reference, args.output_folder, args.threshold, logger)

elif args.command == 'store-vcf':
    
    # Store .npy file with separated maternal and paternal strands back in .vcf file 
    # and, optionally, remove undesired SNPs
    store_allele_data_in_vcf(args.vcf_query, args.npy_query, args.output_folder, args.binary_indexes, logger)
    