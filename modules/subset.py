################################################################################
# Subsets the reference and query datasets 
# to the common markers (i.e. SNPs with identical CHROM, POS, REF, 
# and ALT fields) between both datasets with data for a particular chromosome.
################################################################################

import os
import logging
from typing import List

from utils.io import read_vcf_file, write_vcf_file
from utils.misc import check_chromosome
from utils.vcf_clean import keep_common_markers_single_chr


def subset_common_markers(query_paths: List[str], reference_paths: List[str], output_folder: str, 
                          logger: logging.Logger) -> None:

    """
    Filters the reference and the query to only contain the common markers (i.e.
    SNPs at the same CHROM, POS, REF and ALT) between both datasets.
    
    Args:
        query_paths (List[str]): paths to reference .vcf files with data for a single chromosome each.
        reference_paths (List[str]): paths to query .vcf files with data for a single chromosome each.
        output_folder (str): path to output folder.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """
    
    for query_path, reference_path in zip(query_paths, reference_paths):
        # For each query and reference paths...
        # Read the query .vcf file using scikit-allel
        # with data from a single chromosome
        logger.debug(f'Reading query file {query_path}.')
        query = read_vcf_file(query_path, logger)
        
        # Read the reference .vcf file using scikit-allel
        # with data from a single chromosome
        logger.debug(f'Reading reference file {reference_path}.')
        reference = read_vcf_file(reference_path, logger)
        
        # Ensure the reference and the query contain data for the same chromosome
        # and obtain the chromosome in particular
        chrom = check_chromosome(query, reference)
        logger.debug(f'Searching common markers in chromosome {chrom}...')
    
        # Filter the reference and the query to only contain the common markers
        reference, query, _, _ = keep_common_markers_single_chr(reference, query, logger)
        
        logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} '\
                    'samples in the query after filtering to common markers.')
        logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} '\
                    'samples in the reference after filtering to common markers.')
        
        # Define output name to .vcf file with filtered SNPs for the 
        # query and to .vcf file with filtered SNPs for the reference
        # The .vcf files have the same base name as the input file, but 
        # ending with '_subsetted.vcf'
        # Example: query_chr1_subsetted.vcf
        output_name_query = f'{os.path.basename(query_path)[:-4]}_subsetted.vcf'
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_subsetted.vcf'
        
        logger.debug(f'Writing .vcf data for the query in {output_folder}{output_name_query}.')
        # write_vcf_file(query, output_folder, output_name_query)
        
        # Write cleaned reference and query in output .vcf files
        logger.debug(f'Writing .vcf data for the reference in {output_folder}{output_name_reference}.')
        # write_vcf_file(reference, output_folder, output_name_reference)