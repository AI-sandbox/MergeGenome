################################################################################
# Removes all the common SNPs with a mean absolute difference higher than 
# a specified threshold and stores the result in a new VCF file
################################################################################

import os
import logging
from typing import List

from utils.io import read_vcf_file, write_vcf_file
from utils.misc import check_chromosome
from utils.vcf_clean import remove_snps_different_means


def remove_snps_with_different_means(query_paths: List[str], reference_paths: List[str], output_folder: str, 
                                     threshold: float, logger: logging.Logger) -> None:
    """
    To remove common markers (i.e., SNPs at the same CHROM, POS, REF, 
    and ALT) with a mean absolute difference higher than a threshold. 
    The filtered query and reference are stored in new VCF files.
    
    Args:
        query_paths (List[str]): paths to query .vcf files with data for a 
        single chromosome each.
        reference_paths (List[str]): paths to reference .vcf files with data 
        for a single chromosome each.
        output_folder (str): folder to save the output.
        threshold (float): all common SNPs with a mean absolute difference 
        higher than the threshold will be removed.
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
        # If reference_path = None, reference will also be None
        logger.debug(f'Reading reference file {reference_path}.')
        reference = read_vcf_file(reference_path, logger)
        
        # Ensure the reference and the query contain data for the same chromosome
        # and obtain the chromosome in particular
        chrom = check_chromosome(query, reference)
        logger.debug(f'Cleaning chromosome {chrom}...')
        
        # Remove all SNPs with a mean absolute difference higher than a threshold
        logger.debug(f'Removing all SNPs with a mean absolute difference higher than {threshold}')
        reference, query = remove_snps_different_means(reference, query, threshold, logger)

        logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples '\
                    f'in total in the query after removing SNPs with a mean abs difference > {threshold}.')
        logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples '\
                    f'in total in the reference after removing SNPs with a mean abs difference > {threshold}.')
        
        # Define output name to .vcf file with removed SNPs for the 
        # query and to .vcf file with removed SNPs for the reference
        # The .vcf files have the same base name as the input file, but 
        # ending with '_removed.vcf'
        # Example: query_chr1_removed.vcf
        output_name_query = f'{os.path.basename(query_path)[:-4]}_removed.vcf'
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_removed.vcf'
        
        logger.debug(f'Writing .vcf data for the query in {output_folder}{output_name_query}.')
        # write_vcf_file(query, output_folder, output_name_query)
        
        # Write cleaned reference and query in output .vcf files
        logger.debug(f'Writing .vcf data for the reference in {output_folder}{output_name_reference}.')
        # write_vcf_file(reference, output_folder, output_name_reference)
