################################################################################
# Subsets the reference and query datasets 
# to the common markers
################################################################################

import os
import logging
from typing import List

from utils.io import read_vcf_file, write_vcf_file
from utils.vcf_utils import obtain_chromosomes
from utils.vcf_clean import search_and_keep_common_markers_single_chr


def subset_common_markers(reference_paths: List[str], query_paths: List[str], output_folder: str, logger: logging.Logger) -> None:

    """
    Subsets the reference and the query to the common markeys (i.e.
    SNPs at the same CHROM, POS, REF and ALT).
    
    Args:
        reference_paths (List[str]): list with order paths to all reference .vcf files to clean/preprocess.
        query_paths (List[str]): list with order paths to all query .vcf files to clean/preprocess.
        output_folder (str): folder to save the output.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """
    
    for reference_path, query_path in zip(reference_paths, query_paths):
        
        # Read the reference and query .vcf files using scikit-allel
        logger.debug('Reading the reference and query files.')
        reference = read_vcf_file(reference_path)
        query = read_vcf_file(query_path)
        
        # Obtain the dimensions of the reference and query
        logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples in total in the reference.')
        logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples in total in the query.')
        
        # Ensure the reference and the query contain data for the same chromosome
        chroms_r = obtain_chromosomes(reference)
        chroms_q = obtain_chromosomes(query)
        assert len(chroms_r) == len(chroms_q) == 1, 'The reference or the query contain data for more than one chromosome.'\
        'Use partition command to partition the data into a separate VCF file per chromosome.'
        assert chroms_r[0] == chroms_q[0], 'The reference and the query contain data for a different chromosome. Check the order of the inputs.'
        logger.debug(f'Cleaning chromosome {chroms_r[0]}...')
        
        ## Keep common SNPs of the first and the second datasets, if present in the third dataset
        # The SNPs that are in the first or second dataset but not in the third dataset are removed
        reference, query, _, _ = search_and_keep_common_markers_single_chr(reference, query, logger)
        
        logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples in total in the reference '\
                    'after subtetting to common markers.')
        logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples in total in the query '\
                    'after subtetting to common markers.')
        
        # Define output name to cleaned reference and query .vcf files
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_cleaned.vcf'
        output_name_query = f'{os.path.basename(query_path)[:-4]}_cleaned.vcf'
        
        # Write cleaned reference and query in output .vcf files
        logger.debug(f'Writing subsetted VCF data for the reference.')
        # write_vcf_file(reference, output_folder, output_name_reference)
        
        logger.debug(f'Writing subsetted VCF data for the query.')
        # write_vcf_file(query, output_folder, output_name_query)