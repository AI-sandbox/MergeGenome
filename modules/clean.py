################################################################################
# Cleans genomic sequences by:
# * Removing undesired samples by sample ID.
# * Removing ambiguous SNPs.
# * Correcting SNP flips.
# * Removing SNP mismatches.
# * Renaming missing values (from key = old_value to value = new_value 
# in mapping).
################################################################################

import os
import logging
from typing import List
from itertools import zip_longest

from utils.io import read_vcf_file, write_vcf_file
from utils.misc import check_chromosome
from utils.vcf_utils import filter_samples, rename_missings
from utils.vcf_clean import remove_ambiguous_snps, correct_flips_by_pos, remove_mismatches_by_pos


def clean_genomic_data(query_paths: List[str], reference_paths: List[str], output_folder: str, 
                       remove_sample_ID_query: List[str], remove_sample_ID_reference: List[str], 
                       remove_ambiguous_snps_query: bool, remove_ambiguous_snps_reference: bool,
                       correct_snp_flips: bool, remove_mismatching_snps: bool, rename_map_query: str, 
                       rename_map_reference: str, logger: logging.Logger) -> None:
    """
    Cleans the genomic sequences from the reference and query files provided
    by applying the preprocessing steps specified.
    
    Args:
        query_paths (List[str]): paths to query .vcf files with data for a single chromosome each.
        reference_paths (List[str]): paths to reference .vcf files with data for a single chromosome 
        each.
        output_folder (str): path to output folder.
        remove_sample_ID_query (List[str]): sample IDs or substring of samples IDs to remove from 
        the query.
        remove_sample_ID_reference (List[str]): sample IDs or substring of samples IDs to remove from
        the reference.
        remove_ambiguous_snps_query (bool): to remove (or not) ambiguous SNPs from the query.
        remove_ambiguous_snps_reference (bool): to remove (or not) ambiguous SNPs from the reference.
        correct_snp_flips (bool): to correct (or not) SNP flips in the query with respect to the 
        reference.
        remove_mismatching_snps (bool): to remove (or not) mismatching SNPs between the reference and
        the query.
        rename_map_query (str): mapping from old to new missing values notation for the query.
        rename_map_reference (str): mapping from old to new missing value notation for the reference.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """
    
    for query_path, reference_path in zip_longest(query_paths, reference_paths, fillvalue='None'):
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
        
        if remove_sample_ID_query is not None:
            # Filter undesired samples in remove_sample_ID_query from the query
            logger.debug(f'Removing samples by ID from the query.')
            query = filter_samples(query, remove_sample_ID_query, logger)
            logger.info(f'There are {len(query["variants/ID"])} SNPs and '\
                        f'{len(query["samples"])} samples in total in the query '\
                        f'after removing undesired samples.')
        
        if remove_sample_ID_reference is not None:
            # Filter undesired samples in remove_sample_ID_reference from the reference
            logger.debug(f'Removing samples by ID from the reference.')
            reference = filter_samples(reference, remove_sample_ID_reference, logger)
            logger.info(f'There are {len(reference["variants/ID"])} SNPs and '\
                        f'{len(reference["samples"])} samples in total in the reference after'\
                        f'removing undesired samples.')
        
        if remove_ambiguous_snps_query:
            # Search and remove ambiguous SNPs from the query
            logger.debug(f'Searching and removing ambiguous SNPs in the query.')
            query = remove_ambiguous_snps(query, logger)
            logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])}'\
                        f'samples in total in the query after removing ambiguous SNPs.')
        
        if remove_ambiguous_snps_reference:
            # Search and remove ambiguous SNPs from the reference
            logger.debug(f'Searching and removing ambiguous SNPs in the reference.')
            reference = remove_ambiguous_snps(reference, logger)
            logger.info(f'There are {len(reference["variants/ID"])} SNPs and '
                        f'{len(reference["samples"])} samples in total in the reference after '\
                        f'removing ambiguous SNPs.')
        
        if correct_snp_flips:
            # Correct SNP flips in the query with respect to the reference
            logger.debug('Searching and correcting SNP flips in the query with respect to '\
                         'the reference.')
            query = correct_flips_by_pos(reference, query, logger)
        
        if remove_mismatching_snps:
            # Remove mismtaching SNPs between the reference and the query
            logger.debug('Searching and removing mismatching SNPs between the reference and the '\
                         'query.')
            reference, query = remove_mismatches_by_pos(reference, query, logger)
            
        if rename_map_query is not None:
            # Rename missing values in the query
            logger.debug(f'Renaming missing values in the query from {rename_map_query.keys()} '\
                         f'to {rename_map_query.values()}.')
            query = rename_missings(query, rename_map_query.keys(), rename_map_query.values())
        
        if rename_map_reference is not None:
            # Rename missing values in the reference
            logger.debug(f'Renaming missing values in the reference from '\
                         f'{rename_map_reference.keys()} to {rename_map_reference.values()}.')
            reference = rename_missings(reference, rename_map_reference.keys(),
                                        rename_map_reference.values())
        
        # Define output name to .vcf file with cleaned data for the 
        # query and to .vcf file with cleaned data for the reference
        # The .vcf files have the same base name as the input file, but 
        # ending with '_clenaed.vcf'
        # Example: query_chr1_clenaed.vcf
        output_name_query = f'{os.path.basename(query_path)[:-4]}_cleaned.vcf'
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_cleaned.vcf'
        
        logger.debug(f'Writing .vcf data for the query in {output_folder}{output_name_query}.')
        # write_vcf_file(query, output_folder, output_name_query)
        
        # Write cleaned reference and query in output .vcf files
        logger.debug(f'Writing .vcf data for the reference in {output_folder}{output_name_reference}.')
        # write_vcf_file(reference, output_folder, output_name_reference)
        