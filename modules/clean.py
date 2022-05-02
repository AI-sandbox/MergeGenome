################################################################################
# Function to clean genomic sequences:
# * Remove undesired samples by sample ID
# * Remove ambiguous SNPs
# * Correct SNP flips
# * Remove SNP mismatches
# * Rename missings (from key = old_value to value = new_value)
################################################################################

import os
import logging
from typing import List
from utils.io import read_vcf_file, write_vcf_file
from utils.vcf_utils import obtain_chromosomes, filter_samples, rename_missings
from utils.vcf_clean import search_and_remove_ambiguous_snps, search_and_correct_flips_by_pos, search_and_remove_mismatches_by_pos


def clean_genomic_data(query_paths: List[str], reference_paths: List[str], output_folder: str, 
                       remove_sample_ID_query: List[str], remove_sample_ID_reference: List[str], 
                       remove_ambiguous_snps_query: bool, remove_ambiguous_snps_reference: bool,
                       correct_snp_flips: bool, remove_mismatching_snps: bool, rename_map_query: str, 
                       rename_map_reference: str, logger: logging.Logger) -> None:
    """
    Cleans/preprocesses the genomic sequences from all reference and query files provided
    by applying all the preprocessing steps specified as input.
    
    Args:
        query_paths (List[str]): list with order paths to all query .vcf files to clean/preprocess.
        reference_paths (List[str]): list with order paths to all reference .vcf files to clean/preprocess.
        output_folder (str): folder to save the output.
        remove_sample_ID_query (List[str]): list with sample IDs or substring of samples IDs to remove.
        remove_sample_ID_reference (List[str]): list with sample IDs or substring of samples IDs to remove.
        remove_ambiguous_snps_query (bool): remove (or not) ambiguous SNPs from the query.
        remove_ambiguous_snps_reference (bool): remove (or not) ambiguous SNPs from the reference.
        correct_snp_flips (bool): correct (or not) SNP in the query with respect to the reference.
        remove_mismatching_snps (bool): remove (or not) mismatching SNPs between the reference and the query.
        rename_map_query (str): dictionary with mapping from old to new missing notation for the query.
        rename_map_reference (str): dictionary with mapping from old to new missing notation for the reference.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """
    
    if reference_paths is None:
        # To be able to iterate over empty reference paths
        reference_paths = [None]*len(query_paths)
    
    for query_path, reference_path in zip(query_paths, reference_paths):
        
        # Read the query .vcf file using scikit-allel
        logger.debug('Reading the query file.')
        query = read_vcf_file(query_path)
        logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples in total in the query.')
        
        if reference_path is not None:
            logger.debug('Reading the reference file.')
            reference = read_vcf_file(reference_path)
            logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples in total in the reference.')
        
            # Ensure the reference and the query contain data for the same chromosome
            chroms_r = obtain_chromosomes(reference)
            chroms_q = obtain_chromosomes(query)
            assert len(chroms_r) == len(chroms_q) == 1, 'The reference or the query contain data for more than one chromosome.'\
            'Use partition command to partition the data into a separate VCF file per chromosome.'
            assert chroms_r[0] == chroms_q[0], 'The reference and the query contain data for a different chromosome. Check the order of the inputs.'
            logger.info(f'Cleaning chromosome {chroms_r[0]}...')
        
        if remove_sample_ID_query is not None:
            
            # Subset samples to remove undesired breeds/species with name in substrings list in the query
            logger.debug(f'Removing samples by ID in the query dataset.')
            query = filter_samples(query, remove_sample_ID_query, logger)
            logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples ' \
                        f'in total in the query after removing undesired samples.')
            
        if remove_sample_ID_reference is not None:
        
            # Subset samples to remove undesired breeds/species with name in substrings list in the reference
            logger.debug(f'Removing samples by ID in the reference dataset.')
            reference = filter_samples(reference, remove_sample_ID_reference, logger)
            logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples ' \
                        f'in total in the reference after removing undesired samples.')
            
        if remove_ambiguous_snps_query:
            
            # Search and remove ambiguous SNPs from the query
            logger.debug(f'Searching/removing ambiguous SNPs in the query dataset.')
            query = search_and_remove_ambiguous_snps(query, logger)
            logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples ' \
                        f'in total in the query after removing ambiguous SNPs.')
            
        if remove_ambiguous_snps_reference:
            
            # Search and remove ambiguous SNPs from the reference
            logger.debug(f'Searching/removing ambiguous SNPs in the reference dataset.')
            reference = search_and_remove_ambiguous_snps(reference, logger)
            logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples ' \
                        f'in total in the reference after removing ambiguous SNPs.')
            
        if correct_snp_flips:
            
            # Correct SNP flips in the query with respect to the reference
            logger.debug('Searching/correcting SNP flips in the reference with respect to the query')
            query = search_and_correct_flips_by_pos(reference, query, logger)
        
        if remove_mismatching_snps:
            
            # Remove mismtaching SNPs between the reference and the query
            logger.debug('Searching/removing mismatching SNPs between the reference and the query')
            reference, query = search_and_remove_mismatches_by_pos(reference, query, logger)
            
        if rename_map_query is not None:
            
            # Obtain dict from str
            rename_map_query = eval(rename_map_query)
        
            # Rename missing values in the query
            logger.debug(f'Renaming missing values in the query from {rename_map_query.keys()} to {rename_map_query.values()}.')
            query = rename_missings(query, rename_map_query.keys(), rename_map_query.values())
        
        if rename_map_reference is not None:
            
            # Obtain dict from str
            rename_map_reference = eval(rename_map_reference)
            
            # Rename missing values in the reference
            logger.debug(f'Renaming missing values in the reference from {rename_map_reference.keys()} to {rename_map_reference.values()}.')
            reference = rename_missings(reference, rename_map_reference.keys(), rename_map_reference.values())
        
        # Define output name to cleaned reference and query .vcf files
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_cleaned.vcf'
        output_name_query = f'{os.path.basename(query_path)[:-4]}_cleaned.vcf'
        
        # Write cleaned reference and query in output .vcf files
        logger.debug(f'Writing cleaned VCF data for the reference.')
        # write_vcf_file(reference, output_folder, output_name_reference)
        
        logger.debug(f'Writing cleaned VCF data for the query.')
        # write_vcf_file(query, output_folder, output_name_query)
        