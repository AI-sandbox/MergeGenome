################################################################################
# Searches for the indexes of the SNPs of a query
# dataset that are also present in a reference dataset and viceversa
# The indexes of the common markers are saved in a .npy or .h5 file
################################################################################

import os
import logging
from typing import List
import h5py
import numpy as np

from utils.io import read_vcf_file
from utils.vcf_clean import search_and_keep_common_markers_several_chr


def store_indexes_common_markers(query_path: str, reference_path: str, output_folder: str, file_format: str) -> None:
    """
    Searches for the indexes of the SNPs of a query dataset that are also 
    present in a reference dataset and viceversa. The indexes of the
    common markers are saved in a .npy or .h5 file.
    
    Args:
        query_path (str): path to query .vcf file.
        reference_path (str): path to reference .vcf file.
        output_folder (str): folder to save the output.
        file_format (str): format of the file to store the indexes.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """
    
    # Read the reference and query .vcf files using scikit-allel
    logger.debug('Reading the reference and query files.')
    reference = read_vcf_file(reference_path)
    query = read_vcf_file(query_path)

    # Obtain the dimensions of the reference and query
    logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples in total in the reference.')
    logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples in total in the query.')

    # Obtain indexes of the common markers between the query and the referece
    _, _, indexes_query, indexes_reference = search_and_keep_common_markers_several_chr(query, reference, logger)

    # Define counter of common markers to ensure the found indexes 
    # are common markers
    count = 0
    for i in range (len(indexes_query)):
        if(query['variants/REF'][indexes_query[i]] != reference['variants/REF'][indexes_reference[i]] or 
           query['variants/ALT'][indexes_query[i]][0] != reference['variants/ALT'][indexes_reference[i]][0]):
            count += 1
    
    # Ensure the found indexes are common markers
    assert count == len(indexes_query) == len(indexes_reference)
    
    if file_format == '.npy':
        
        # Define output name to .npy file with query and reference indexes
        output_name_query = f'{os.path.basename(query_path)[:-4]}_indexes.npy'
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_indexes.npy'
        
        # Define entire output path to query and reference indexes
        output_path_query = output_folder + output_name_query
        output_path_reference = output_folder + output_name_reference
        
        # Save the snps data in .npy format
        logger.debug(f'Saving data in .npy format in {output_path_query}.')
        np.save(output_path_query, indexes_query)
        logger.debug(f'Saving data in .npy format in {output_path_reference}.')
        np.save(output_path_reference, indexes_reference)

    elif file_format == '.h5':

        # Define output name to .npy file with query and reference indexes
        output_name_query = f'{os.path.basename(query_path)[:-4]}_indexes.h5'
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_indexes.h5'
        
        # Define entire output path to query and reference indexes
        output_path_query = output_folder + output_name_query
        output_path_reference = output_folder + output_name_reference
        
        # Save the snps data in .h5 format
        logger.debug(f'Saving data in .h5 format in {output_path_query}.')
        h5f = h5py.File(output_path_query, 'w')
        h5f.create_dataset(name=output_path_query, data=indexes_query)
        logger.debug(f'Saving data in .h5 format in {output_path_reference}.')
        h5f = h5py.File(output_path_reference, 'w')
        h5f.create_dataset(name=output_path_reference, data=indexes_reference)
