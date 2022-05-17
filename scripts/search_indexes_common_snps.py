################################################################################
# Searches for the indexes of the common markers between a query and 
# reference datasets.
# The indexes of the common markers are saved in a .npy or .h5 file.
################################################################################

import os
import logging
from typing import List
import h5py
import numpy as np

from utils.io import read_vcf_file
from utils.vcf_clean import keep_common_markers_several_chr


def store_indexes_common_markers(query_path: str, reference_path: str, output_folder: str, 
                                 file_format: str, logger: logging.Logger) -> None:
    """
    Searches the indexes of the common markers (i.e., SNPs at the same CHROM, 
    POS, REF, and ALT) between the query and the reference. The indexes of the
    common markers are saved in a .npy or .h5 file.
    
    Args:
        query_path (str): path to query .vcf file with data for a single or 
        multiple chromosomes.
        reference_path (str): path to reference .vcf file with data for a 
        single or multiple chromosomes.
        output_folder (str): folder to save the output.
        file_format (str): format of the file to store the indexes.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """
    
    # Read the query .vcf file using scikit-allel
    # with data from a single or multiple chromosomes
    logger.debug(f'Reading query file {query_path}.')
    query = read_vcf_file(query_path, logger)
    
    # Read the reference .vcf file using scikit-allel
    # with data from a single or multiple chromosomes
    logger.debug(f'Reading reference file {reference_path}.')
    reference = read_vcf_file(reference_path, logger)
    
    # Obtain indexes of the common markers between the query and the referece
    _, _, idxs_reference, idxs_query = keep_common_markers_several_chr(query, reference, logger)
    
    # Ensure the found indexes have the same length in the 
    # query than in the reference
    assert len(idxs_query) == len(idxs_reference)
    
    if file_format == '.npy':
        
        # Define output name to .npy file with query and reference indexes
        output_name_query = f'{os.path.basename(query_path)[:-4]}_indexes.npy'
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_indexes.npy'
        
        # Define entire output path to .npy files with query and reference indexes
        # The .npy files have the same base name as the input file, but 
        # ending with '_indexes.npy'
        # Example: query_chr1_indexes.npy
        output_path_query = output_folder + output_name_query
        output_path_reference = output_folder + output_name_reference
        
        # Save the snps data in .npy format
        logger.debug(f'Saving data in .npy format in {output_path_query}.')
        np.save(output_path_query, idxs_query)
        logger.debug(f'Saving data in .npy format in {output_path_reference}.')
        np.save(output_path_reference, idxs_reference)

    elif file_format == '.h5':

        # Define entire output path to .h5 files with query and reference indexes
        # The .h5 files have the same base name as the input file, but 
        # ending with '_indexes.h5'
        # Example: reference_chr1_indexes.h5
        output_name_query = f'{os.path.basename(query_path)[:-4]}_indexes.h5'
        output_name_reference = f'{os.path.basename(reference_path)[:-4]}_indexes.h5'
        
        # Define entire output path to query and reference indexes
        output_path_query = output_folder + output_name_query
        output_path_reference = output_folder + output_name_reference
        
        # Save the snps data in .h5 format
        logger.debug(f'Saving data in .h5 format in {output_path_query}.')
        h5f = h5py.File(output_path_query, 'w')
        h5f.create_dataset(name=output_path_query, data=idxs_query)
        logger.debug(f'Saving data in .h5 format in {output_path_reference}.')
        h5f = h5py.File(output_path_reference, 'w')
        h5f.create_dataset(name=output_path_reference, data=idxs_reference)
