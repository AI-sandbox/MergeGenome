################################################################################
# Convert VCF file to .npy or .h5.
################################################################################

import os
import logging
from typing import List
import h5py
import numpy as np
from utils.io import read_vcf_file, read_npy_file, write_vcf_file
from utils.vcf_utils import combine_chrom_strands


def store_allele_data_in_npy_or_h5(query_path: str, output_folder: str, data_format: str, 
                                   file_format: str, logger: logging.Logger) -> None:
    
    """
    Separates or averages the maternal and paternal strands and 
    stores the allel data in calldata/GT in .npy or .h5 format.
    
    Args:
        query_path (str): path to query .vcf file.
        output_folder (str): path to output folder.
        data_format (str): specify the format of the data for the analysis. 
            * 'separted': split maternal and paternal strands into separate samples. 
            * 'averaged': combine maternal and  paternal strands by averaging them.
        file_format (str): format of the file to store the allele data.
            * '.npy': store data in numpy format.
            * '.h5': store data in h5py format.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        (None)
    
    """
    
    # Read the .vcf file using scikit-allel
    logger.debug('Reading the provided genomic file.')
    data = read_vcf_file(query_path, logger)    
    
    # Obtain SNPs data
    snps = data['calldata/GT']
        
    # Obtain the dimensions of the data
    num_snps, num_samples, num_strands = snps.shape
    
    # Convert SNPs from MxNx2 to Mx2N shape
    logger.debug('Separating maternal and paternal strands.')
    snps = snps.reshape(num_snps, num_samples*2).T
    num_samples, num_snps = snps.shape
    
    if data_format == 'averaged':
        
        # Convert SNPs from Mx2N shape to MxN shape
        # Average the maternal and paternal strands
        logger.debug('Averaging the maternal and paternal strands.')
        snps = combine_chrom_strands(snps)
        
    logger.info(f'There are {num_snps} SNPs and {num_samples} samples in total after "{data_format}" reshaping.')
    
    if file_format == '.npy':
        
        # Define output name to .npy file with .npy SNPs data
        # The .npy file has the same base name as the input file, but 
        # ending with '_separated.npy' or '_averaged.npy'
        output_name = f'{os.path.basename(query_path)[:-4]}_{data_format}.npy'
        
        # Define entire output path
        output_path = output_folder + output_name
        
        # Save the snps data in .npy format
        logger.debug(f'Writing data in .npy format in {output_path}.')
        np.save(output_path, snps)

    elif file_format == '.h5':

        # Define output name to .h5 file with .h5 SNPs data
        # The .h5 file has the same base name as the input file, but 
        # ending with '_separated.h5' or '_averaged.h5'
        output_name = f'{os.path.basename(query_path)[:-4]}_{data_format}.h5'
        
        # Define entire output path
        output_path = output_folder + output_name
        
        # Save the snps data in .h5 format
        logger.debug(f'Writing data in .h5 format in {output_path}.')
        h5f = h5py.File(output_path, 'w')
        h5f.create_dataset(name=output_path, data=snps)
        
        
def store_allele_data_in_vcf(vcf_path: str, npy_path: str, output_folder: str, binary_indexes_path: np.array,
                             logger: logging.Logger) -> None:
    
    """
    Reads .npy with separated maternal and paternal strands and 
    stores the allel data in calldata/GT from a given VCF file.
    If --binary_indexes, removes all SNPs with value zero in the provided
    numpy array. Stores the result in a VCF file.
    
    Args:
        vcf_path (str): path to .vcf file.
        npy_path (str): path to .npy file.
        output_folder (str): folder to save the output.
        binary_indexes (np.array): path to binary indexes where 1 indicates the 
        SNP is correct and 0 otherwise. Incorrects SNPs will be removed from output
        VCF file.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        (None)
    
    """
    
    # Read the .vcf file using scikit-allel
    logger.debug('Reading the provided genomic file.')
    data = read_vcf_file(vcf_path, logger)
    
    # Obtain SNPs data
    snps = data['calldata/GT']
    
    # Obtain the dimensions of the data
    num_snps, num_samples, num_strands = snps.shape
    
    # Read the .npy data with Mx2N shape
    snps_flattened = read_npy_file(npy_path, logger)
    
    # Obtain array of all ones that will be filled with .npy data
    new_snps = np.ones((num_snps, num_samples, 2))
    
    # Store .npy allele data in calldata/GT VCF allele data
    # Even strands go to position 0 in the 3rd dimension and
    # odd strands go to position 1 in the 3rd dimension
    snps[:,:,0] = snps_flattened[::2,:].T
    snps[:,:,1] = snps_flattened[1::2,:].T
    
    # Read the .npy with the binary indexes
    binary_indexes = read_npy_file(binary_indexes_path, logger)
    
    if binary_indexes is not None:
        # Convert 1 to True and 0 to False
        binary_indexes = binary_indexes.astype(bool)
        
        # Keep correct SNPs and remove incorrect ones
        data = select_snps(data, binary_indexes)
        
    # Define output name to .vcf file with .npy allele data
    # The .vcf file has the same base name as the input file, but 
    # ending with '_from_npy.vcf'
    output_name = f'{os.path.basename(vcf_path)[:-4]}_from_npy.vcf'
    
    # Write query data with .npy allele data in .vcf file
    # with name output_name inside folder output_folder
    logger.debug(f'Writing .vcf data with .npy allele data in {output_folder}{output_name}.')
    write_vcf_file(data, output_folder, output_name)

