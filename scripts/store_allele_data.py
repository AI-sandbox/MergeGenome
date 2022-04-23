################################################################################
# Convert VCF file to .npy
################################################################################


import os
import logging
from typing import List
import h5py
import numpy as np
from utils.io import read_vcf_file
from utils.vcf_utils import combine_chrom_strands


def store_allele_data_in_npy_or_h5(input_path: str, output_folder: str, data_format: str, file_format: str, logger: logging.Logger) -> None:
    
    """
    Separates or averages the maternal and paternal strands and 
    stores the allel data in calldata/GT in .npy or .h5 format.
    
    Args:
        input_path (str): path to .vcf file.
        output_folder (str): folder to save the output.
        data_format (str): specify the format of the data for the analysis. 
            * 'separted': split maternal and paternal strands into separate samples. 
            * 'averaged': combine maternal and  paternal strands by averaging them.
        file_format: format of the file to store the allele data.
            * '.npy': store data in numpy format.
            * '.h5': store data in h5py format.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        (None)
    
    """
    
    # Read the .vcf file using scikit-allel
    logger.debug('Reading the provided genomic file.')
    data = read_vcf_file(input_path)
        
    # Obtain SNPs data
    snps = data['calldata/GT']
        
    # Obtain the dimensions of the data
    logger.debug('Separating maternal and paternal strands.')
    num_snps, num_samples, num_strands = snps.shape
    logger.info(f'There are {num_snps} SNPs and {num_samples} samples in total.')
    
    # Convert SNPs from MxNx2 to Mx2N shape
    snps = snps.reshape(num_snps, num_samples*2).T
    num_samples, num_snps = snps.shape
    
    if data_format == 'averaged':
        
        # Convert SNPs from Mx2N shape to MxN shape
        # Average the maternal and paternal strands
        logger.debug('Averaging the maternal and paternal strands.')
        snps = combine_chrom_strands(snps)
        
    logger.info(f'There are {num_snps} SNPs and {num_samples} samples in total after "{data_format}" reshaping.')
    
    if file_format == '.npy':
        
        # Define output name to .npy file
        output_name = f'{os.path.basename(input_path)[:-4]}_{data_format}.npy'
        
        # Save the snps data in .npy format
        logger.debug('Saving data in .npy format.')
        np.save(output_folder + output_name, snps)

    elif file_format == '.h5':

        # Define output name to .npy file
        output_name = f'{os.path.basename(input_path)[:-4]}_{data_format}.h5'
        
        # Save the snps data in .h5 format
        logger.debug('Saving data in .h5 format.')
        h5f = h5py.File(output_path + output_file, 'w')
        h5f.create_dataset(name=output_folder + output_name, data=snps)
    
    
    
    
    
    
    
    
    

