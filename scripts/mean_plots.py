################################################################################
# Plot the SNP means of the common markers between the reference and 
# the query datasets
################################################################################

import os
import numpy as np
import pandas as pd
from typing import Dict
import matplotlib.pyplot as plt

from utils.vcf_utils import read_vcf_file, combine_chrom_strands
from utils.vcf_clean import search_and_keep_common_markers_single_chr
from utils.plot_utils import snp_means_plot


def plot_snp_means(query_paths: List[str], reference_paths: List[str], plot_dict: Dict, output_folder: str, logger: logging.Logger) -> None:
    """
    Plots the SNP means of all the common markers between the reference and 
    the query datasets. If data for more than one chromosome is provided,
    it will be concatenated and displayed in a single plot. The mean of 
    each SNP for each dataset is computed by first averaging the maternal
    and the paternal strands.
    
    Args:
        query_paths (List[str]): list with order paths to all query .vcf files to clean/preprocess.
        reference_paths (List[str]): list with order paths to all reference .vcf files to clean/preprocess.
        plot_dict[Dict]: dictionary with the configuration parameters of the plot.
        output_folder (str): folder to save the output.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """

    for query_path, reference_path in zip(query_paths, reference_paths):
        
        # Read the query .vcf file using scikit-allel
        logger.debug('Reading the query file.')
        query = read_vcf_file(query_path)
        logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples in total in the query.')
        
        logger.debug('Reading the reference file.')
        reference = read_vcf_file(reference_path)
        logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples in total in the reference.')
        
        # Ensure the reference and the query contain data for the same chromosome
        chroms_r = obtain_chromosomes(reference)
        chroms_q = obtain_chromosomes(query)
        assert len(chroms_r) == len(chroms_q) == 1, 'The reference or the query contain data for more than one chromosome.'\
        'Use partition command to partition the data into a separate VCF file per chromosome.'
        assert chroms_r[0] == chroms_q[0], 'The reference and the query contain data for a different chromosome. Check the order of the inputs.'
        logger.debug(f'Cleaning chromosome {chroms_r[0]}...')
        
        # Common markers
        reference, query, _, _ = search_and_keep_common_markers_single_chr(reference, query, logger)
        
        # Obtain snps for a specific chromosome
        query_snps = query['calldata/GT']
        reference_snps = reference['calldata/GT']
        
        # Convert snps from MxNx2 to MxN shape (averaged)
        length_query, num_dogs_query, num_strands_query = query_snps.shape
        query_snps = query_snps.reshape(length_query, num_dogs_query*2).T
        query_snps = combine_chrom_strands(query_snps)

        length_reference, num_dogs_reference, num_strands_reference = reference_snps.shape
        reference_snps = reference_snps.reshape(length_reference, num_dogs_reference*2).T
        reference_snps = combine_chrom_strands(reference_snps)
        
        # Concatenate snps for all chromosomes
        if i == 1:
            all_query_snps = query_snps
            all_reference_snps = reference_snps
        else:
            all_query_snps = np.concatenate((all_query_snps, query_snps), axis=1)
            all_reference_snps = np.concatenate((all_reference_snps, reference_snps), axis=1)
        
        logger.info(f'{all_query_snps.shape[1]} common SNPs in total between the query and the reference datasets.')
    
    # Plot the SNP means for the SNPs that are common markers between all_query_snps and all_reference_snps
    snp_means_plot(all_query_snps, all_reference_snps, plot_dict, output_folder)
