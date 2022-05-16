################################################################################
# Makes plots to evaluate the quality of the merged dataset or identify outlers:
# * Plots the SNP means of the common markers between the reference and 
# the query datasets.
# * Plots the PCA of the common markers between the reference and 
# the query datasets.
# * Plots the PCA of the query.
################################################################################

import os
import numpy as np
import pandas as pd
from typing import Dict, List
import logging
import matplotlib.pyplot as plt

from utils.io import read_vcf_file
from utils.vcf_utils import combine_chrom_strands
from utils.vcf_clean import keep_common_markers_several_chr
from utils.misc import check_chromosome
from utils.plot_utils import (snp_means_plot, PCA_trained_and_projected_on_query, 
                              PCA_trained_on_query_projected_on_both, PCA_trained_and_projected_on_both)


def plot_snp_means(query_paths: List[str], reference_paths: List[str], plot_dict: Dict, 
                   output_folder: str, logger: logging.Logger) -> None:
    """
    Plots the SNP means of all the common markers between the query and the reference 
    datasets. If data for more than one .vcf file is provided, the data for all 
    chromosomes will be concatenated and displayed in the plot. First, the common 
    markers (i.e., SNPs at the same CHROM, POS, REF, and ALT) are found and the rest 
    of the SNPs are discarded. Then, the maternal and paternal strands are averaged. 
    Afterward, the mean of each SPN is computed for both datasets and, finally, 
    the SNP means are plotted.
    
    Args:
        query_paths (List[str]): paths to reference .vcf files with data for a 
        single or multiple chromosomes each.
        reference_paths (List[str]): paths to query .vcf files with data for a 
        single or multiple chromosomes each.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """

    # Define counter
    i = 1
    
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
        chrom = check_chromosome(query, reference, single=False)
        
        # Filter the reference and the query to only contain the common markers
        if type(chrom) is not list: logger.debug(f'Searching common markers in chromosome {chrom}...')
        else: logger.debug(f'Searching common markers in all chromosomes...')
        reference, query, _, _ = keep_common_markers_several_chr(reference, query, logger)
        
        # Obtain SPNs data
        query_snps = query['calldata/GT']
        reference_snps = reference['calldata/GT']
        
        # Convert SNPs data from MxNx2 to MxN shape in the query
        # and average maternal and paternal strands
        logger.debug(f'Separating maternal and paternal strands in the query '\
                    'and averaging both strands...')
        length_query, num_dogs_query, num_strands_query = query_snps.shape
        query_snps = query_snps.reshape(length_query, num_dogs_query*2).T
        query_snps = combine_chrom_strands(query_snps)
        
        # Convert SNPs data from MxNx2 to MxN shape in the reference
        # and average maternal and paternal strands
        logger.debug(f'Separating maternal and paternal strands in the reference '\
                    'and averaging both strands...')
        length_reference, num_dogs_reference, num_strands_reference = reference_snps.shape
        reference_snps = reference_snps.reshape(length_reference, num_dogs_reference*2).T
        reference_snps = combine_chrom_strands(reference_snps)
        
        # Concatenate SNPs data
        if i == 1:
            all_query_snps = query_snps
            all_reference_snps = reference_snps
        else:
            all_query_snps = np.concatenate((all_query_snps, query_snps), axis=1)
            all_reference_snps = np.concatenate((all_reference_snps, reference_snps), axis=1)
        
    logger.info(f'{all_query_snps.shape[1]} common SNPs in total between the query and the reference datasets.')
    
    # Plot the SNP means for the common markers between the query and the reference
    snp_means_plot(all_query_snps, all_reference_snps, plot_dict, output_folder, logger)

    
def plot_pca(query_paths: List[str], reference_paths: List[str], train_both: bool, 
             PCA_plot_dict: Dict, output_folder: str, logger: logging.Logger) -> None:
    """
    Applies reduction of dimensionality through Principal Component Analysis (PCA) 
    to the SNPs data, transforming the high-dimensional SNP sequences to 2D points.
    If data for more than one chromosome is provided, the SNPs are concatenated.
    The maternal and paternal strands are averaged before applying the PCA.
    
    If only data from the query is provided, the PCA is trained and projected on 
    that data. If both data from the query and the reference are provided, the PCA
    is trained and projected on both datasets, except when train_both = False, then
    the data is trained on the query and projected on both the query and the reference. 
    
    Args:
        query_paths (List[str]): paths to reference .vcf files with data for a 
        single chromosome each.
        reference_paths (List[str]): paths to query .vcf files with data for a 
        single chromosome each.
        train_both (bool): train on both datasets or only on the query data.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
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
            logger.debug(f'Cleaning chromosome {chroms_r[0]}...')
            
            # Common markers
            reference, query, _, _ = search_and_keep_common_markers_single_chr(reference, query, logger)
        
            # Obtain snps for a specific chromosome
            reference_snps = reference['calldata/GT']
        
            # Convert snps from MxNx2 to MxN shape (averaged)
            length_reference, num_dogs_reference, num_strands_reference = reference_snps.shape
            reference_snps = reference_snps.reshape(length_reference, num_dogs_reference*2).T
            reference_snps = combine_chrom_strands(reference_snps)
            
            # Concatenate snps for all chromosomes
            if i == 1:
                all_reference_snps = reference_snps
            else:
                all_reference_snps = np.concatenate((all_reference_snps, reference_snps), axis=1)
        
        # Obtain snps for a specific chromosome
        query_snps = query['calldata/GT']
        
        # Convert snps from MxNx2 to MxN shape (averaged)
        length_query, num_dogs_query, num_strands_query = query_snps.shape
        query_snps = query_snps.reshape(length_query, num_dogs_query*2).T
        query_snps = combine_chrom_strands(query_snps)
        
        # Concatenate snps for all chromosomes
        if i == 1:
            all_query_snps = query_snps
        else:
            all_query_snps = np.concatenate((all_query_snps, query_snps), axis=1)
        
        logger.info(f'{all_query_snps.shape[1]} SNPs in total.')
    
    # Plot PCA
    if reference_paths is None:
        # Trained and projected on the query
        PCA_trained_and_projected_on_query(all_snps_query, PCA_plot_dict, output_folder)
    elif not train_both:
        # Trained on the query and projected on both the query and the reference
        PCA_trained_on_query_projected_on_both(all_snps_query, all_snps_reference, PCA_plot_dict, output_folder)
    else:
        # Trained and projected on both the query and the reference
        PCA_trained_and_projected_on_both(all_snps_query, all_snps_reference, PCA_plot_dict, output_folder)