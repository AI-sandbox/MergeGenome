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
from itertools import zip_longest
import matplotlib.pyplot as plt

from utils.io import read_vcf_file
from utils.vcf_utils import combine_chrom_strands
from utils.vcf_clean import keep_common_markers_several_chr
from utils.misc import check_chromosome
from utils.plot_utils import (
    snp_means_plot, 
    PCA_trained_and_projected_on_query, 
    PCA_trained_on_query_projected_on_both, 
    PCA_trained_and_projected_on_both
)

def plot_snp_means(query_paths: List[str], reference_paths: List[str], indexes_path: str, 
                   plot_dict: Dict, output_folder: str, logger: logging.Logger) -> None:
    """
    Plots the SNP means of all the common markers between the query and the reference 
    datasets. If data for more than one .vcf file is provided, the data for all 
    chromosomes will be concatenated and displayed in the plot. First, the common 
    markers (i.e., SNPs at the same CHROM, POS, REF, and ALT) are found and the rest 
    of the SNPs are discarded. Then, the maternal and paternal strands are averaged. 
    Afterward, the mean of each SPN is computed for both datasets and, finally, 
    the SNP means are plotted. If indexes_path is not None, then the SNPs at the
    specified indexes will be plotted in a different color.
    
    Args:
        query_paths (List[str]): paths to reference .vcf files with data for a 
        single or multiple chromosomes each.
        reference_paths (List[str]): paths to query .vcf files with data for a 
        single or multiple chromosomes each.
        indexes_path (str): path to indexes to be plotted in a different color.
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
        
        # Ensure the reference and the query contain data for the same chromosome/s
        # and obtain the chromosome/s in particular
        chrom = check_chromosome(query, reference, single=False)
        
        # Filter the reference and the query to only contain the common markers
        if not isinstance(chrom, list): logger.debug(f'Searching common markers in chromosome {chrom}...')
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
        
        i += 1
        
    logger.info(f'{all_query_snps.shape[1]} common SNPs in total between the query and the '\
                'reference datasets.')
    
    # Read indexes of SNPs that will be plotted in a different color
    indexes = np.load(indexes_path)
    
    # Plot the SNP means for the common markers between the query and the reference
    snp_means_plot(all_query_snps, all_reference_snps, indexes, plot_dict, output_folder, logger)


def plot_pca(query_paths: List[str], reference_paths: List[str], train_query: bool, 
             PCA_plot_dict: Dict, output_folder: str, logger: logging.Logger) -> None:
    """
    Plots the first two components of the Principal Component Analysis (PCA) of the 
    data, transforming the high-dimensional SNP sequences to 2D points. If the query 
    is provided but not the reference, the PCA is applied to the query. If both the 
    query and the reference are provided, the PCA is trained and projected on both datasets. 
    If train_both = False, then the data is trained on the query and projected on both.
    If data for more than one .vcf file is provided, the data for all chromosomes will 
    be concatenated and displayed in the plot. First, the maternal and paternal strands 
    are averaged. Afterward, the PCA is computed and, finally, the first two components 
    are plotted.
    
    Args:
        query_paths (List[str]): paths to reference .vcf files with data for a 
        single or multiple chromosomes each.
        reference_paths (List[str]): paths to query .vcf files with data for a 
        single or multiple chromosomes each.
        train_query (bool): to train the PCA the PCA only on the query instead of 
        on both datasets.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """
    
    # Define counter
    i = 1

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
        
        # Ensure the reference and the query contain data for the same chromosome/s
        # and obtain the chromosome/s in particular
        chrom = check_chromosome(query, reference, single=False)
        
        # Filter the reference and the query to only contain the common markers
        if not isinstance(chrom, list): logger.debug(f'Searching common markers in chromosome {chrom}...')
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
        
        if reference_path is not None:
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
            if reference_path is not None:
                all_reference_snps = reference_snps
        else:
            all_query_snps = np.concatenate((all_query_snps, query_snps), axis=1)
            if reference_path is not None:
                all_reference_snps = np.concatenate((all_reference_snps, reference_snps), axis=1)
        
        i += 1
        
    logger.info(f'{all_query_snps.shape[1]} common SNPs in total between the query and the '\
                'reference datasets.')
        
    # Plot PCA
    if reference_paths is None:
        # Train and project PCA on the query
        PCA_trained_and_projected_on_query(all_query_snps, PCA_plot_dict, output_folder, logger)
    elif train_query:
        # Train PCA on the query and project on both the query and the reference
        PCA_trained_on_query_projected_on_both(
            all_query_snps, all_reference_snps, PCA_plot_dict, output_folder, logger
        )
    else:
        # Train and project PCA on both the query and the reference
        PCA_trained_and_projected_on_both(
            all_query_snps, all_reference_snps, PCA_plot_dict, output_folder, logger
        )