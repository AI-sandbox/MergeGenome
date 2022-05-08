################################################################################
# Functions to assist with the manipulation of .vcf files:
# * Obtain chromosomes with available data
# * Filtering data by chromosome
# * Rename chromosome nomenclature
# * Obtaining percentage of SNPs with some missing value (encoded as -1)
################################################################################

import os
import numpy as np
import pandas as pd
import random
from typing import List
import logging

def obtain_chromosomes(vcf_data: dict) -> List:
    """
    Obtains the name of the chromosomes with available data.
    
    Args:
        vcf_data (dict): allel.read_vcf output.
    
    Returns:
        (List): list with the unique chromosomes in 'variants/CHROM'
        in order of appearance.

    """
    
    # Obtain the chromosomes in data as unique keys where order is preserved
    chroms = dict.fromkeys(vcf_data['variants/CHROM'])
    
    # Convert to list
    return list(chroms)


def obtain_renamed_chrom(rename_chr: bool, actual_chrom: str, rename_map: dict) -> str:
    """
    Obtains the new notation for a specific chromosome. If rename_chr=True and 
    rename_map=None, the variants/CHROM format will change from '<chrom_number>' 
    to 'chr<chrom_number>' or vice-versa. If the notation of a chromosome is in 
    neither of these formats, its notation will remain the same. If rename_chr=True
    and rename_map is specified, the renaming will be done with a mapping for all 
    the chromosomes to modify.
    
    Args:
        rename_chr (bool): rename (or not) the chromosome nomenclature.
        actual_chrom (str): chromosome notation before renaming.
        rename_map (str): dictionary with mapping from old to new chromosome notation.
                          The keys are the old notations and the values are the new notations.
    
    Returns:
        new_chrom (str): new chromosome notation after renaming.

    """
    
    # Initially, set the new chromosome notation as the old one
    new_chrom = actual_chrom
    
    if rename_chr:
        
        if (rename_map is not None):
            
            # Extract dict from str
            rename_map = eval(rename_map)
            
            if actual_chrom in rename_map.keys():
                # Change through mapping
                new_chrom = rename_map[actual_chrom]
        
        else:
            if actual_chrom.startswith('chr'):
                if actual_chrom[3:].isdigit():
                    # Change from 'chr<chrom_number>' to '<chrom_number>'
                    new_chrom = actual_chrom[3:]
                    
            elif actual_chrom.isdigit():
                # Change from ''<chrom_number>' to 'chr<chrom_number>'
                new_chrom = f'chr{actual_chrom}'

    return new_chrom


def filter_by_chromosome(vcf_data: dict, chrom: str) -> dict:
    """
    Filters data to keep the info of a specified chromosome.
    
    Args:
        vcf_data (dict): allel.read_vcf output.
        chrom (str): chromosome in 'variants/CHROM'.
    
    Returns:
        vcf_data (dict): vcf_data with filtered data for the desired chromosome.

    """
    
    # Find the indexes of the SNPs of the particular chromosome
    indexes_chr = vcf_data['variants/CHROM'] == str(chrom)
    
    # Select data for the SNPs of the particular chromosome
    for key in vcf_data.keys():
        if key != 'samples':
            vcf_data[key] = vcf_data[key][indexes_chr]
                
    return vcf_data


def rename_chrom_field(vcf_data: dict, actual_chrom: str, new_chrom: str) -> dict:
    """
    Renames variants/CHROM in vcf_data.
    
    Args:
        vcf_data (dict): allel.read_vcf output.
        actual_chrom (str): chromosome notation before renaming.
        new_chrom (str): new chromosome notation after renaming.

    Returns:
        vcf_data (dict): vcf_data with renamed variants/CHROM notation.
    
    """

    # Rename variants/CHROM notation
    vcf_data['variants/CHROM'] = np.where(vcf_data['variants/CHROM'] == actual_chrom, new_chrom, vcf_data['variants/CHROM'])
    
    return vcf_data


def filter_samples(vcf_data: dict, substrings: List[str], logger: logging.Logger):
    """
    Subsets samples to remove undesired breeds/species with name in substrings list.
    
    Args:
        vcf_data (dict): allel.read_vcf output.
        substrings (List[str]): substrings of the samples to be removed. Example: ['wolf', 'fox', 'coyote', 'dhole', 'GDJK_GDJK_24316']. 
                                In this case, the sample with ID "IrishWolfhound01" would be removed because "wolf" is a substring of 
                                the sample id. Note that the substrings can be in upper or lowercase since they are standardized prior to
                                filtering.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        vcf_data (dict): filtered vcf_data to only keep the desired dog breeds/species.
    
    """
    
    # Standardize substrings and sample names to lowercase
    substrings = list(map(lambda x: x.lower(), substrings))
    sample_names = list(map(lambda x: x.lower(), vcf_data['samples']))
    remove_samples = []
    
    # Define list of booleans with the samples to be kept (True) and to be removed (False).
    binary_samples = []
    
    # For each sample ID, see if it contains any substring and store results in binary_samples to later filter those samples
    for sample_name in sample_names:
        
        # True if any substring is contained in sample name, false otherwise
        substring_found = any(substring in sample_name for substring in substrings)
        
        # Store result
        # Note: samples with true values will be kept. Otherwise, they will be removed
        binary_samples.append(not substring_found)
        
        # Write sample names that will be removed in track file
        if substring_found:
            remove_samples.append(sample_name)
    
    # Obtain total number of removed samples
    removed_samples = len(binary_samples) - sum(binary_samples)
    
    # Write total number of removed samples
    logger.info(f'--> {removed_samples} samples removed in total')
    
    # Filter to remove the samples corresponding to undesired breeds/species
    vcf_data['samples'] = vcf_data['samples'][binary_samples]
    vcf_data['calldata/GT'] = vcf_data['calldata/GT'][:,binary_samples,:]
    
    return vcf_data


def search_percentage_SNPs_with_missings(vcf_data: dict, miss_code: [float, str]) -> float:
    """
    Computes the percentage of SNPs with some missing value in a certain codification.
    
    Args:
        vcf_data: allel.read_vcf output.
    
    Returns:
        perc_missings: % of SNPs with some missing value in calldata/GT encoded as "-1".
    
    """
    
    ## Obtain npy matrix with all the SNPs
    npy = vcf_data['calldata/GT']
    length, num_dogs, num_strands = npy.shape
    npy = npy.reshape(length, num_dogs*2).T
    
    ## Obtain number of missing values encoded as "-1" in each of the SNPs
    num_missings = np.sum(npy == -1, axis=0)
    
    ## Obtain the % of SNPs with some missing value
    perc_missings = sum(num_missings>0)/len(num_missings)*100
    
    return perc_missings


def missings_in_each_SNP(vcf_data):
    '''
    Objective:
        - Print the % of SNPs that contain some missing value encoded as "-1". If your missings are encoded differently, 
          change the -1 in the function for the alternative missings nomenclature.
    Input:
        - vcf_data: allel.read_vcf output.
    '''
    
    ## Obtain npy matrix with SNPs
    npy = vcf_data['calldata/GT']
    length, num_dogs, num_strands = npy.shape
    npy = npy.reshape(length, num_dogs*2).T
    
    ## Obtain number of samples with some missing value encoded as "-1" for each SNP
    num_neg = np.sum(npy == -1, axis=0)

    ## Print the % of SNPs that contain some missing value
    print(sum(num_neg>0)/len(num_neg)*100)


def rename_missings(vcf_data, before, after):
    '''
    Objective: rename missings in calldata/GT in vcf_data, from {before} to {after}.
    Input:
        - vcf_data: allel.read_vcf output.
        - before: missing value before renaming. Example: "-1".
        - after: missing value after renaming. Example: ".".
    Output:
        - vcf_data: vcf_data with renamed missings in calldata/GT.
    '''
    
    ## Rename variants/CHROM from {before} to {after}
    vcf_data['calldata/GT'] = np.where(vcf_data['calldata/GT'] == before, after, vcf_data['calldata/GT']) 

    return vcf_data


def combine_chrom_strands(chrom_data):
    '''
    Input:
        - chrom_data: matrix of chromosone data (SNPs) in the shape (num_dogs*2, length). The maternal and paternal strands should be split and
          consecutively placed in the matrix.
    Output:
        - chrom_data_combined: numpy matrix of chromosone data (SNPs) with combined maternal and paternal strands. 
          Shape should be (num_dogs, length).
    '''

    ## Determine the dimensions of the chromosone data (SNPs)
    # length corresponds to the length of the chromosone strands (the number of SNPs).
    # num_samples corresponds to the number of dogs available for analysis times 2 because
    # there are two rows in the matrix for each dog, one for maternal and one for paternal data
    num_samples, length = chrom_data.shape

    ## Combine the information from the maternal and paternal strands by averaging
    chrom_data_combined = (chrom_data[::2] + chrom_data[1::2])/2

    ## Determine dimensions of reshaped data
    num_dogs, length_new = chrom_data_combined.shape

    ## Check that the data was reformatted properly
    assert (length == length_new), "Number of snps per sample changed. Check reshaping function."
    assert (num_samples == num_dogs*2), "Number of samples changed. Check reshaping function."

    return chrom_data_combined