################################################################################
# Reads .vcf file
# Writes .vcf file
# And performs small transformations to .vcf data: 
#   - Filtering by chromosome
#   - Renaming chromosome nomenclature
#   - Obtaining percentage of SNPs with some missing (encoded as -1)
#   - etc.
################################################################################

import os
import numpy as np
import pandas as pd
import random

from utils.track import track


def filter_by_chromosome(vcf_data, chrom):
    '''
    Objective:
        - Filter vcf_data to keep the SNPs of the specified chromosome.
    Input:
        - vcf_data: allel.read_vcf output.
        - chrom: CHROM to select the SNPs. Usually chr{chromosome_number} or {chromosome_number}.
    Output:
        - vcf_data: vcf_data containing only the selected SNPs of the desired chromosome.
    '''
    
    ## Find the indexes of the SNPs that correspond to the chromosome in particular
    indexes_chr = vcf_data['variants/CHROM'] == str(chrom)
    
    ## Select data for the SNPs of the chromosome in particular
    for key in vcf_data.keys():
        if key != 'samples':
            vcf_data[key] = vcf_data[key][indexes_chr]
                
    return vcf_data


def rename_chromosome(vcf_data, before, after):
    '''
    Objective: rename variants/CHROM in vcf_data. Example: from "1" to "chr1".
    Input:
        - vcf_data: allel.read_vcf output.
        - before: chromosome name before renaming. Example: "1".
        - after: chromosome name after renaming. Example: "chr1".
    Output:
        - vcf_data: vcf_data with renamed variants/CHROM.
    '''
    
    ## Rename variants/CHROM
    vcf_data['variants/CHROM'] = np.where(vcf_data['variants/CHROM'] == before, after, vcf_data['variants/CHROM']) 
    
    return vcf_data


def search_percentage_SNPs_with_missings(vcf_data):
    '''
    Objective:
        - Return the % of SNPs with some missing value. Missing values are assumed to be encoded as "-1". If your missings are encoded differently, 
          change the -1 in the function for the alternative missings nomenclature.
    Input:
        - vcf_data: allel.read_vcf output.
    Output:
        - perc_missings: % of SNPs with some missing value in calldata/GT encoded as "-1".
    '''
    
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


def filter_samples(vcf_data, substrings, track_name, verbose=False):
    '''
    Objective: subsets samples to remove undesired breeds/species with name in substrings list.
    Input:
        - vcf_data: allel.read_vcf output.
        - substrings: substrings of the samples to be removed. 
          Example: ['wolf', 'fox', 'coyote', 'dhole', 'GDJK_GDJK_24316']. In this case, the sample with ID "IrishWolfhound01" would be removed because
          "wolf" is a substring of the sample id. Note that the substrings can be in upper or lowercase since they are standardized prior to filtering.
        - verbose: if True, the names of the sample ID that are removed are written in track file.
    Output:
        - vcf_data: filtered vcf_data to only keep the desired dog breeds/species.
    '''
    
    ## Standardize substrings and sample names to lowercase
    substrings = list(map(lambda x: x.lower(), substrings))
    sample_names = list(map(lambda x: x.lower(), vcf_data['samples']))
    
    ## Define list of booleans with the samples to be kept (True) and to be removed (False).
    samples_keep_remove = []
        
    if verbose:
        track('Removed samples:', track_name)
    
    ## For each sample ID, see if it contains any substring and store results in samples_keep_remove to later filter those samples
    for sample_name in sample_names:
        ## True if any substring is contained in sample name, false otherwise
        substring_found = any(substring in sample_name for substring in substrings)
        
        ## Store result
        # Samples with true values will be kept. Otherwise, they will be removed
        samples_keep_remove.append(not substring_found)
        
        ## Write sample names that will be removed in track file
        if verbose and substring_found:
            track('-> ' + sample_name, track_name)
    
    ## Obtain total number of removed samples
    removed_samples = len(samples_keep_remove) - sum(samples_keep_remove)
    
    if verbose and removed_samples == 0:
        track('None', track_name)
    
    ## Write total number of removed samples
    track('--> {} samples removed in total'.format(removed_samples), track_name)
    
    ## Filter to remove the samples corresponding to undesired breeds/species
    vcf_data['samples'] = vcf_data['samples'][samples_keep_remove]
    vcf_data['calldata/GT'] = vcf_data['calldata/GT'][:,samples_keep_remove,:]
    
    return vcf_data


def train_test_split(vcf_data, train_percentage):
    '''
    Objective: perform train/test splitting. A set of samples are randomly selected to go to the train set (as much as train_percentage).
    The rest of samples go to the train set.
    Input:
        - vcf_data: allel.read_vcf output.
        - train_percentage: percentage of samples that will go to the train set. 
    Output:
        - vcf_data_train: vcf_data with the samples of the train set.
        - vcf_data_test: vcf_data with the samples of the test set.
    '''
    ## Obtain number of samples and SNPs in vcf_data
    samples = len(vcf_data['samples'])
    snps = len(vcf_data['variants/ID'])
    
    ## Define vector with all the sample indexes (from 0 to samples-1)
    row_idxs = list(range(0,samples))

    ## Define number of samples that will go to the train and test sets
    n_train = round(train_percentage*samples)
    n_test = samples - n_train

    ## Randomly sample n_train numbers from row_idxs without replacement
    # Set seed for reproducibility
    random.seed(0)
    train_idxs = random.sample(row_idxs, n_train)

    ## Ensure all numbers in train_idxs are unique
    assert len(np.unique(train_idxs)) == n_train, 'There are repeated row indexes for the train set. Check again.'

    ## Obtain the indexes of the samples that will go to the test set
    test_idxs = list(set(row_idxs) - set(train_idxs))

    ## Ensure all numbers in test_idxs are unique
    assert len(np.unique(test_idxs)) == n_test, 'There are repeated row indexes for the test set. Check again.'
    
    ## Make two copies of vcf_data
    vcf_data_train = vcf_data.copy()
    vcf_data_test = vcf_data.copy()
    
    ## Modify the first copy to filter by the samples that go to the train set
    vcf_data_train['samples'] = vcf_data_train['samples'][train_idxs]
    vcf_data_train['calldata/GT'] = vcf_data_train['calldata/GT'][:,train_idxs,:]
    
    ## Modify the second copy to filter by the samples that go to the test set
    vcf_data_test['samples'] = vcf_data_test['samples'][test_idxs]
    vcf_data_test['calldata/GT'] = vcf_data_test['calldata/GT'][:,test_idxs,:]
    
    assert len(vcf_data_train['samples']) + len(vcf_data_test['samples']) == len(vcf_data['samples']), 'The number of samples in the train and test sets does not match the original number of samples'
    
    return vcf_data_train, vcf_data_test
    

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
