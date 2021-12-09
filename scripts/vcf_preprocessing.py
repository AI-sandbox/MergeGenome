################################################################################
# Removes undesired samples by sample ID (breed name) 
# Removes ambiguous SNPs
# Corrects SNP flips
# Removes SNP mismatches
# Renames missings (from -1 to ".")
################################################################################

import sys
import os
sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')
from utils.vcf_utils import read_vcf_file, filter_samples, rename_missings, write_vcf_file
from utils.vcf_preprocessing import search_and_remove_ambiguous_snps, search_and_correct_flips_by_pos, search_and_remove_mismatches_by_pos
from utils.track import track

################################################################################

## USER INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define path to input .vcf files with data for each chromosome to be preprocessed
# PATH1 is the path to dataset 1 (short number of SNPs), PATH2 is the path to dataset 2 (large number of SNPs)
PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/whole_genome/chr{}_unfiltered_phased.vcf'
PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/All_Pure_150k_chr{}.vcf'

## Define path to output .vcf files with preprocessed data for each chromosome
OUTPUT_PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/whole_genome/cleaned/chr{}_unfiltered_phased_cleaned.vcf'
OUTPUT_PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/cleaned/All_Pure_150k_chr{}_cleaned.vcf'

## Define name given to datasets 1 and 2 to be preprocessed
dataset1_name = 'whole_genome'
dataset2_name = 'array'

## Define name of .txt that will contain the tracking of the script (saved prints)
# Note: track_name should end in .txt
track_name = 'preprocessing_{}_{}.txt'.format(dataset1_name, dataset2_name)

## Define path to output .txt file that will store the tracking of the script
track_path = '/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/output/{}'.format(track_name)

## Define list with sample IDs or substring of sample IDs to be removed 
# Note: the strings can be in uppercase or lowercase, they will be standardized latter
# Example: substrings = ['wolf', 'fox', 'coyote', 'dhole', 'GDJK_GDJK_24316', 'GDJK_GDJK_24589', 'GoldenJackal01', 'WO001_895', 'WO002_732', 'WO003_636']
substrings = ['wolf', 'fox', 'coyote', 'dhole', 'GDJK_GDJK_24316', 'GDJK_GDJK_24589', 'GoldenJackal01', 'WO001_895', 'WO002_732', 'WO003_636']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Remove content of {TRACK_PATH} if it already exists
if os.path.exists(track_path):
    os.remove(track_path)
    
for i in range(1, 39):
    ## For each chromosome number...
    track('\n------------------------- Preprocessing chromosome {} -------------------------\n'.format(i), track_path)
    
    ## Define path to input .vcf files for chromosome i of datasets 1 and 2
    path1 = PATH1.format(i)
    path2 = PATH2.format(i)
    
    ## Define paths to output .vcf files with preprocessed data for chromosome i of datasets 1 and 2
    output_path_1 = OUTPUT_PATH1.format(i)
    output_path_2 = OUTPUT_PATH2.format(i)
    
    ## Read input .vcf files for chromosome i of datasets 1 and 2
    data1 = read_vcf_file(path1)
    data2 = read_vcf_file(path2)
    
    track('{} SNPs and {} samples in {} dataset before preprocessing'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
    track('{} SNPs and {} samples in {} dataset before preprocessing\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_path)
    
    ## Subset samples to remove undesired breeds/species with name in substrings list
    # For dataset 1
    track('Subsetting samples in {} dataset'.format(dataset1_name), track_path)
    data1 = filter_samples(data1, substrings, track_path)
    track('{} SNPs and {} samples in {} dataset after subsetting samples\n'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
    
    # For dataset 2
    track('Subsetting samples samples in {} dataset'.format(dataset2_name), track_path)
    data2 = filter_samples(data2, substrings, track_path)
    track('{} SNPs and {} samples in {} dataset after subsetting samples\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_path)
    
    ## Remove ambiguous SNPs
    # For dataset 1
    track('Searching and removing ambiguous SNPs in {} dataset'.format(dataset1_name), track_path)
    data1 = search_and_remove_ambiguous_snps(data1, track_path)
    track('{} SNPs and {} samples in {} dataset after removing ambiguous SNPs\n'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
    
    # For dataset 2
    track('Searching and removing ambiguous SNPs in {} dataset'.format(dataset2_name), track_path)
    data2 = search_and_remove_ambiguous_snps(data2, track_path)
    track('{} SNPs and {} samples in {} dataset after removing ambiguous SNPs\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_path)
    
    ## Correct SNP flips in dataset 2 with respect to the dataset 1
    track('Searching and correcting SNP flips in {} dataset'.format(dataset2_name), track_path)
    data2 = search_and_correct_flips_by_pos(data1, data2, track_path)
    
    ## Remove SNPs with mismatches between datasets 1 and 2
    track('\nSearching and removing SNPs with mismatches between {} and {} datasets'.format(dataset1_name, dataset2_name), track_path)
    data1, data2 = search_and_remove_mismatches_by_pos(data1, data2, track_path)
    track('{} SNPs and {} samples in {} dataset after removing SNPs with mismatches'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
    track('{} SNPs and {} samples in {} dataset after removing SNPs with mismatches\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_path)
    
    ## Rename missing values from -1 to "." in dataset 2
    # Note: in this example, dataset 1 does not contain any missing value
    track('Renaming missings in calldata/GT from -1 to "." in {} dataset\n'.format(dataset2_name), track_path)
    data2 = rename_missings(data2, -1, '.')
    
    track('{} SNPs and {} samples in {} dataset after preprocessing'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
    track('{} SNPs and {} samples in {} dataset after preprocessing\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_path)
    
    ## Write preprocessed data in .vcf file saved in paths output_path_1 (for dataset 1) and output_path_2 (for dataset 2)
    write_vcf_file(data1, output_path_1)
    write_vcf_file(data2, output_path_2)
    
    track('Finished writing preprocessed data'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
