################################################################################
# Removes all the SNPs with a mean absolute difference higher than 
# a given threshold
# Renames missings (from -1 to ".")
################################################################################

import sys
import os
sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')
from utils.vcf_utils import read_vcf_file, combine_chrom_strands, rename_missings, write_vcf_file
from utils.vcf_preprocessing import search_and_remove_snps_different_means
from utils.track import track

################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## USER INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define path to input .vcf files with data for each chromosome to be preprocessed
# PATH1 is the path to dataset 1 (large number of SNPs), PATH2 is the path to dataset 2 (short number of SNPs)
PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/whole_genome/cleaned/chr{}_unfiltered_phased_cleaned.vcf'
PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/cleaned/All_Pure_150k_chr{}_cleaned.vcf'

## Define path to output .vcf files with removed SNPs with a mean absolute difference higher than a threshold for each chromosome
OUTPUT_PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/whole_genome/cleaned_rm_01/chr{}_unfiltered_phased_cleaned_rm_01.vcf'
OUTPUT_PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/cleaned_rm_01/All_Pure_150k_chr{}_cleaned_rm_01.vcf'

## Define name given to datasets 1 and 2 to be preprocessed
dataset1_name = 'whole_genome_cleaned'
dataset2_name = 'array_cleaned'

## Define threshold
# All common SNPs between datasets 1 and 2 with a mean absolute difference higher this threshold will be removed
threshold = 0.1

## Define name of .txt that will contain the tracking of the script (saved prints)
# Note: track_name should end in .txt
track_name = 'removing_snps_with_different_means_{}_{}.txt'.format(dataset1_name, dataset2_name)

## Define path to output .txt file that will store the tracking of the script
track_path = '/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/output/{}'.format(track_name)

## Remove content of {TRACK_PATH} if it already exists
if os.path.exists(track_path):
    os.remove(track_path)
    
for i in range(27, 39):
    ## For each chromosome number...
    track('\n------------------------- Reading chromosome {} -------------------------\n'.format(i), track_path)
    
    ## Define path to input .vcf files for chromosome i of datasets 1 and 2
    path1 = PATH1.format(i)
    path2 = PATH2.format(i)
    
    ## Define paths to output .vcf files with removed SNPs for chromosome i of datasets 1 and 2
    output_path_1 = OUTPUT_PATH1.format(i)
    output_path_2 = OUTPUT_PATH2.format(i)
    
    ## Read input .vcf files for chromosome i of datasets 1 and 2
    data1 = read_vcf_file(path1)
    data2 = read_vcf_file(path2)
        
    track('{} SNPs and {} samples in {} dataset'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
    track('{} SNPs and {} samples in {} dataset\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_path)
        
    ## Remove all SNPs with a mean absolute difference higher than the threshold
    track('Removing all SNPs with a mean absolute difference higher than {}'.format(dataset1_name), track_path)
    data1, data2 = search_and_remove_snps_different_means(data1, data2, threshold, track_path)
    
    track('\n{} SNPs and {} samples in {} dataset after removing SNPs with a mean absolute difference higher than {}'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name, threshold), track_path)
    track('{} SNPs and {} samples in {} dataset after removing SNPs with a mean absolute difference higher than {}\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name, threshold), track_path)
    
    ## Rename missing values from -1 to "." in dataset 2
    # Note: in this example, dataset 1 does not contain any missing value
    track('Renaming missings in calldata/GT from -1 to "." in {} dataset'.format(dataset2_name), track_path)
    data2 = rename_missings(data2, -1, '.')
    track('All missings renamed\n'.format(dataset2_name), track_path)
    
    ## Write preprocessed data in .vcf file saved in paths output_path_1 (for dataset 1) and output_path_2 (for dataset 2)
    write_vcf_file(data1, output_path_1)
    write_vcf_file(data2, output_path_2)
    
    track('Finished writing data with removed SNPs'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
    
