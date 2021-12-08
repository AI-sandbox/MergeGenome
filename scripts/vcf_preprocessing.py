################################################################################
# Removes undesired samples by sample ID (breed name) 
# Removes ambiguous SNPs
# Corrects SNP flips
# Removes SNP mismatches
# Renames missings (from "-1" to ".")
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

## Define path to input .vcf files for each chromosome to be preprocessed
PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/whole_genome/chr{}_unfiltered_phased.vcf'
PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/All_Pure_150k_chr{}.vcf'

## Define path to output .vcf files with preprocessed data for each chromosome
OUTPUT_PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/whole_genome/cleaned/chr{}_unfiltered_phased_cleaned.vcf'
OUTPUT_PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/cleaned/All_Pure_150k_chr{}_cleaned.vcf'

## Define name of datasets 1 and 2 to be preprocessed
dataset1_name = 'whole_genome'
dataset2_name = 'array'

## Define name of .txt that will contain the tracking of the script
track_name = 'preprocessing_{}_{}.txt'.format(dataset1_name, dataset2_name)

## Define samples or substrings of samples to be removed (in uppercase or lowercase)
substrings = ['wolf', 'fox', 'coyote', 'dhole', 'GDJK_GDJK_24316', 'GDJK_GDJK_24589', 'GoldenJackal01', 'WO001_895', 'WO002_732', 'WO003_636']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Remove content of .txt with name track_name if it already exists
if os.path.exists('../output/{}'.format(track_name)):
    os.remove('../output/{}'.format(track_name))
    
## For every chromosome
for i in range(1, 39):
    track('\n------------------------- Preprocessing chromosome {} -------------------------\n'.format(i), track_name)
    
    ## Define paths to input .vcf files containing data for chromosome i to be preprocessed
    path1 = PATH1.format(i)
    path2 = PATH2.format(i)
    
    ## Define paths to output .vcf files with preprocessed data for chromosome i
    output_path_1 = OUTPUT_PATH1.format(i)
    output_path_2 = OUTPUT_PATH2.format(i)
    
    ## Read input .vcf files for chromosome i
    data1 = read_vcf_file(path1)
    data2 = read_vcf_file(path2)
    
    track('{} SNPs and {} samples in {} dataset before preprocessing'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_name)
    track('{} SNPs and {} samples in {} dataset before preprocessing\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_name)
    
    ## Subset samples to remove strange breeds or species
    track('Subsetting samples in {} dataset'.format(dataset1_name), track_name)
    data1 = filter_samples(data1, substrings, track_name)
    track('{} SNPs and {} samples in {} dataset after subsetting samples\n'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_name)
    
    track('Subsetting samples samples in {} dataset'.format(dataset2_name), track_name)
    data2 = filter_samples(data2, substrings, track_name)
    track('{} SNPs and {} samples in {} dataset after subsetting samples\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_name)
    
    ## Remove ambiguous SNPs
    track('Searching and removing ambiguous SNPs in {} dataset'.format(dataset1_name), track_name)
    data1 = search_and_remove_ambiguous_snps(data1, track_name)
    track('{} SNPs and {} samples in {} dataset after removing ambiguous SNPs\n'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_name)
    
    track('Searching and removing ambiguous SNPs in {} dataset'.format(dataset2_name), track_name)
    data2 = search_and_remove_ambiguous_snps(data2, track_name)
    track('{} SNPs and {} samples in {} dataset after removing ambiguous SNPs\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_name)
    
    ## Correct SNP flips in the second dataset with respect to the first dataset
    track('Searching and correcting SNP flips in {} dataset'.format(dataset2_name), track_name)
    data2 = search_and_correct_flips_by_pos(data1, data2, track_name)
    
    ## Remove SNPs with mismatches between the two datasets
    track('\nSearching and removing SNPs with mismatches between {} and {} datasets'.format(dataset1_name, dataset2_name), track_name)
    data1, data2 = search_and_remove_mismatches_by_pos(data1, data2, track_name)
    track('{} SNPs and {} samples in {} dataset after removing SNPs with mismatches'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_name)
    track('{} SNPs and {} samples in {} dataset after removing SNPs with mismatches\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_name)
    
    track('Renaming missings in calldata/GT from -1 to "." in {} dataset\n'.format(dataset2_name), track_name)
    data2 = rename_missings(data2, -1, '.')
    
    track('{} SNPs and {} samples in {} dataset after preprocessing'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_name)
    track('{} SNPs and {} samples in {} dataset after preprocessing\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_name)
    
    ## Write output in .vcf file
    write_vcf_file(data1, output_path_1)
    write_vcf_file(data2, output_path_2)
    
    track('Finished writing preprocessed data'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_name)
