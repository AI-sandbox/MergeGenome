################################################################################
# Searches for the indexes of the SNPs of a dataset
# that are also present in another dataset
# The indexes of the common markers are saved in a numpy array
################################################################################

import sys
import os
import numpy as np
import pandas as pd

sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')
from utils.vcf_utils import read_vcf_file, filter_by_chromosome
from utils.vcf_preprocessing import search_and_keep_common_markers_several_chr
from utils.track import track

################################################################################

## USER INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define base path to .vcf files with all chromosomes data
#PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/merged/merged_filteredsamples_allchr.vcf'
#PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/cleaned/All_Pure_150k_allchr_cleaned.vcf'

PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/embark/embark_allchr.vcf'
PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/whole_to_array/X_whole_to_array_allchr.vcf'

## Define path were the numpy array with the indexes of the common SNPs will be saved
PATH_INDEXES = '/home/users/miriambt/my_work/dog-gen-to-phen/imputation/data/indexes_{}_in_{}'

## Define name of datasets 1 and 2
#dataset1_name = 'merged'
#dataset2_name = 'array'

dataset1_name = 'embark'
dataset2_name = 'whole_to_array'

## Define name of .txt to contain the tracking of the script
track_name = 'searching_snps_of_{}_in_{}.txt'.format(dataset1_name, dataset2_name)

## Define path to output .txt file that will store the tracking of the script
track_path = '/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/output/{}'.format(track_name)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Remove content of .txt with name track_name if it already exists
if os.path.exists(track_path):
    os.remove(track_path)
    
## Read input .vcf files for chromosome i of datasets 1 and 2
data1 = read_vcf_file(PATH1)
data2 = read_vcf_file(PATH2)

track('{} SNPs and {} samples in {} dataset'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_path)
track('{} SNPs and {} samples in {} dataset\n'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_path)

## Obtain indexes of the SNPs in data1 that are also in data2
_, _, indexes1, indexes2 = search_and_keep_common_markers_several_chr(data1, data2, track_path)

## Save the found indexes
# np.save(PATH_INDEXES.format(dataset1_name, dataset2_name), indexes1)

## Ensure the indexes were correctly found
count = 0
for i in range (len(indexes1)):
    if(data1['variants/REF'][indexes1[i]] != data2['variants/REF'][indexes2[i]] or 
       data1['variants/ALT'][indexes1[i]][0] != data2['variants/ALT'][indexes2[i]][0]):
        count += 1

assert count == 0, 'All the indexes found do not correspond to common SNPs. Check the script.'
