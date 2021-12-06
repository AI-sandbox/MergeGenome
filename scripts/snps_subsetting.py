import sys
import os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')

from utils.vcf_utils import read_vcf_file, write_vcf_file
from utils.vcf_preprocessing import search_and_keep_common_markers
from utils.track import track

#############################################
### INPUT
#############################################

## Define base path to .vcf files containing the SNPs to be subsetted, for each chromosome
PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/whole_genome/cleaned/chr*_unfiltered_phased_cleaned.vcf'
PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/imputed/All_Pure_150k_chr*_cleaned_imputed.vcf'

## Define base path to .vcf files containing the SNPs that we want to keep (if present), for each chromosome
PATH3 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/embark_old/subset_chr*.vcf'

## Define base path to .vcf files to save the subsetted data for each chromosome
OUTPUT_PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/embark/embark_chr*.vcf'
OUTPUT_PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array_imputed_subset/imputed_array_subset_chr*.vcf'

## Define name of datasets 1 and 2
dataset1_name = 'whole_genome'
dataset2_name = 'imputed_array'
dataset3_name = 'embark'
dataset4_name = 'imputed_array_subsetted'

## Define name of .txt to contain the tracking of the script
track_name = 'trash_subsetting_{}_{}_to_{}_{}.txt'.format(dataset1_name, dataset2_name, dataset3_name, dataset4_name)

## Remove content of .txt with name track_name if it already exists
if os.path.exists('../output/{}'.format(track_name)):
    os.remove('../output/{}'.format(track_name))

## For every chromosome
for i in range(20, 39):
    track('\n------------------------- Chromosome{} -------------------------\n'.format(i), track_name)
    
    ## Define paths for chromosome i
    path1 = PATH1.replace('*', str(i))
    path2 = PATH2.replace('*', str(i))
    path3 = PATH3.replace('*', str(i))
    
    output_path_1 = OUTPUT_PATH1.replace('*', str(i))
    output_path_2 = OUTPUT_PATH2.replace('*', str(i))
    
    ## Read .vcf files for chromosome i using skikit-allel
    data1 = read_vcf_file(path1)
    data2 = read_vcf_file(path2)
    data3 = read_vcf_file(path3)
    
    track('{} SNPs and {} samples in {} dataset'.format(len(data1['variants/ID']), len(data1['samples']), dataset1_name), track_name)
    track('{} SNPs and {} samples in {} dataset'.format(len(data2['variants/ID']), len(data2['samples']), dataset2_name), track_name)
    track('{} SNPs and {} samples in {} dataset'.format(len(data3['variants/ID']), len(data3['samples']), dataset3_name), track_name)
    
    ## Keep common SNPs of the first and the second datasets, if present in the third dataset
    # The SNPs that are in the first or second dataset but not in the third dataset are removed
    data1, data3 = search_and_keep_common_markers(data1, data3, track_name)
    data2, data3 = search_and_keep_common_markers(data2, data3, track_name)
    
    ## Write subsetted data in output .vcf file
    write_vcf_file(data1, output_path_1)
    write_vcf_file(data2, output_path_2)
    
    track('{} SNPs and {} samples in {} dataset'.format(len(data1['variants/ID']), len(data1['samples']), dataset3_name), track_name)
    track('{} SNPs and {} samples in {} dataset'.format(len(data2['variants/ID']), len(data2['samples']), dataset4_name), track_name)