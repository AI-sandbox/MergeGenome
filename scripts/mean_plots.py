################################################################################
# Makes SNP plot and save it in output path
################################################################################

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')
from utils.vcf_utils import read_vcf_file, combine_chrom_strands
from utils.vcf_preprocessing import search_and_keep_common_markers_single_chr
from utils.plot_utils import snp_means
from utils.track import track

################################################################################

## USER INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define base name to paths
PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/embark/embark_chr{}.vcf'
PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array_imputed_subset/imputed_array_subset_chr{}.vcf'

## Define path to directory that will contain the folder with the SNPs means plot
output_path = '/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/output/SNP_means_plots/'

## Define name of datasets 1 and 2
dataset1_name = 'Embark'
dataset2_name = 'Imputed Array Subsetted'

## Define name of .txt that will contain the tracking of the script (saved prints)
# Note: track_name should end in .txt
track_name = 'plotting_snp_means_{}_and_{}.txt'.format(dataset1_name, dataset2_name)

## Define path to output .txt file that will store the tracking of the script
track_path = '/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/output/{}'.format(track_name)

## Define plot characteristics
plot_params = {
        'FONTSIZE' : 25,
        'FIG_WIDTH' : 26,
        'FIG_HEIGHT' : 15,
        's' : 0.1,
        'alpha' : 1.0,
        'color' : '#306998'
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## If it does not already exist, create directory that will contain the folder with PCA plots
if not os.path.exists(output_path):
    os.makedirs(output_path)

## Remove content of .txt with name track_name if it already exists
if os.path.exists(track_path):
    os.remove(track_path)
    
for i in range(1, 39):
    ## For each chromosome number...
    track('\n------------------------- Reading data of chromosome {} -------------------------\n'.format(i), track_path)
    
    ## Define path to input .vcf files for chromosome i of datasets 1 and 2
    path1 = PATH1.format(i)
    path2 = PATH2.format(i)
    
    ## Read input .vcf files for chromosome i of datasets 1 and 2
    data1 = read_vcf_file(path1)
    data2 = read_vcf_file(path2)
    
    ## Common markers
    data1_chr, data2_chr, _, _ = search_and_keep_common_markers_single_chr(data1, data2, track_path)
    
    ## Obtain snps for chromosome i
    snps_1 = data1_chr['calldata/GT']
    snps_2 = data2_chr['calldata/GT']
        
    ## Convert snps from MxNx2 to MxN shape (averaged)
    length_1, num_dogs_1, num_strands_1 = snps_1.shape
    snps_1 = snps_1.reshape(length_1, num_dogs_1*2).T
    snps_1 = combine_chrom_strands(snps_1)

    length_2, num_dogs_2, num_strands_2 = snps_2.shape
    snps_2 = snps_2.reshape(length_2, num_dogs_2*2).T
    snps_2 = combine_chrom_strands(snps_2)
       
    ## Concatenate snps for all chromosomes
    if i == 1:
        all_snps_1 = snps_1
        all_snps_2 = snps_2
    else:
        all_snps_1 = np.concatenate((all_snps_1, snps_1), axis=1)
        all_snps_2 = np.concatenate((all_snps_2, snps_2), axis=1)
        
track('\n {} common SNPs in total between {} and {} datasets'.format(all_snps_1.shape[1], dataset1_name, dataset2_name), track_path)

## Plot the SNP means for the SNPs that are common markers between all_snps_1 and all_snps_2
snp_means(all_snps_1, all_snps_2, plot_params, dataset1_name, dataset2_name, output_path)
