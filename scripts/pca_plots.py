################################################################################
# Makes PCA plots...
# Trained and projected on the first dataset
# Trained and projected on the second dataset
# Trained on the first dataset and projected on both the first and the second dataset
# Trained on the second dataset and projected on both the first and the second dataset
# Trained and projected on both the first and the second dataset
################################################################################

import sys
import os
import numpy as np

sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')
from utils.vcf_utils import read_vcf_file, combine_chrom_strands
from utils.plot_utils import PCA_2D_trained_and_projected_on_dataset, PCA_2D_trained_on_dataset1_projected_on_both, PCA_2D_trained_and_projected_on_both
from utils.track import track

################################################################################

## USER INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define base path to .vcf files for preprocessed data for each chromosome
PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/embark/embark_chr*.vcf'
PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array_imputed_subset/imputed_array_subset_chr*.vcf'

## Define path to directory that will contain the folder with the PCA plots
output_path = '/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/output/PCA_plots/'

## Define plot characteristics
plot_params = {
        'FONTSIZE' : 16,
        'FIG_WIDTH' : 12,
        'FIG_HEIGHT' : 10,
        's' : 15,
        'alpha' : 0.4
}

## Define name of datasets 1 and 2
dataset1_name = 'Embark'
dataset2_name = 'Imputed Array Subsetted'

## Define color of data points for datasets 1 and 2
color1 = '#259988' 
color2 = '#EBD0A1'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define path that will contain the PCA plots
output_path2 = output_path + dataset1_name.replace(' ', '_') + '_and_' + dataset2_name.replace(' ', '_') + '/'

## If it does not already exist, create directory that will contain the folder with PCA plots
if not os.path.exists(output_path):
    os.makedirs(output_path)

## If it does not already exist, create directory that will contain the PCA plots
if not os.path.exists(output_path2):
    os.makedirs(output_path2)
    
for i in range(1, 39):  
    ## For each chromosome number...
    ## Define path to input .vcf files for chromosome i of datasets 1 and 2
    path1 = PATH1.replace('*', str(i))
    path2 = PATH2.replace('*', str(i))
    
   ## Read input .vcf files for chromosome i of datasets 1 and 2
    data1 = read_vcf_file(path1)
    data2 = read_vcf_file(path2)
        
    ## Obtain SNPs data for chromosome i of datasets 1 and 2
    snps_1 = data1['calldata/GT']
    snps_2 = data2['calldata/GT']
        
    ## Convert SNPs from MxNx2 to MxN shape (averaged)
    # For dataset 1...
    length_1, num_dogs_1, num_strands_1 = snps_1.shape
    snps_1 = snps_1.reshape(length_1, num_dogs_1*2).T
    snps_1 = combine_chrom_strands(snps_1)
    # And for dataset 2...
    length_2, num_dogs_2, num_strands_2 = snps_2.shape
    snps_2 = snps_2.reshape(length_2, num_dogs_2*2).T
    snps_2 = combine_chrom_strands(snps_2)
        
    ## Concatenate SNPs for all chromosomes
    if i == 1:
        all_snps_1 = snps_1
        all_snps_2 = snps_2
    else:
        all_snps_1 = np.concatenate((all_snps_1, snps_1), axis=1)
        all_snps_2 = np.concatenate((all_snps_2, snps_2), axis=1)
        

## Make PCA plots... 

# 1. Trained and projected on the SNPs of the first dataset
PCA_2D_trained_and_projected_on_dataset(all_snps_1, plot_params, dataset1_name, color1, output_path2)

# 2. Trained and projected on the SNPs of the second dataset
PCA_2D_trained_and_projected_on_dataset(all_snps_2, plot_params, dataset2_name, color2, output_path2)

# 3. Trained on the SNPs of the first dataset and projected on both
PCA_2D_trained_on_dataset1_projected_on_both(all_snps_1, all_snps_2, plot_params, dataset1_name, dataset2_name, color1, color2, output_path2)

# 4. Trained on the SNPs of the second dataset and projected on both
PCA_2D_trained_on_dataset1_projected_on_both(all_snps_2, all_snps_1, plot_params, dataset2_name, dataset1_name, color2, color1, output_path2)

# 5. Trained and projected on the concatenation of the SNPs of both datasets
PCA_2D_trained_and_projected_on_both(all_snps_1, all_snps_2, plot_params, dataset1_name, dataset2_name, color1, color2, output_path2)
