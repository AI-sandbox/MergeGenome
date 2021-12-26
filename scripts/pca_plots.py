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
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from umap import UMAP
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')
from utils.vcf_utils import read_vcf_file, combine_chrom_strands
from utils.track import track

################################################################################

## USER INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define base path to .vcf files for preprocessed data for each chromosome
PATH1 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/embark/embark_chr*.vcf'
PATH2 = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array_imputed_subset/imputed_array_subset_chr*.vcf'

## Define name of datasets 1 and 2
dataset1_name = 'Embark'
dataset2_name = 'Imputed Array Subsetted'

## Define path to directory that will contain the folder with the PCA plots
output_path = '/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/output/PCA_plots/'

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
        

################################################################################

## Standardize the SNPs data to have 0 mean and 1 std
# For dataset 1...
all_snps_1_scaled = StandardScaler().fit_transform(all_snps_1)
# And for dataset 2...
all_snps_2_scaled = StandardScaler().fit_transform(all_snps_2)

## Define fontsize, figure width and height of PCA plots
FONTSIZE = 16
FIG_WIDTH = 12
FIG_HEIGHT = 10

# ------------------- Make PCA plot trained and projected on the SNPs of the first dataset -------------------

## Define PCA object with 2 components
pca = PCA(n_components=2)

## Fit and transform the PCA model on the SNPs of the first dataset
principal_components = pca.fit_transform(all_snps_1_scaled)

## Plot the result
plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.scatter(principal_components[:,0], principal_components[:,1], s=10, c='#259988', alpha=0.3)
plt.title("{} PCA".format(dataset1_name), fontsize = FONTSIZE)
plt.xlabel('\nPrincipal Component 1', fontsize = FONTSIZE)
plt.ylabel('\nPrincipal Component 2', fontsize = FONTSIZE)

## Save figure in output path
plt.savefig(output_path+'trained_and_projected_on_{}'.format(dataset1_name.replace(' ', '_')),  bbox_inches='tight')

# ------------------- Make PCA plot trained and projected on the SNPs of the first dataset -------------------

## Define PCA object with 2 components
pca = PCA(n_components=2)

## Fit and transform the PCA model on the SNPs of the second dataset
principal_components = pca.fit_transform(all_snps_2_scaled)

## Plot the result
plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.scatter(principal_components[:,0], principal_components[:,1], s=10, c='#EBD0A1', alpha=0.3)
plt.title("{} PCA".format(dataset2_name), fontsize = FONTSIZE)
plt.xlabel('\nPrincipal Component 1', fontsize = FONTSIZE)
plt.ylabel('\nPrincipal Component 2', fontsize = FONTSIZE)

## Save figure in output path
plt.savefig(output_path+'trained_and_projected_on_{}'.format(dataset2_name.replace(' ', '_')),  bbox_inches='tight')

# ------------------- Make PCA plot trained on the SNPs of the first dataset and projected on both -------------------

## Define PCA object with 2 components
pca = PCA(n_components=2)

## Fit the PCA model on the SNPs of the first dataset
pca = pca.fit(all_snps_1_scaled)

## Transform on the SNPs of the first and second dataset (separately)
principal_components1 = pca.transform(all_snps_1_scaled)
principal_components2 = pca.transform(all_snps_2_scaled)

## Concatenate results
principal_components_concat = np.concatenate((principal_components1, principal_components2), axis=0, out=None)

## Plot the result
plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.scatter(principal_components_concat[:,0], principal_components_concat[:,1], s=10, alpha=0.3, c=['#259988']*principal_components1.shape[0] + ['#EBD0A1']*principal_components2.shape[0])
plt.title("Trained on {} projected on {} and {} PCA".format(dataset1_name, dataset1_name, dataset2_name), fontsize = FONTSIZE)
plt.xlabel('\nPrincipal Component 1', fontsize = FONTSIZE)
plt.ylabel('\nPrincipal Component 2', fontsize = FONTSIZE)

## Save figure in output path
plt.savefig(output_path+'trained_on_{}_projected_on_{}_and_{}'.format(dataset1_name.replace(' ', '_'), 
                                                                      dataset1_name.replace(' ', '_'), 
                                                                      dataset2_name.replace(' ', '_'), bbox_inches='tight'))

# ------------------- Make PCA plot trained on the SNPs of the first dataset and projected on both -------------------

## Define PCA object with 2 components
pca = PCA(n_components=2)

## Fit the PCA model on the SNPs of the second dataset
pca = pca.fit(all_snps_2_scaled)

## Transform on the SNPs of the first and second dataset (separately)
principal_components1 = pca.transform(all_snps_1_scaled)
principal_components2 = pca.transform(all_snps_2_scaled)

## Concatenate results
principal_components_concat = np.concatenate((principal_components1, principal_components2), axis=0, out=None)

## Plot the result
plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.scatter(principal_components_concat[:,0], principal_components_concat[:,1], s=10, c=['#259988']*principal_components1.shape[0] + ['#EBD0A1']*principal_components2.shape[0], alpha=0.3)
plt.title("Trained on {} projected on {} and {} PCA".format(dataset2_name, dataset1_name, dataset2_name), fontsize = FONTSIZE)
plt.xlabel('\nPrincipal Component 1', fontsize = FONTSIZE)
plt.ylabel('\nPrincipal Component 2', fontsize = FONTSIZE)

## Save figure in output path
plt.savefig(output_path+'trained_on_{}_projected_on_{}_and_{}'.format(dataset2_name.replace(' ', '_'), 
                                                                      dataset1_name.replace(' ', '_'), 
                                                                      dataset2_name.replace(' ', '_'), bbox_inches='tight'))

# ------------------- Make PCA plot trained and projected on the SNPs of both datasets -------------------

## Concatenate the SNPs data
concat = np.concatenate((all_snps_1, all_snps_2), axis=0, out=None)

## Standardize the SNPs data to have 0 mean and 1 std
concat_scaled = StandardScaler().fit_transform(concat)

## Define PCA object with 2 components
pca = PCA(n_components=2)

## Fit and transform the PCA model on the SNPs of both datasets
pca = pca.fit(concat_scaled)
principal_components_concat = pca.transform(concat_scaled)

## Plot the result
plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.scatter(principal_components_concat[:,0], principal_components_concat[:,1], s=10, c=['#259988']*all_snps_1.shape[0] + ['#EBD0A1']*all_snps_2.shape[0], alpha=0.3)
plt.title("Trained and projected on {} and {} PCA".format(dataset1_name, dataset2_name), fontsize = FONTSIZE)
plt.xlabel('\nPrincipal Component 1', fontsize = FONTSIZE)
plt.ylabel('\nPrincipal Component 2', fontsize = FONTSIZE)

## Save figure in output path
plt.savefig(output_path+'trained_on_{}_projected_on_{}_and_{}'.format(dataset2_name.replace(' ', '_'), 
                                                                      dataset1_name.replace(' ', '_'), 
                                                                      dataset2_name.replace(' ', '_'), bbox_inches='tight'))