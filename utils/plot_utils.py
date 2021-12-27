################################################################################
# Useful functions for PCA, UMAP and tSNE plots
################################################################################

import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from umap import UMAP
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt


def PCA_2D_trained_and_projected_on_dataset(snps, plot_params, dataset_name, color, output_path):
    '''
    Objective: make PCA plot trained and projected on snps.
    Input:
        - snps: SNPs data.
        - plot_params: dictionary containing plot parameters.
        - dataset_name: name given to dataset and that will appear in the title.
        - color: of the data points.
        - output_path: path to output that will save the plots.
    '''

    ## Standardize the SNPs data to have 0 mean and 1 std
    snps_scaled = StandardScaler().fit_transform(snps)
    
    ## Define PCA object with 2 components
    pca = PCA(n_components=2)
    
    ## Fit and transform on snps
    princ_comp = pca.fit_transform(snps_scaled)

    ## Define plot figure
    plt.figure(figsize=(plot_params['FIG_WIDTH'], plot_params['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    plt.scatter(princ_comp[:,0], princ_comp[:,1], s=plot_params['s'], alpha=plot_params['alpha'], c=color, label=dataset_name)
    
    plt.title("{} PCA".format(dataset_name), fontsize = plot_params['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_params['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_params['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_params['FONTSIZE']})

    ## Save figure in output path
    plt.savefig(output_path+'trained_and_projected_on_{}'.format(dataset_name.replace(' ', '_')),  bbox_inches='tight')
    
    
def PCA_2D_trained_on_dataset1_projected_on_both(snps1, snps2, plot_params, dataset1_name, dataset2_name, color1, color2, output_path):
    '''
    Objective: make PCA plot trained on snps1 and projected on both snps1 and snps2.
    Input:
        - snps1: SNPs data of dataset 1.
        - snps2: SNPs data of dataset 2.
        - plot_params: dictionary containing plot parameters.
        - dataset1_name: name given to dataset 1 and that will appear in the title.
        - dataset2_name: name given to dataset 2 and that will appear in the title.
        - color1: of the data points of dataset 1.
        - color2: of the data points of dataset 1.
        - output_path: path to output that will save the plots.
    '''
    
    ## Standardize the SNPs data to have 0 mean and 1 std
    snps_scaled1 = StandardScaler().fit_transform(snps1)
    snps_scaled2 = StandardScaler().fit_transform(snps2)

    ## Define PCA object with 2 components
    pca = PCA(n_components=2)

    ## Fit on snps1
    pca = pca.fit(snps_scaled1)

    ## Projecte on both snps1 and snps2
    princ_comp1 = pca.transform(snps_scaled1)
    princ_comp2 = pca.transform(snps_scaled2)

    ## Define plot figure
    plt.figure(figsize=(plot_params['FIG_WIDTH'], plot_params['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    plt.scatter(princ_comp1[:,0], princ_comp1[:,1], s=plot_params['s'], alpha=plot_params['alpha'], c=color1, label=dataset1_name)
    plt.scatter(princ_comp2[:,0], princ_comp2[:,1], s=plot_params['s'], alpha=plot_params['alpha'], c=color2, label=dataset2_name)
    
    plt.title("Trained on {} projected on {} and {} PCA".format(dataset1_name, dataset1_name, dataset2_name), fontsize = plot_params['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_params['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_params['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_params['FONTSIZE']})

    ## Save figure in output path
    plt.savefig(output_path+'trained_on_{}_projected_on_{}_and_{}'.format(dataset1_name.replace(' ', '_'), dataset1_name.replace(' ', '_'), 
                                                                          dataset2_name.replace(' ', '_'), bbox_inches='tight'))


def PCA_2D_trained_and_projected_on_both(snps1, snps2, plot_params, dataset1_name, dataset2_name, color1, color2, output_path):
    '''
    Objective: make PCA plot trained and projected on concatenation of snps1 and snps2.
    Input:
        - snps1: SNPs data of dataset 1.
        - snps2: SNPs data of dataset 2.
        - plot_params: dictionary containing plot parameters.
        - dataset1_name: name given to dataset 1 and that will appear in the title.
        - dataset2_name: name given to dataset 2 and that will appear in the title.
        - color1: of the data points of dataset 1.
        - color2: of the data points of dataset 1.
        - output_path: path to output that will save the plots.
    '''
    
    ## Concatenate the SNPs data
    concat = np.concatenate((all_snps_1, all_snps_2), axis=0, out=None)

    ## Standardize the SNPs data to have 0 mean and 1 std
    concat_scaled = StandardScaler().fit_transform(concat)

    ## Define PCA object with 2 components
    pca = PCA(n_components=2)

    ## Fit and transform the PCA model on the SNPs of both datasets
    pca = pca.fit(concat_scaled)
    princ_comp_concat = pca.transform(concat_scaled)

    ## Define plot figure
    plt.figure(figsize=(plot_params['FIG_WIDTH'], plot_params['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    plt.scatter(princ_comp_concat[:snps1.shape[0],0], princ_comp_concat[:snps1.shape[0],1], s=plot_params['s'], alpha=plot_params['alpha'], c=color1, label=dataset1_name)
    plt.scatter(princ_comp_concat[snps1.shape[0]:,0], princ_comp_concat[snps1.shape[0]:,1], s=plot_params['s'], alpha=plot_params['alpha'], c=color2, label=dataset2_name)
    
    plt.title("Trained and projected on {} and {} PCA".format(dataset1_name, dataset2_name), fontsize = plot_params['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_params['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_params['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_params['FONTSIZE']})

    ## Save figure in output path
    plt.savefig(output_path2+'trained_and_projected_on_{}_and_{}'.format(dataset1_name.replace(' ', '_'), dataset2_name.replace(' ', '_'), bbox_inches='tight'))
