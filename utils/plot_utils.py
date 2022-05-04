################################################################################
# Useful functions for PCA, UMAP and tSNE plots
################################################################################

import numpy as np
from typing import Dict
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


def plot_snp_means(all_query_snps: np.array, all_reference_snps: np.array, plot_dict: Dict, output_folder: str):
    """
    Plot a scatter plot with the SNP means of the sorted 
    common markers between two datasets.
    
    Args:
        all_query_snps (List[str]): concatenation of common markers in the query.
        all_reference_snps (List[str]): concatenation of common markers in the reference.
        plot_dict: dictionary with the configuration parameters of the plot.
        output_folder (str): folder to save the output.
        
    Returns:
        (None)
    
    """
    
    # Ensure all_query_snps and all_reference_snps contain the same number of SNPs
    assert all_query_snps.shape[1] == all_reference_snps.shape[1], 'The number of SNPs in the query and reference do not match.' \
    'Compute the common markers.'
    
    # Compute the mean of each SNP
    mean_query = np.mean(all_query_snps, axis=0)
    mean_reference = np.mean(all_reference_snps, axis=0)
    
    # Define plot figure
    plt.figure(figsize=(plot_dict['FIG_WIDTH'], plot_dict['FIG_HEIGHT']))
    plt.scatter(mean_query, mean_reference, s=plot_dict['s'], alpha=plot_dict['alpha'], color=plot_dict['color'])
    plt.title('SNP means', fontsize = plot_dict['FONTSIZE'])
    plt.xlabel(f'\nSNP mean in {plot_dict["x_axis_name"]}', fontsize = plot_dict['FONTSIZE'])
    plt.ylabel(f'\nSNP mean in {plot_dict["y_axis_name"]}', fontsize = plot_dict['FONTSIZE'])
    plt.xticks(ticks = [0, .2, .4, .6, .8, 1.], fontsize = plot_dict['FONTSIZE']-3)
    plt.yticks(ticks = [0, .2, .4, .6, .8, 1.], fontsize = plot_dict['FONTSIZE']-3)
    
    # Save figure in output path
    plt.savefig(f'{output_path}snp_means_{plot_dict["x_axis_name"]}_and_{plot_dict["y_axis_name"]}', bbox_inches='tight')


def PCA_2D_trained_and_projected_on_dataset(snps, plot_dict, dataset_name, color, output_path):
    '''
    Objective: make PCA plot trained and projected on snps.
    Input:
        - snps: SNPs data.
        - plot_dict: dictionary containing plot parameters.
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
    plt.figure(figsize=(plot_dict['FIG_WIDTH'], plot_dict['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    plt.scatter(princ_comp[:,0], princ_comp[:,1], s=plot_dict['s'], alpha=plot_dict['alpha'], c=color, label=dataset_name)
    
    plt.title("{} PCA".format(dataset_name), fontsize = plot_dict['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_dict['FONTSIZE']})

    ## Save figure in output path
    plt.savefig(output_path+'trained_and_projected_on_{}'.format(dataset_name.replace(' ', '_')),  bbox_inches='tight')
    
    
def PCA_2D_trained_on_dataset1_projected_on_both(all_query_snps, all_reference_snps, plot_dict, dataset1_name, dataset2_name, color1, color2, output_path):
    '''
    Objective: make PCA plot trained on all_query_snps and projected on both all_query_snps and all_reference_snps.
    Input:
        - all_query_snps: SNPs data of dataset 1.
        - all_reference_snps: SNPs data of dataset 2.
        - plot_dict: dictionary containing plot parameters.
        - dataset1_name: name given to dataset 1 and that will appear in the title.
        - dataset2_name: name given to dataset 2 and that will appear in the title.
        - color1: of the data points of dataset 1.
        - color2: of the data points of dataset 1.
        - output_path: path to output that will save the plots.
    '''
    
    ## Standardize the SNPs data to have 0 mean and 1 std
    snps_scaled1 = StandardScaler().fit_transform(all_query_snps)
    snps_scaled2 = StandardScaler().fit_transform(all_reference_snps)

    ## Define PCA object with 2 components
    pca = PCA(n_components=2)

    ## Fit on all_query_snps
    pca = pca.fit(snps_scaled1)

    ## Projecte on both all_query_snps and all_reference_snps
    princ_comp1 = pca.transform(snps_scaled1)
    princ_comp2 = pca.transform(snps_scaled2)

    ## Define plot figure
    plt.figure(figsize=(plot_dict['FIG_WIDTH'], plot_dict['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    plt.scatter(princ_comp1[:,0], princ_comp1[:,1], s=plot_dict['s'], alpha=plot_dict['alpha'], c=color1, label=dataset1_name)
    plt.scatter(princ_comp2[:,0], princ_comp2[:,1], s=plot_dict['s'], alpha=plot_dict['alpha'], c=color2, label=dataset2_name)
    
    plt.title("Trained on {} projected on {} and {} PCA".format(dataset1_name, dataset1_name, dataset2_name), fontsize = plot_dict['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_dict['FONTSIZE']})

    ## Save figure in output path
    plt.savefig(output_path+'trained_on_{}_projected_on_{}_and_{}'.format(dataset1_name.replace(' ', '_'), dataset1_name.replace(' ', '_'), 
                                                                          dataset2_name.replace(' ', '_'), bbox_inches='tight'))


def PCA_2D_trained_and_projected_on_both(all_query_snps, all_reference_snps, plot_dict, dataset1_name, dataset2_name, color1, color2, output_path):
    '''
    Objective: make PCA plot trained and projected on concatenation of all_query_snps and all_reference_snps.
    Input:
        - all_query_snps: SNPs data of dataset 1.
        - all_reference_snps: SNPs data of dataset 2.
        - plot_dict: dictionary containing plot parameters.
        - dataset1_name: name given to dataset 1 and that will appear in the title.
        - dataset2_name: name given to dataset 2 and that will appear in the title.
        - color1: of the data points of dataset 1.
        - color2: of the data points of dataset 1.
        - output_path: path to output that will save the plots.
    '''
    
    ## Concatenate the SNPs data
    concat = np.concatenate((all_query_snps, all_reference_snps), axis=0, out=None)

    ## Standardize the SNPs data to have 0 mean and 1 std
    concat_scaled = StandardScaler().fit_transform(concat)

    ## Define PCA object with 2 components
    pca = PCA(n_components=2)

    ## Fit and transform the PCA model on the SNPs of both datasets
    pca = pca.fit(concat_scaled)
    princ_comp_concat = pca.transform(concat_scaled)

    ## Define plot figure
    plt.figure(figsize=(plot_dict['FIG_WIDTH'], plot_dict['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    plt.scatter(princ_comp_concat[:all_query_snps.shape[0],0], princ_comp_concat[:all_query_snps.shape[0],1], s=plot_dict['s'], alpha=plot_dict['alpha'], c=color1, label=dataset1_name)
    plt.scatter(princ_comp_concat[all_query_snps.shape[0]:,0], princ_comp_concat[all_query_snps.shape[0]:,1], s=plot_dict['s'], alpha=plot_dict['alpha'], c=color2, label=dataset2_name)
    
    plt.title("Trained and projected on {} and {} PCA".format(dataset1_name, dataset2_name), fontsize = plot_dict['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_dict['FONTSIZE']})

    ## Save figure in output path
    plt.savefig(output_path+'trained_and_projected_on_{}_and_{}'.format(dataset1_name.replace(' ', '_'), dataset2_name.replace(' ', '_'), bbox_inches='tight'))
