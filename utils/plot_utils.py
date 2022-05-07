################################################################################
# Useful functions to plot a scatter plot with SNP means
# or 2D data from Principal Component Analysis (PCA)
################################################################################

import numpy as np
from typing import Dict
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


def plot_snp_means(all_query_snps: np.array, all_reference_snps: np.array, plot_dict: Dict, output_folder: str) -> None:
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
    
    # Save figure in output folder
    plt.savefig(f'{output_folder}snp_means_{plot_dict["x_axis_name"]}_and_{plot_dict["y_axis_name"]}', bbox_inches='tight')


def PCA_trained_and_projected_on_query(all_query_snps: np.array, plot_dict: Dict, output_folder: str) -> None:
    """
    Train and project Principal Component Analysis (PCA) on 
    the query SNPs data and plot a scatter plot with the 
    2D points. First, data is standardized to have 0 mean
    and 1 standard devision (std).
    
    Args:
        all_query_snps (List[str]): concatenation of SNPs in the query.
        plot_dict: dictionary with the configuration parameters of the plot.
        output_folder (str): folder to save the output.
        
    Returns:
        (None)
    
    """

    # Standardize the SNPs data to have 0 mean and 1 std
    snps_scaled = StandardScaler().fit_transform(all_query_snps)
    
    # Define PCA with two components
    pca = PCA(n_components=2)
    
    # Fit and transform PCA
    princ_comp = pca.fit_transform(snps_scaled)

    # Define plot figure
    plt.figure(figsize=(plot_dict['FIG_WIDTH'], plot_dict['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    # Plot the 2D PCA points
    plt.scatter(princ_comp[:,0], princ_comp[:,1], s=plot_dict['s'], alpha=plot_dict['alpha'], c=plot_dict["color_query"], label="Query")
    
    plt.title("PCA", fontsize = plot_dict['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_dict['FONTSIZE']})

    # Save figure in output folder
    plt.savefig(f'{output_folder}trained_and_projected_on_query', bbox_inches='tight')
    
    
def PCA_trained_on_query_projected_on_both(all_query_snps: np.array, all_reference_snps: np.array, plot_dict: Dict, output_folder: str) -> None:
    """
    Train Principal Component Analysis (PCA) on 
    the query SNPs data, project on both the query and the 
    reference and plot a scatter plot with the 2D points. First, data 
    is standardized to have 0 mean and 1 standard devision (std) 
    (based on the query).
    
    Args:
        all_query_snps (List[str]): concatenation of common markers in the query.
        all_reference_snps (List[str]): concatenation of common markers in the reference.
        plot_dict: dictionary with the configuration parameters of the plot.
        output_folder (str): folder to save the output.
        
    Returns:
        (None)
    
    """
    
    # Standardize SNPs data to have 0 mean and 1 std
    standardizer = StandardScaler().fit(all_query_snps)
    
    # Transform
    snps_reference = StandardScaler().transform(all_query_snps)
    snps_reference = StandardScaler().transform(all_reference_snps)
    
    # Define PCA with two components
    pca = PCA(n_components=2)

    # Fit PCA on the query
    pca = pca.fit(snps_query)

    # Project PCA on both the query and the reference
    princ_query = pca.transform(snps_query)
    princ_reference = pca.transform(snps_reference)

    # Define plot figure
    plt.figure(figsize=(plot_dict['FIG_WIDTH'], plot_dict['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    # Plot the query 2D PCA points
    plt.scatter(princ_query[:,0], princ_query[:,1], s=plot_dict['s'], alpha=plot_dict['alpha'], c=plot_dict["color_query"], label="query")
    
    # Plot the reference 2D PCA points
    plt.scatter(princ_reference[:,0], princ_reference[:,1], s=plot_dict['s'], alpha=plot_dict['alpha'], c=plot_dict["color_reference"], label="reference")
    
    plt.title("PCA trained on the query projected on both", fontsize = plot_dict['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_dict['FONTSIZE']})

    # Save figure in output folder
    plt.savefig(f'{output_folder}trained_on_query_projected_on_both', bbox_inches='tight')
    

def PCA_trained_and_projected_on_both(all_query_snps: np.array, all_reference_snps: np.array, plot_dict: Dict, output_folder: str) -> None:
    """
    Train Principal Component Analysis (PCA) on 
    the query and reference SNPs data and plot a scatter plot with the 2D points. 
    First, data is standardized to have 0 mean and 1 standard devision (std).
    
    Args:
        all_query_snps (List[str]): concatenation of common markers in the query.
        all_reference_snps (List[str]): concatenation of common markers in the reference.
        plot_dict: dictionary with the configuration parameters of the plot.
        output_folder (str): folder to save the output.
        
    Returns:
        (None)
    
    """
    
    # Concatenate SNPs data
    concat = np.concatenate((all_query_snps, all_reference_snps), axis=0, out=None)

    # Standardize SNPs data to have 0 mean and 1 std
    concat_scaled = StandardScaler().fit_transform(concat)

    # Define PCA with two components
    pca = PCA(n_components=2)

    # Fit and transform PCA model on both the query and the reference
    pca = pca.fit(concat_scaled)
    princ_comp_concat = pca.transform(concat_scaled)

    # Define plot figure
    plt.figure(figsize=(plot_dict['FIG_WIDTH'], plot_dict['FIG_HEIGHT']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    
    # Plot the query 2D PCA points
    plt.scatter(princ_comp_concat[:all_query_snps.shape[0],0], princ_comp_concat[:all_query_snps.shape[0],1], s=plot_dict['s'], 
                alpha=plot_dict['alpha'], c=plot_dict["color_query"], label="query")
    
    # Plot the reference 2D PCA points
    plt.scatter(princ_comp_concat[all_query_snps.shape[0]:,0], princ_comp_concat[all_query_snps.shape[0]:,1], s=plot_dict['s'], 
                alpha=plot_dict['alpha'], c=color2, label="reference")
    
    plt.title("PCA trained and projecred on both the query the reference", fontsize = plot_dict['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_dict['FONTSIZE']})

    # Save figure in output folder
    plt.savefig(f'{output_folder}trained_and_projected_on_both', bbox_inches='tight')
