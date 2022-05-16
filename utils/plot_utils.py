################################################################################
# Useful functions to plot a scatter plot with SNP means
# or 2D data from Principal Component Analysis (PCA)
################################################################################

import numpy as np
from typing import Dict, List
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import logging


def snp_means_plot(query_snps: np.array, reference_snps: np.array, 
                   plot_dict: Dict, output_folder: str, logger: logging.Logger) -> None:
    """
    Generates a scatter plot with the SNP means of the common markers 
    between the query and the reference datasets. Stores the plot in
    output folder.
    
    Args:
        query_snps (List[str]): concatenation of all SNPs in the query
        that are also present in the reference.
        reference_snps (List[str]): concatenation of all SNPs in the 
        reference that are also present in the query.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        (None)
    
    """
  
    # Compute the mean of each SNP in the query
    mean_query = np.mean(query_snps, axis=0)
    
    # Compute the mean of each SNP in the reference
    mean_reference = np.mean(reference_snps, axis=0)
    
    # Define plot figure
    plt.figure(figsize=(plot_dict['fig_width'], plot_dict['fig_height']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams.update({'font.size': plot_dict['fontsize']})
    
    # Plot SNP means
    plt.scatter(mean_query, mean_reference, s=plot_dict['s'], color=plot_dict['color'])
    
    # Define plot title, x and y axis
    plt.title('SNP means')
    plt.xlabel(f'\nSNP mean in {plot_dict["x_axis_name"]}')
    plt.ylabel(f'\nSNP mean in {plot_dict["y_axis_name"]}')
    plt.xticks(ticks = [0, .2, .4, .6, .8, 1.])
    plt.yticks(ticks = [0, .2, .4, .6, .8, 1.])
    
    # Define output name to .png image with SNP means plot 
    # The .png image will have the name snp_means_{x_axis_name}_{y_axis_name}.png
    # Example: snp_means_array_embark.png
    output_name = f'snp_means_{plot_dict["x_axis_name"]}_and_{plot_dict["y_axis_name"]}.png'
    
    # Save figure in output folder
    logger.debug(f'Saving plot with SNP means in {output_folder}{output_name}.')
    plt.savefig(output_folder+output_name, bbox_inches='tight')


def PCA_trained_and_projected_on_query(query_snps: np.array, plot_dict: Dict, 
                                       output_folder: str) -> None:
    """
    Train and project Principal Component Analysis (PCA) on 
    the query SNPs data and plot a scatter plot with the 
    2D points. First, data is standardized to have 0 mean
    and 1 standard devision (std).
    
    Args:
        query_snps (List[str]): concatenation of SNPs in the query.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        
    Returns:
        (None)
    
    """

    # Standardize the SNPs data to have 0 mean and 1 std
    snps_scaled = StandardScaler().fit_transform(query_snps)
    
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
    
    
def PCA_trained_on_query_projected_on_both(query_snps: np.array, reference_snps: np.array, 
                                           plot_dict: Dict, output_folder: str) -> None:
    """
    Train Principal Component Analysis (PCA) on 
    the query SNPs data, project on both the query and the 
    reference and plot a scatter plot with the 2D points. First, data 
    is standardized to have 0 mean and 1 standard devision (std) 
    (based on the query).
    
    Args:
        query_snps (List[str]): concatenation of common markers in the query.
        reference_snps (List[str]): concatenation of common markers in the reference.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        
    Returns:
        (None)
    
    """
    
    # Standardize SNPs data to have 0 mean and 1 std
    standardizer = StandardScaler().fit(query_snps)
    
    # Transform
    snps_reference = StandardScaler().transform(query_snps)
    snps_reference = StandardScaler().transform(reference_snps)
    
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
    

def PCA_trained_and_projected_on_both(query_snps: np.array, reference_snps: np.array, 
                                      plot_dict: Dict, output_folder: str) -> None:
    """
    Train Principal Component Analysis (PCA) on 
    the query and reference SNPs data and plot a scatter plot with the 2D points. 
    First, data is standardized to have 0 mean and 1 standard devision (std).
    
    Args:
        query_snps (List[str]): concatenation of common markers in the query.
        reference_snps (List[str]): concatenation of common markers in the reference.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        
    Returns:
        (None)
    
    """
    
    # Concatenate SNPs data
    concat = np.concatenate((query_snps, reference_snps), axis=0, out=None)

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
    plt.scatter(princ_comp_concat[:query_snps.shape[0],0], princ_comp_concat[:query_snps.shape[0],1], s=plot_dict['s'], 
                alpha=plot_dict['alpha'], c=plot_dict["color_query"], label="query")
    
    # Plot the reference 2D PCA points
    plt.scatter(princ_comp_concat[query_snps.shape[0]:,0], princ_comp_concat[query_snps.shape[0]:,1], s=plot_dict['s'], 
                alpha=plot_dict['alpha'], c=color2, label="reference")
    
    plt.title("PCA trained and projecred on both the query the reference", fontsize = plot_dict['FONTSIZE'])
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['FONTSIZE'])
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['FONTSIZE'])
    plt.legend(loc='upper right', prop={'size': plot_dict['FONTSIZE']})

    # Save figure in output folder
    plt.savefig(f'{output_folder}trained_and_projected_on_both', bbox_inches='tight')
