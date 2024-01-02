################################################################################
# Useful functions to plot a scatter plot with: 
# * SNP means comparison.
# * First two components of Principal Component Analysis (PCA).
# For the common markers between the query and reference.
################################################################################

import numpy as np
from typing import Dict, List
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import logging


def snp_means_plot(query_snps: np.array, reference_snps: np.array, indexes: np.array,
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
        indexes (List[int]): indexes of the SNPs that will be plotted in a 
        different color.
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
    
    # Define vector with all the SNP indexes from [0, n_snps-1]
    snp_idxs = list(range(len(mean_query)))
    
    # Define vector with all SNP indexes that are not in indexes
    not_indexes = list(set(snp_idxs) - set(indexes))
    
    # Define plot figure
    plt.figure(figsize=(plot_dict['fig_width'], plot_dict['fig_height']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams.update({'font.size': plot_dict['fontsize']})
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.tick_params('both', length=16, width=2, which='major')
    plt.tick_params(axis='x', which='major', pad=8)
    plt.tick_params(axis='y', which='major', pad=8)
    plt.rcParams.update({'patch.linewidth' : '1.5'})
    for spine in ['top', 'right', 'bottom', 'left']:
        plt.gca().spines[spine].set_linewidth(2)
    
    # Plot SNP means
    plt.scatter(mean_reference[indexes], mean_query[indexes], s=plot_dict['s'], alpha=plot_dict['alpha'],
                color=plot_dict['color_points_indexes'], label=plot_dict['legend_points_indexes'])
    plt.scatter(mean_reference[not_indexes], mean_query[not_indexes], s=plot_dict['s'], alpha=plot_dict['alpha'],
                color=plot_dict['color_points'], label=plot_dict['legend_points'])
    
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Define x and y axis
    plt.xlabel(f'\nSNP mean in {plot_dict["x_axis_name"]}', fontsize = plot_dict['fontsize'], labelpad=-30)
    plt.ylabel(f'\nSNP mean in {plot_dict["y_axis_name"]}', fontsize = plot_dict['fontsize'])
    plt.xticks(ticks = [0, .2, .4, .6, .8, 1.])
    plt.yticks(ticks = [0, .2, .4, .6, .8, 1.])
    
    # Add legend
    plt.legend(loc='lower right', markerscale=20)
    plt.grid()
    
    # Define output name to .png image with SNP means plot 
    # The .png image will have the name snp_means_{x_axis_name}_{y_axis_name}.png
    # Example: snp_means_array_embark.png
    output_name = f'snp_means_{plot_dict["x_axis_name"]}_and_{plot_dict["y_axis_name"]}.png'
    
    # Save figure in output folder
    logger.debug(f'Saving plot with SNP means in {output_folder}{output_name}.')
    plt.savefig(output_folder+output_name, bbox_inches='tight', dpi=300)
    plt.close()
    

def PCA_trained_and_projected_on_query(query_snps: np.array, plot_dict: Dict, output_folder: str, 
                                       logger: logging.Logger) -> None:
    """
    Applies Principal Component Analysis (PCA) on the query and genenrates a 
    scatter plot with the first two components. First, data is standardized to 
    have 0 mean and 1 standard devision (std).
    
    Args:
        query_snps (List[str]): concatenation of all SNPs in the query.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        
    Returns:
        (None)
    
    """
    # Standardize SNPs to have 0 mean and 1 std
    snps_scaled = StandardScaler().fit_transform(query_snps)
    
    # Define PCA with two components
    pca = PCA(n_components=2)
    
    # Fit and transform PCA
    princ_comp = pca.fit_transform(snps_scaled)

    # Define plot figure
    plt.figure(figsize=(plot_dict['fig_width'], plot_dict['fig_height']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams.update({'font.size': plot_dict['fontsize']})
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.tick_params('both', length=16, width=2, which='major')
    plt.tick_params(axis='x', which='major', pad=8)
    plt.tick_params(axis='y', which='major', pad=8)
    plt.rcParams.update({'patch.linewidth' : '1.5'})
    for spine in ['top', 'right', 'bottom', 'left']:
        plt.gca().spines[spine].set_linewidth(2)
    
    # Plot the first two components
    plt.scatter(princ_comp[:,0], princ_comp[:,1], s=plot_dict['s'], c=plot_dict["color_query"], label=plot_dict['legend_query'])
    
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Define x and y axis, and legend
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['fontsize'], labelpad=-30)
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['fontsize'])
    plt.legend(loc='upper right', markerscale=20)

    # Define output name to .png image with SNP means plot 
    # The .png image will have the name 'trained_query.png'
    output_name = 'trained_query.png'
    
    # Save figure in output folder
    logger.debug(f'Saving plot with SNP means in {output_folder}{output_name}.')
    plt.savefig(output_folder+output_name, bbox_inches='tight', dpi=300)
    plt.close()
    
    
def PCA_trained_on_query_projected_on_both(query_snps: np.array, reference_snps: np.array, plot_dict: Dict,
                                           output_folder: str, logger: logging.Logger) -> None:
    """
    Applies Principal Component Analysis (PCA) on the query and obtains
    the first two components of the projected query and reference. Genenrates a 
    scatter plot with the first two components. First, data is standardized to 
    have 0 mean and 1 standard devision (std) based on the query.
    
    Args:
        query_snps (List[str]): concatenation of all SNPs in the query
        that are also present in the reference.
        reference_snps (List[str]): concatenation of all SNPs in the 
        reference that are also present in the query.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        
    Returns:
        (None)
    
    """
    # Standardize SNPs to have 0 mean and 1 std
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
    plt.figure(figsize=(plot_dict['fig_width'], plot_dict['fig_height']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams.update({'font.size': plot_dict['fontsize']})
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.tick_params('both', length=16, width=2, which='major')
    plt.tick_params(axis='x', which='major', pad=8)
    plt.tick_params(axis='y', which='major', pad=8)
    plt.rcParams.update({'patch.linewidth' : '1.5'})
    for spine in ['top', 'right', 'bottom', 'left']:
        plt.gca().spines[spine].set_linewidth(2)
    
    # Plot the reference 2D PCA points
    plt.scatter(princ_reference[:,0], princ_reference[:,1], s=plot_dict['s'], 
                alpha=plot_dict['alpha'], c=plot_dict["color_reference"], label=plot_dict['legend_reference'])
    
    # Plot the query 2D PCA points
    plt.scatter(princ_query[:,0], princ_query[:,1], s=plot_dict['s'], 
                alpha=plot_dict['alpha'], c=plot_dict["color_query"], label=plot_dict['legend_query'])
    
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Define x and y axis, and legend
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['fontsize'], labelpad=-30)
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['fontsize'])
    plt.legend(loc='upper right', markerscale=1)

    # Define output name to .png image with SNP means plot 
    # The .png image will have the name 'trained_query_projected_both.png'
    output_name = 'trained_query_projected_both.png'
    
    # Save figure in output folder
    logger.debug(f'Saving plot with SNP means in {output_folder}{output_name}.')
    plt.savefig(output_folder+output_name, bbox_inches='tight', dpi=300)
    plt.close()


def PCA_trained_and_projected_on_both(query_snps: np.array, reference_snps: np.array, plot_dict: Dict, 
                                      output_folder: str, logger: logging.Logger) -> None:
    """
    Train Principal Component Analysis (PCA) on 
    the query and reference SNPs data and plot a scatter plot with the 2D points. 
    First, data is standardized to have 0 mean and 1 standard devision (std).
    
    Args:
        query_snps (List[str]): concatenation of all SNPs in the query
        that are also present in the reference.
        reference_snps (List[str]): concatenation of all SNPs in the 
        reference that are also present in the query.
        plot_dict[Dict]: configuration parameters of the plot.
        output_folder (str): path to output folder.
        
    Returns:
        (None)
    
    """
    # Concatenate SNPs from the query and the reference
    concat = np.concatenate((query_snps, reference_snps), axis=0, out=None)

    # Standardize SNPs to have 0 mean and 1 std
    concat_scaled = StandardScaler().fit_transform(concat)

    # Define PCA with two components
    pca = PCA(n_components=2)

    # Fit and transform PCA model on both the query and the reference
    pca = pca.fit(concat_scaled)
    princ_comp_concat = pca.transform(concat_scaled)
    
    # Define plot figure
    plt.figure(figsize=(plot_dict['fig_width'], plot_dict['fig_height']))
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams.update({'font.size': plot_dict['fontsize']})
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.tick_params('both', length=16, width=2, which='major')
    plt.tick_params(axis='x', which='major', pad=8)
    plt.tick_params(axis='y', which='major', pad=8)
    plt.rcParams.update({'patch.linewidth' : '1.5'})
    for spine in ['top', 'right', 'bottom', 'left']:
        plt.gca().spines[spine].set_linewidth(2)
    
    # Plot the reference 2D PCA points
    plt.scatter(princ_comp_concat[query_snps.shape[0]:,0], princ_comp_concat[query_snps.shape[0]:,1], 
                s=plot_dict['s'], alpha=plot_dict['alpha'], c=plot_dict["color_reference"], label=plot_dict['legend_reference'])
    
    # Plot the query 2D PCA points
    plt.scatter(princ_comp_concat[:query_snps.shape[0],0], princ_comp_concat[:query_snps.shape[0],1], 
                s=plot_dict['s'], alpha=plot_dict['alpha'], c=plot_dict["color_query"], label=plot_dict['legend_query'])
    
    # Define x and y axis, and legend
    plt.xlabel('\nPrincipal Component 1', fontsize = plot_dict['fontsize'], labelpad=-30)
    plt.ylabel('\nPrincipal Component 2', fontsize = plot_dict['fontsize'])
    plt.legend(loc='upper right', markerscale=3.5)
    plt.grid()
    
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    # Define output name to .png image with SNP means plot 
    # The .png image will have the name 'trained_both_projected_both.png'
    output_name = 'trained_both_projected_both.png'
    
    # Save figure in output folder
    logger.debug(f'Saving plot with SNP means in {output_folder}{output_name}.')
    plt.savefig(output_folder+output_name, bbox_inches='tight', dpi=300)
    plt.close()
