################################################################################
# Reading/writing .vcf files and performing small transformations
################################################################################

import os
import allel
import numpy as np
import pandas as pd

from utils.track import track

def read_vcf_file(path_to_vcf_file):
    '''
    Objective:
        - Read .vcf file and return the content of the .vcf file.
    Input:
        - path_to_vcf_file: path to .vcf file.
    Output:
        - vcf_data: dictionary with the content of .vcf file.
    '''
    
    ## Check the path to the .vcf file exists
    assert os.path.isfile(path_to_vcf_file), f'{path_to_vcf_file} not found. Is it a file?'
    
    ## Read the .vcf file using skikit-allel
    vcf_data = allel.read_vcf(path_to_vcf_file)
    
    return vcf_data


def write_vcf_file(vcf_data, output_path_to_vcf_file):
    
    '''
    Objective:
        - Write SNPs data in .vcf file.
    Input:
        - vcf_data: modified allel.read_vcf output to be saved.
        - output_path_to_vcf_file: str output vcf path.       
    '''
    
    if output_path_to_vcf_file.split(".")[-1] not in ["vcf", "bcf"]:
        output_path_to_vcf_file += ".vcf"
    
    ## Obtain npy matrix with SNPs
    npy = vcf_data['calldata/GT']
    length, num_dogs, num_strands = npy.shape
    npy = npy.reshape(length, num_dogs*2).T
    
    ## Obtain metadata from .vcf file
    data = vcf_data
    
    print(np.unique(npy))
    
    ## Infer chromosome length and number of samples
    npy = npy.astype(str)  ## int (original)
    chmlen, _, _ = data["calldata/GT"].shape
    h, c = npy.shape
    n = h//2
    
    ## Keep sample names if appropriate
    if "samples" in list(data.keys()) and len(data["samples"]) == n:
        data_samples = data["samples"]
    else:
        data_samples = [get_name() for _ in range(n)]
          
    ## Metadata 
    df = pd.DataFrame()
    df["CHROM"]  = data["variants/CHROM"]
    df["POS"]    = data["variants/POS"]
    df["ID"]     = data["variants/ID"]
    df["REF"]    = data["variants/REF"]
    df["ALT"]    = data["variants/ALT"][:,0]  # ONLY THE FIRST SINCE WE ONLY CARE ABOUT BI-ALLELIC SNPS HERE FOR NOW
    df["QUAL"]   = data["variants/QUAL"]
    df["FILTER"] = ["PASS"]*chmlen
    df["INFO"]   = ["."]*chmlen
    df["FORMAT"] = ["GT"]*chmlen
    
    ## Genotype data for each sample
    for i in range(n):
        # Get that particular individual's maternal and paternal snps
        maternal = npy[i*2,:].astype(str) # maternal is the first
        paternal = npy[i*2+1,:].astype(str) # paternal is the second

        # Create "maternal|paternal"
        lst = [maternal, ["|"]*chmlen, paternal]
        genotype_dog = list(map(''.join, zip(*lst)))
        df[data_samples[i]] = genotype_dog

    # Write header
    with open(output_path_to_vcf_file,"w") as f:
        f.write("##fileformat=VCFv4.1\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n')
        f.write("#"+"\t".join(df.columns)+"\n") # Mandatory header
    
    # Genotype data
    df.to_csv(output_path_to_vcf_file,sep="\t",index=False,mode="a",header=False)
    
    return

def filter_by_chromosome(vcf_data, chrom):
    '''
    Objective:
        - Filter vcf_data to keep the SNPs of the specified chromosome.
    Input:
        - vcf_data: allel.read_vcf output.
        - chrom: name/number of the chromosome for which we want to select the SNPs.
    Output:
        - vcf_data: vcf_data containing only the selected SNPs of the desired chromosome.
    '''
    
    ## Find the indexes of the SNPs of the desired chromosomme
    indexes_chr = vcf_data['variants/CHROM'] == str(chrom)
    
    ## Select data for the desired SNPs
    for key in vcf_data.keys():
        if key != 'samples':
            vcf_data[key] = vcf_data[key][indexes_chr]
                
    return vcf_data


def rename_chromosome(vcf_data, before, after):
    '''
    Objective: rename variants/CHROM in vcf_data.
    Input:
        - vcf_data: allel.read_vcf output.
        - before: chromosome name before renaming.
        - after: chromosome name after renaming.
    Output:
        - vcf_data: vcf_data with renamed variants/CHROM.
    '''
    
    ## Rename variants/CHROM
    vcf_data['variants/CHROM'] = np.where(vcf_data['variants/CHROM'] == before, after, vcf_data['variants/CHROM']) 

    return vcf_data

def search_percentage_SNPs_with_missings(vcf_data):
    '''
    Objective:
        - Return the % of SNPs with some missing value. Missing values are assumed to be encoded as "-1".
    Input:
        - vcf_data: allel.read_vcf output.
    Output:
        - perc_missings: % of SNPs with some missing value in calldata/GT.
    '''
    
    ## Obtain npy matrix with SNPs
    npy = vcf_data['calldata/GT']
    length, num_dogs, num_strands = npy.shape
    npy = npy.reshape(length, num_dogs*2).T
    
    ## Obtain number of missing values encoded as "-1" in each SNP
    num_missings = np.sum(npy == -1, axis=0)
    
    ## Obtain the % of SNPs with some missing value
    perc_missings = sum(num_missings>0)/len(num_missings)*100
    
    return perc_missings


def rename_missings(vcf_data, before, after):
    '''
    Objective: rename missings in calldata/GT in vcf_data.
    Input:
        - vcf_data: allel.read_vcf output.
        - before: missing value before renaming.
        - after: missing value after renaming.
    Output:
        - vcf_data: vcf_data with renamed missings in calldata/GT.
    '''
    
    ## Rename variants/CHROM
    vcf_data['calldata/GT'] = np.where(vcf_data['calldata/GT'] == before, after, vcf_data['calldata/GT']) 

    return vcf_data


def filter_samples(vcf_data, substrings, track_name, verbose=False):
    '''
    Objective: filter the data to only keep the samples where substrings are not in sample name.
    Input:
        - vcf_data: allel.read_vcf output.
        - substrings: substrings of the samples to be removed. 
          Example: the substring "wolf" would remove a sample with name "IrishWolfhound01".
          Note: names can be in uppercase or lowercase, they will be standardized before filtering.
        - verbose: if True, removed sample names will be written in track file.
    Output:
        - vcf_data: filtered vcf_data to only keep the desired dog breeds.
    '''
    
    ## Standardize substrings and sample names to lowercase
    substrings = list(map(lambda x: x.lower(), substrings))
    sample_names = list(map(lambda x: x.lower(), vcf_data['samples']))
    
    ## Define list of booleans with the samples to be kept (as True) and to be removed (as False).
    samples_keep_remove = []
        
    if verbose:
        track('Removed samples:', track_name)
    
    ## For each sample name, see if it contains any substring and store results in samples_keep_remove to later filter the samples
    for sample_name in sample_names:
        ## True if any substring is contained in sample name. Otherwise, false
        substring_found = any(substring in sample_name for substring in substrings)
        
        ## Store result
        # Samples with true values will be kept. Otherwise, they will be removed
        samples_keep_remove.append(not substring_found)
        
        ## Write sample names that will be removed in track file
        if verbose and substring_found:
            track('-> ' + sample_name, track_name)
    
    ## Obtain total number of removed samples
    removed_samples = len(samples_keep_remove) - sum(samples_keep_remove)
    
    if verbose and removed_samples == 0:
        track('None', track_name)
    
    ## Write total number of removed samples
    track('--> {} samples removed in total'.format(removed_samples), track_name)
    
    ## Filter the samples
    vcf_data['samples'] = vcf_data['samples'][samples_keep_remove]
    vcf_data['calldata/GT'] = vcf_data['calldata/GT'][:,samples_keep_remove,:]
    
    return vcf_data


def combine_chrom_strands(chrom_data):
    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Inputs:
        - chrom_data: matrix of chromosone data (SNPs) in
              the shape (num_dogs*2, length). The maternal
              and paternal strands should be split and
              consecutively placed in the matrix

    Outputs:
        - chrom_data_combined: numpy matrix of chromosone data
              (SNPs) with combined maternal and paternal
              strands. Shape should be (num_dogs, length)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    '''

    ## Determine the dimensions of the chromosone data (SNPs)
    # length corresponds to the length of the chromosone strands (the number of SNPs).
    # num_samples corresponds to the number of dogs available for analysis times 2 because
    # there are two rows in the matrix for each dog, one for maternal and one for paternal data
    num_samples, length = chrom_data.shape

    ## Combine the information from the maternal and paternal strands by averaging
    chrom_data_combined = (chrom_data[::2] + chrom_data[1::2])/2

    ## Determine dimensions of reshaped data
    num_dogs, length_new = chrom_data_combined.shape

    ## Check that the data was reformatted properly
    assert (length == length_new), "Number of snps per sample changed. Check reshaping function."
    assert (num_samples == num_dogs*2), "Number of samples changed. Check reshaping function."

    return chrom_data_combined


def missings_in_each_SNP(vcf_data):
    '''
    Objective:
        - Print the % of SNPs that contain some missing value encoded as "-1".
    Input:
        - vcf_data: allel.read_vcf output.
    '''
    
    ## Obtain npy matrix with SNPs
    npy = vcf_data['calldata/GT']
    length, num_dogs, num_strands = npy.shape
    npy = npy.reshape(length, num_dogs*2).T
    
    ## Obtain number of samples with some missing value encoded as "-1" for each SNP
    num_neg = np.sum(npy == -1, axis=0)

    ## Print the % of SNPs that contain some missing value
    print(sum(num_neg>0)/len(num_neg)*100)