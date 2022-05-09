################################################################################
# Util functions to read and write .vcf files
################################################################################

import os
import allel
import logging
import pandas as pd
from typing import Dict


def read_vcf_file(input_path: str, logger: logging.Logger) -> dict:
    """
    Reads data .vcf file.
    
    Args:
        input_path (str): path to .vcf file.
        logger (logging.Logger): debug/information tracker.
        
    Returns:
        data (Dict): dictionary with the content of .vcf file. If the input_path
        is None, then returns None.

    """
    
    if input_path is not None:
        # Read the .vcf file using skikit-allel
        data = allel.read_vcf(input_path)
        
        # Obtain dimensions of the data
        # That is, the amount of samples and of SNPs
        logger.info(f'There are {len(data["variants/ID"])} SNPs and {len(data["samples"])} '\
                    f'samples in {input_path}.')
    else:
        data = None
    
    return data


def write_vcf_file(vcf_data: Dict, output_path: str, output_name: str) -> None:
    """
    Write vcf_data in .vcf file in the specified output path.
    
    Args:
        vcf_data (Dict): dictionary with the content of .vcf file to be saved.
        output_path (str): folder to save the output.
        output_name (str): name of the file to be saved.
    
    Returns:
        (None)

    """
    
    if output_name.split(".")[-1] not in ["vcf", "bcf"]:
        output_name += ".vcf"
    
    ## Obtain npy matrix with SNPs
    npy = vcf_data['calldata/GT']
    length, num_dogs, num_strands = npy.shape
    npy = npy.reshape(length, num_dogs*2).T
    
    ## Obtain metadata from .vcf file
    data = vcf_data
        
    ## Infer chromosome length and number of samples
    npy = npy.astype(str)
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
    
    # Genotype data for each sample
    for i in range(n):
        # Get that particular individual's maternal and paternal snps
        maternal = npy[i*2,:].astype(str) # maternal is the first
        paternal = npy[i*2+1,:].astype(str) # paternal is the second

        # Create "maternal|paternal"
        lst = [maternal, ["|"]*chmlen, paternal]
        genotype_dog = list(map(''.join, zip(*lst)))
        df[data_samples[i]] = genotype_dog
    
    # If the directory does not already exist
    # create the directory that will contain the result
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    # Write header
    with open(output_path+output_name,"w") as f:
        f.write("##fileformat=VCFv4.1\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n')
        f.write("#"+"\t".join(df.columns)+"\n") # Mandatory header
    
    # Genotype data
    df.to_csv(output_path+output_name, sep="\t", index=False, mode="a", header=False)
