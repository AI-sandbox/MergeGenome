import allel

def read_vcf_file(path_to_vcf_file):
    '''
    Objective:
        - Read .vcf file and return the content of the .vcf file.
    Input:
        - path_to_vcf_file: path to .vcf file.
    Output:
        - vcf_data: dictionary with the content of .vcf file.
    '''
    
    ## Read the .vcf file using skikit-allel
    vcf_data = allel.read_vcf(path_to_vcf_file)
    
    return vcf_data


def write_vcf_file(vcf_data, output_path_to_vcf_file):
    
    '''
    Objective:
        - Write vcf_data in .vcf file.
    Input:
        - vcf_data: allel.read_vcf output to be saved.
        - output_path_to_vcf_file: str to .vcf path.       
    '''
    
    if output_path_to_vcf_file.split(".")[-1] not in ["vcf", "bcf"]:
        output_path_to_vcf_file += ".vcf"
    
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