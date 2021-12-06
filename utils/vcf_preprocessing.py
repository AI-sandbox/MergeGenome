################################################################################
# Renaming chromosomes, removing ambigous SNPs, correcting SNP flips, 
# and removing mismatches between SNPs in two .vcf files
################################################################################

import numpy as np
from utils.track import track
    
def select_snps(vcf_data, indexes):
    '''
    Objective: select the SNPs at the specified indexes (in indexes param).
    Input:
        - vcf_data: allel.read_vcf output.
        - indexes: list with the indexes to keep.
    Output:
        - vcf_data: vcf_data containing the SNPs at the specified indexes.
    '''

    ## Keep all the SNPs at the specified index positions for all the keys different from samples
    for key in vcf_data.keys():
        if key != 'samples':
            vcf_data[key] = vcf_data[key][indexes]
            
    return vcf_data

def remove_snps(vcf_data, indexes):
    '''
    Objective: remove the SNPs at the specified indexes (in indexes param).
    Input:
        - vcf_data: allel.read_vcf output.
        - indexes: list with the indexes to remove.
    Output:
        - vcf_data: vcf_data containing the SNPs at the specified indexes.
    '''
    
    ## Remove all the SNPs at the specified index positions for all the keys different from samples
    for key in vcf_data.keys():
        if key != 'samples':
            vcf_data[key] = np.delete(vcf_data[key], indexes, axis=0)
            
    return vcf_data


def count_snps_same_position(vcf_data_1, vcf_data_2, dataset_name_1, dataset_name_2, track_name):
    '''
    Objective: count the number of SNPs at the same position between the two datasets.
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2.
    '''
    ## Define iterators over dataset 1 and 2
    index1 = 0
    index2 = 0
    
    ## Define counters of SNPs that coincide in POS
    n_coincidences = 0
    while index1 < len(vcf_data_1['variants/POS']) and index2 < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/POS'][index1] < vcf_data_2['variants/POS'][index2]:
            index1 += 1
        elif vcf_data_1['variants/POS'][index1] > vcf_data_2['variants/POS'][index2]:
            index2 += 1
        ## Found a SNP at the same position
        else:
            n_coincidences += 1
            index1 += 1
            index2 += 1
    
    track('--> {} SNPs found at the same position between {} and {} datasets'.format(n_coincidences, dataset_name_1, dataset_name_2), track_name)


def search_and_remove_ambiguous_snps(vcf_data, track_name='preprocessing.txt'):
    '''
    Objective:
        - Search and remove non-ambiguous SNPs from vcf data.
    Input:
        - vcf_data: allel.read_vcf output. Might contain ambiguous SNPs to be removed.
          Note: ambiguous SNPs can be A-T, T-A, C-G, or G-C pairs where the first element is the reference 
          and the second the alternate.
        - track_name: name of .txt file to write results and findings.
    Output:
        - vcf_data: contains dictionary with vcf data only for non-ambiguous SNPs (ambiguous SNPs are removed).
    '''

    ## Obtain REF and ALT from vcf_data
    REF = vcf_data['variants/REF']
    ALT = vcf_data['variants/ALT']
    
    ## Create list that will store the indexes of non-ambiguous SNPs
    non_ambiguous_idx = []
    
    ## Create counters of the number of ambiguities of type A-T, T-A, C-G, and G-C in this order
    A_T_count = 0
    T_A_count = 0
    C_G_count = 0
    G_C_count = 0

    ## Seach for ambiguous SNPs
    for i in range (0, len(REF)):
        if REF[i] == 'A' and ALT[i][0] == 'T':    # Ambiguity of A-T type
            A_T_count += 1
        elif REF[i] == 'T' and ALT[i][0] == 'A':  # Ambiguity of T-A type
            T_A_count += 1
        elif REF[i] == 'C' and ALT[i][0] == 'G':  # Ambiguity of C-G type
            C_G_count += 1
        elif REF[i] == 'G' and ALT[i][0] == 'C':  # Ambiguity of G-C type
            G_C_count += 1
        else:
            non_ambiguous_idx.append(i)

    ## Write how many SNP ambiguities of each type were found in .txt file with name track_name
    track('--> {} ambiguities found of A-T type'.format(A_T_count), track_name)
    track('--> {} ambiguities found of T-A type'.format(T_A_count), track_name)
    track('--> {} ambiguities found of C-G type'.format(C_G_count), track_name)
    track('--> {} ambiguities found of G-C type'.format(G_C_count), track_name)
    track('--> {} ambiguous SNPs in total'.format(A_T_count+T_A_count+C_G_count+G_C_count), track_name)
    
    ## Keep the SNPs that are not ambiguous 
    vcf_data = select_snps(vcf_data, non_ambiguous_idx)
            
    track('--> All SNPs with ambiguities removed correctly', track_name)
    
    return vcf_data


def correct_flips(vcf_data, indexes):
    '''
    Objective:
        - Correct flips in the dataset by swapping the reference and the alternate, and changing 0's by 1's and 1's by 0's.
    Input:
        - vcf_data: allel.read_vcf output.
        - indexes: list with the indexes of the SNPs with flips to be corrected.
    Output:
        - vcf_data: contains dictionary with vcf data with corrected flips.
    '''
    if len(indexes) > 0:
        ## Swapping the reference and the alternate
        ref_aux = vcf_data['variants/REF'][indexes]
        alt_aux = vcf_data['variants/ALT'][indexes,0]
        vcf_data['variants/REF'][indexes] = vcf_data['variants/ALT'][indexes,0]
        vcf_data['variants/ALT'][indexes,0] = ref_aux
        
        ## Ensure that the reference and the alternate were correctly swapped
        assert (vcf_data['variants/ALT'][indexes,0] == ref_aux).all(), 'The reference and the alternate were not swapped correctly'
        assert (vcf_data['variants/REF'][indexes] == alt_aux).all(), 'The reference and the alternate were not swapped correctly'
        
        ## Changing 0's by 1's and 1's by 0's
        snps = vcf_data['calldata/GT'][indexes,:,:]
        vcf_data['calldata/GT'][indexes,:,:] = np.where((snps==0)|(snps==1), snps^1, snps)
        
        ## Ensure the 0's and 1's where correctly swapped by veryfing that the number of 1's before the swap
        # Is the same as the number of 0' after the swap
        number_ones_before = np.sum(snps == 1)
        number_zeros_after = np.sum(vcf_data['calldata/GT'][indexes,:,:] == 0)
        assert number_ones_before == number_zeros_after, 'The zeros and ones where not swapped correctly'
    
    return vcf_data


def search_and_correct_flips_by_pos(vcf_data_1, vcf_data_2, track_name='preprocessing.txt'):
    '''
    Objective:
        - Search for possible flips in the reference and alternate for SNPs at the same position 
          in the two datasets (vcf_data_1 and vcf_data_2).
          Note: a SNP is flipped when the reference in the first dataset is the alternate in the second dataset, 
          and the alternate in the first dataset is the reference in the second dataset. 
        - Correct the flips by swapping the reference and the alternate of the second dataset.
          Also, the 0's are changed by 1's, and viceversa for the flippped SNPs.
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2.
        - track_name: name of .txt file to write results and findings.
    '''
    
    ## Define counters of SNPs that coincide in POS
    n_coincidences = 0
    
    ## Create list that will store the indexes of SNPs that are flipped (in dataset 2)
    flip_idx2 = []
    
    ## Define iterators over dataset 1 and 2
    index1 = 0
    index2 = 0
    while index1 < len(vcf_data_1['variants/POS']) and index2 < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/POS'][index1] < vcf_data_2['variants/POS'][index2]:
            index1 += 1
        elif vcf_data_1['variants/POS'][index1] > vcf_data_2['variants/POS'][index2]:
            index2 += 1
        ## Found a SNP at the same position
        else:   
            ## Search if there is a SNP flip
            if (vcf_data_1['variants/REF'][index1] == vcf_data_2['variants/ALT'][index2,0] and 
            vcf_data_1['variants/ALT'][index1,0] == vcf_data_2['variants/REF'][index2]):
                ## Found a SNP flip
                flip_idx2.append(index2)
            
            n_coincidences += 1
            index1 += 1
            index2 += 1
    
    track('--> {} SNPs found at the same position in the two datasets'.format(n_coincidences), track_name)
    track('--> {} flips found in total'.format(len(flip_idx2)), track_name)
    
    vcf_data_2 = correct_flips(vcf_data_2, flip_idx2)
    
    track('--> All SNPs with flips corrected correctly', track_name)
    
    return vcf_data_2


def search_and_remove_mismatches_by_pos(vcf_data_1, vcf_data_2, track_name='preprocessing.txt'):
    '''
    Objective:
        - Search and remove mismatches between vcf_data_1 and vcf_data_2.
          Note: a mismatch is a difference in the reference, the alternate, the identifier or the position of the same 
          SNP in the two datasets.
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2.
        - track_name: name of .txt file to write results and findings.
    '''
    
    ## Define counters of SNPs that coincide in POS
    n_coincidences = 0
    
    ## Create lists that will store all the indexes of the SNPs that contain a mismatch at the same position
    mismatch_idx1 = []
    mismatch_idx2 = []
    
    ## Define iterators over dataset 1 and 2
    index1 = 0
    index2 = 0
    while index1 < len(vcf_data_1['variants/POS']) and index2 < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/POS'][index1] < vcf_data_2['variants/POS'][index2]:
            index1 += 1
        elif vcf_data_1['variants/POS'][index1] > vcf_data_2['variants/POS'][index2]:
            index2 += 1
        ## Found a SNP at the same position
        else:
            ## Search if there is a mismatch
            if (vcf_data_1['variants/REF'][index1] != vcf_data_2['variants/REF'][index2] or 
                vcf_data_1['variants/ALT'][index1,0] != vcf_data_2['variants/ALT'][index2,0]):
                mismatch_idx1.append(index1)
                mismatch_idx2.append(index2)
            n_coincidences += 1
            index1 += 1
            index2 += 1
    
    track('--> {} SNPs found at the same position in the two datasets'.format(n_coincidences), track_name)
    track('--> {} mismatches found in total'.format(len(mismatch_idx1)), track_name)
  
    ## Remove the SNPs that have a mismatch
    vcf_data_1 = remove_snps(vcf_data_1, mismatch_idx1)
    vcf_data_2 = remove_snps(vcf_data_2, mismatch_idx2)
    
    track('--> All SNPs with mismatches removed correctly', track_name)
    
    return vcf_data_1, vcf_data_2
    

def search_and_keep_common_markers(vcf_data_1, vcf_data_2, track_name):
    '''
    Objective:
        - Search and keep common markers from vcf data.
          Note: common markers are SNPs with same POS, REF, and ALT.
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2.
        - track_name: name of .txt file to write results and findings.
    Output:
        - vcf_data_1:
    '''
    i = 0
    j = 0
    indexes_1 = []
    indexes_2 = []

    while i < len(vcf_data_1['variants/POS']) and j < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/POS'][i] < vcf_data_2['variants/POS'][j]:
            i += 1
        elif vcf_data_1['variants/POS'][i] > vcf_data_2['variants/POS'][j]:
            j += 1
        else: 
            if (vcf_data_1['variants/REF'][i] == vcf_data_2['variants/REF'][j] and 
                vcf_data_1['variants/ALT'][i][0] == vcf_data_2['variants/ALT'][j][0]):
                indexes_1.append(i)
                indexes_2.append(j)
            i += 1
            j += 1
    
    ## Keep the SNPs that are common markers
    vcf_data_1 = select_snps(vcf_data_1, indexes_1)
    vcf_data_2 = select_snps(vcf_data_2, indexes_2)
     
    track('{} common markers in total'.format(len(indexes_1)), track_name)
    
    return vcf_data_1, vcf_data_2