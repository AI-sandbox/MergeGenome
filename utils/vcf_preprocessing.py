################################################################################
# Performs small preprocessing to .vcf data: 
#   - Removing ambiguous SNPs in a .vcf file
#   - Correcting SNP flips between SNPs in two .vcf files
#   - Removing SNP mismatches between SNPs in two .vcf files
#   - etc.
################################################################################

import numpy as np
from utils.vcf_utils import combine_chrom_strands
from utils.track import track
    
def select_snps(vcf_data, indexes):
    '''
    Objective: select the SNPs at the specified index positions (in indexes param).
    Input:
        - vcf_data: allel.read_vcf output.
        - indexes: list with the indexes of the SNPs to be kept. The SNPs that are not in this list will be removed.
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
    Objective: remove the SNPs at the specified index positions (in indexes param).
    Input:
        - vcf_data: allel.read_vcf output.
        - indexes: list with the indexes to be removed. The SNPs that are not in this list will be kept.
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
    Objective: count the number of SNPs at the same position between the two datasets (vcf_data_1 and vcf_data_2).
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2.
        - dataset_name_1: name name given to datasets 1.
        - dataset_name_2: name given to datasets 2.
    '''
    ## Define iterators over dataset 1 and 2
    index1 = 0
    index2 = 0
    
    ## Define counters of SNPs that are at the same position (same POS)
    n_coincidences = 0
    
    while index1 < len(vcf_data_1['variants/POS']) and index2 < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/POS'][index1] < vcf_data_2['variants/POS'][index2]:
            ## If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            ## Increase the iterator of dataset 1
            index1 += 1
        elif vcf_data_1['variants/POS'][index1] > vcf_data_2['variants/POS'][index2]:
            ## If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            ## Increase the iterator of dataset 2
            index2 += 1
        else:
            ## Found a SNP at the same position
            ## Increase the counter of coincidences and the iterator of datasets 1 and 2
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
          Note: ambiguous SNPs can be A/T, T/A, C/G, or G/C pairs where the first element is the reference (REF) and the second the alternate (ALT).
        - track_name: name of .txt file to write results.
    Output:
        - vcf_data: allel.read_vcf output only with non-ambiguous SNPs (ambiguous SNPs are removed).
    '''

    ## Obtain REF and ALT from vcf_data
    REF = vcf_data['variants/REF']
    ALT = vcf_data['variants/ALT']
    
    ## Create list that will store the index positions of non-ambiguous SNPs
    non_ambiguous_idx = []
    
    ## Create counters of the number of ambiguities of type A/T, T/A, C/G, and G/C in this order
    A_T_count = 0
    T_A_count = 0
    C_G_count = 0
    G_C_count = 0

    ## Seach for ambiguous SNPs
    for i in range (0, len(REF)):
        if REF[i] == 'A' and ALT[i][0] == 'T':
            ## If the reference is an A and the alternate is a T...
            ## Found ambiguity of A/T type
            A_T_count += 1
        elif REF[i] == 'T' and ALT[i][0] == 'A':
            ## If the reference is an T and the alternate is a A...
            ## Found ambiguity of T/A type
            T_A_count += 1
        elif REF[i] == 'C' and ALT[i][0] == 'G':
            ## If the reference is an C and the alternate is a G...
            ## Found ambiguity of C/G type
            C_G_count += 1
        elif REF[i] == 'G' and ALT[i][0] == 'C':
            ## If the reference is an G and the alternate is a C...
            ## Found ambiguity of G/C type
            G_C_count += 1
        else:
            # The SNP at position i is not ambiguous
            non_ambiguous_idx.append(i)

    ## Write how many ambiguities of each type were found in .txt file with name track_name
    track('--> {} ambiguities found of A-T type'.format(A_T_count), track_name)
    track('--> {} ambiguities found of T-A type'.format(T_A_count), track_name)
    track('--> {} ambiguities found of C-G type'.format(C_G_count), track_name)
    track('--> {} ambiguities found of G-C type'.format(G_C_count), track_name)
    track('--> {} ambiguous SNPs in total'.format(A_T_count+T_A_count+C_G_count+G_C_count), track_name)
    
    ## Keep the SNPs that are not ambiguous
    # Ambigous SNPs are removed next
    vcf_data = select_snps(vcf_data, non_ambiguous_idx)
            
    track('--> All SNPs with ambiguities removed correctly', track_name)
    
    return vcf_data


def correct_flips(vcf_data, indexes):
    '''
    Objective:
        - Correct flips in the dataset by swapping the reference (REF) and the alternate (ALT), and also changing 0's by 1's and 1's by 0's.
          The index positions of the SNPs with flips to be corrected are found in indexes.
    Input:
        - vcf_data: allel.read_vcf output.
        - indexes: list with the index positions of the SNPs that are flipped.
    Output:
        - vcf_data: allel.read_vcf output with corrected SNP flips.
    '''
    
    if len(indexes) > 0:
        ## If there is some SNP flip...
        ## Swap the reference and the alternate
        ref_aux = vcf_data['variants/REF'][indexes]
        alt_aux = vcf_data['variants/ALT'][indexes,0]
        vcf_data['variants/REF'][indexes] = vcf_data['variants/ALT'][indexes,0]
        vcf_data['variants/ALT'][indexes,0] = ref_aux
        
        ## Ensure that the reference and the alternate were correctly swapped
        assert (vcf_data['variants/ALT'][indexes,0] == ref_aux).all(), 'The reference and the alternate were not swapped correctly'
        assert (vcf_data['variants/REF'][indexes] == alt_aux).all(), 'The reference and the alternate were not swapped correctly'
        
        ## Change 0's by 1's and 1's by 0's
        snps = vcf_data['calldata/GT'][indexes,:,:]
        vcf_data['calldata/GT'][indexes,:,:] = np.where((snps==0)|(snps==1), snps^1, snps)
        
        ## Ensure the 0's and 1's where correctly swapped, by checking that the number of 1's before the swap
        # is the same as the number of 0' after the swap
        number_ones_before = np.sum(snps == 1)
        number_zeros_after = np.sum(vcf_data['calldata/GT'][indexes,:,:] == 0)
        assert number_ones_before == number_zeros_after, 'The zeros and ones where not swapped correctly'
    
    return vcf_data


def search_and_correct_flips_by_pos(vcf_data_1, vcf_data_2, track_name='preprocessing.txt'):
    '''
    Objective:
        - Search for possible SNP flips in the reference and alternate, for the SNPs that are at the same position
          between the two datasets (vcf_data_1 and vcf_data_2).
          Note: a SNP is flipped when the reference (REF) in the first dataset (vcf_data_1) is the alternate (ALT) in the second dataset (vcf_data_2), 
          and the alternate in the former is the reference in the latter. 
        - Correct SNPs flips by swapping the reference and the alternate of the second dataset.
          Also, the 0's are changed by 1's and viceversa for the those SNPs with a flip.
    Input:
        - vcf_data_2: allel.read_vcf output with corrected SNP flips.
        - track_name: name of .txt file to write results.
    '''
    
    ## Define counters of SNPs that are at the same position (same POS)
    n_coincidences = 0
    
    ## Create list that will store the index positions of the SNPs with a flip (in dataset 2)
    flip_idx2 = []
    
    ## Define iterators over the SNPs in dataset 1 and 2
    index1 = 0
    index2 = 0
    while index1 < len(vcf_data_1['variants/POS']) and index2 < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/POS'][index1] < vcf_data_2['variants/POS'][index2]:
            ## If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            ## Increase the iterator of dataset 1
            index1 += 1
        elif vcf_data_1['variants/POS'][index1] > vcf_data_2['variants/POS'][index2]:
            ## If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            ## Increase the iterator of dataset 2
            index2 += 1
        else:   
            ## Found a SNP at the same position...
            ## Search if there is a SNP flip at this position
            if (vcf_data_1['variants/REF'][index1] == vcf_data_2['variants/ALT'][index2,0] and 
            vcf_data_1['variants/ALT'][index1,0] == vcf_data_2['variants/REF'][index2]):
                ## Found a SNP flip. The reference in dataset 1 is the alternate in dataset 2, and the alternate in dataset 1 is the reference in dataset 2...
                ## Therefore, append the index position of the SNP with a flip in dataset 2 to flip_idx2
                flip_idx2.append(index2)
            
            ## Increase the counter of coincidences and the iterator of datasets 1 and 2
            n_coincidences += 1
            index1 += 1
            index2 += 1
    
    track('--> {} SNPs found at the same position in the two datasets'.format(n_coincidences), track_name)
    track('--> {} flips found in total'.format(len(flip_idx2)), track_name)
    
    ## Correct SNPs flips at the indexes found in flip_idx2
    vcf_data_2 = correct_flips(vcf_data_2, flip_idx2)
    
    track('--> All SNPs with flips corrected correctly', track_name)
    
    return vcf_data_2


def search_and_remove_mismatches_by_pos(vcf_data_1, vcf_data_2, track_name='preprocessing.txt'):
    '''
    Objective:
        - Search and remove mismatches between vcf_data_1 and vcf_data_2.
          Note: there is a mismatch when there is a difference in the reference (REF) or alternate (ALT) between SNPs at the same position 
          in two datasets (vcf_data_1 and vcf_data_2). Note: it is recommended to first correct the SNP flips prior to removing SNP mismatches.
          You can correct the SNP flips with search_and_correct_flips_by_pos.
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1 without SNPs with mismatches.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2 without removed SNPs with mismatches.
        - track_name: name of .txt file to write results.
    '''
    
    ## Define counters of SNPs that are at the same position (same POS)
    n_coincidences = 0
    
    ## Create lists that will store all the indexes of the SNPs that contain a mismatch at the same position
    mismatch_idx1 = []
    mismatch_idx2 = []
    
    ## Define iterators over dataset 1 and 2
    index1 = 0
    index2 = 0
    while index1 < len(vcf_data_1['variants/POS']) and index2 < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/POS'][index1] < vcf_data_2['variants/POS'][index2]:
            ## If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            ## Increase the iterator of dataset 1
            index1 += 1
        elif vcf_data_1['variants/POS'][index1] > vcf_data_2['variants/POS'][index2]:
            ## If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            ## Increase the iterator of dataset 2
            index2 += 1
        else:
            ## Found a SNP at the same position...
            ## Search if there is a mismatch
            if (vcf_data_1['variants/REF'][index1] != vcf_data_2['variants/REF'][index2] or 
                vcf_data_1['variants/ALT'][index1,0] != vcf_data_2['variants/ALT'][index2,0]):
                ## If the reference or the alternate between datasets 1 and 2 do not match...
                # Save the index positions of the SNPs that present a mismatch in each dataset
                mismatch_idx1.append(index1)
                mismatch_idx2.append(index2)
            
            ## Increase the counter of coincidences and the iterator of datasets 1 and 2
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
    

def search_and_keep_common_markers_single_chr(vcf_data_1, vcf_data_2, track_name):
    '''
    Objective:
        - Search and keep common markers between two datasets (vcf_data_1 and vcf_data_2).
          Note: common markers are SNPs with same chromosome (CHROM), position (POS), reference (REF), and alernate (ALT). This function assumes
          vcf_data_1 and vcf_data_2 contain data for a single chromosome (the same in both datasets).
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2.
        - track_name: name of .txt file to write results.
    Output:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1 for the common markers.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2 for the common markers.
        - indexes_1: indexes of the SNPs in the .vcf file 1 that are common markers with the .vcf file 2.
        - indexes_2: indexes of the SNPs in the .vcf file 2 that are common markers with the .vcf file 1.
    '''
    
    ## Define iterators over dataset 1 and 2
    i = 0
    j = 0
    
    ## Create lists that will store all the indexes of the SNPs that are common markers
    indexes_1 = []
    indexes_2 = []

    while i < len(vcf_data_1['variants/POS']) and j < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/POS'][i] < vcf_data_2['variants/POS'][j]:
            ## If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            ## Increase the iterator of dataset 1
            i += 1
        elif vcf_data_1['variants/POS'][i] > vcf_data_2['variants/POS'][j]:
            ## If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            ## Increase the iterator of dataset 2
            j += 1
        else:
            ## Found a SNP at the same position...
            ## Search if it corresponds to a common marker
            if (vcf_data_1['variants/REF'][i] == vcf_data_2['variants/REF'][j] and 
                vcf_data_1['variants/ALT'][i][0] == vcf_data_2['variants/ALT'][j][0]):
                ## If the reference and the alternate between datasets 1 and 2 match...
                # Save the index positions of the SNPs that are a common marker in each dataset
                indexes_1.append(i)
                indexes_2.append(j)
            
            ## Increase the iterator of datasets 1 and 2
            i += 1
            j += 1
    
    ## Keep the SNPs that are common markers
    # The SNPs that differe in position, reference or alternate are removed
    vcf_data_1 = select_snps(vcf_data_1, indexes_1)
    vcf_data_2 = select_snps(vcf_data_2, indexes_2)
     
    track('{} common markers in total'.format(len(indexes_1)), track_name)
    
    return vcf_data_1, vcf_data_2, indexes_1, indexes_2


def search_and_keep_common_markers_several_chr(vcf_data_1, vcf_data_2, track_name):
    '''
    Objective:
        - Search and keep common markers between two datasets (vcf_data_1 and vcf_data_2).
          Note: common markers are SNPs with same chromosome (CHROM), position (POS), reference (REF), and alernate (ALT). This function accepts
          vcf_data_1 and vcf_data_2 data for multiple chromosomes.
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2.
        - track_name: name of .txt file to write results.
    Output:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1 for the common markers.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2 for the common markers.
        - indexes_1: indexes of the SNPs in the .vcf file 1 that are common markers with the .vcf file 2.
        - indexes_2: indexes of the SNPs in the .vcf file 2 that are common markers with the .vcf file 1.
    '''
    
    ## Define iterators over dataset 1 and 2
    i = 0
    j = 0

    ## Create lists that will store all the indexes of the SNPs that are common markers
    indexes_1 = []
    indexes_2 = []

    while i < len(vcf_data_1['variants/POS']) and j < len(vcf_data_2['variants/POS']):
        if vcf_data_1['variants/CHROM'][i][3:] < vcf_data_2['variants/CHROM'][j][3:]:
            i += 1
        elif vcf_data_1['variants/CHROM'][i][3:] > vcf_data_2['variants/CHROM'][j][3:]:
            j += 1
        elif vcf_data_1['variants/POS'][i] < vcf_data_2['variants/POS'][j]:
            ## If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            ## Increase the iterator of dataset 1
            i += 1
        elif vcf_data_1['variants/POS'][i] > vcf_data_2['variants/POS'][j]:
            ## If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            ## Increase the iterator of dataset 2
            j += 1
        else:
            ## Found a SNP at the same position...
            ## Search if it corresponds to a common marker
            if (vcf_data_1['variants/REF'][i] == vcf_data_2['variants/REF'][j] and 
                vcf_data_1['variants/ALT'][i][0] == vcf_data_2['variants/ALT'][j][0]):
                ## If the reference and the alternate between datasets 1 and 2 match...
                # Save the index positions of the SNPs that are a common marker in each dataset
                indexes_1.append(i)
                indexes_2.append(j)

                ## Increase the iterator of datasets 1 and 2
                i += 1
                j += 1

    track('{} common markers in total'.format(len(indexes_1)), track_name)
    
    return vcf_data_1, vcf_data_2, indexes_1, indexes_2


def search_and_remove_snps_different_means(vcf_data_1, vcf_data_2, threshold, track_name):
    '''
    Objective:
        - Search and remove common SNPs between vcf_data_1 and vcf_data_2 with a mean absolute difference > threshold.
          Note: common SNPs are SNPs with same chromosome (CHROM), position (POS), reference (REF), and alernate (ALT). This function assumes
          vcf_data_1 and vcf_data_2 contain data for a single chromosome (the same in both datasets).
    Input:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2.
        - threshold: value between 0.0 and 1.0. All SNPs with a mean absolute difference higher than this value will be removed.
        - track_name: name of .txt file to write results.
    Output:
        - vcf_data_1: allel.read_vcf output of the .vcf file 1 without the SNPs with a mean absolute difference > threshold.
        - vcf_data_2: allel.read_vcf output of the .vcf file 2 without the SNPs with a mean absolute difference > threshold.
    '''
        
    ## Obtain SNPs data for datasets 1 and 2
    snps_1 = vcf_data_1['calldata/GT']
    snps_2 = vcf_data_2['calldata/GT']
        
    ## Convert SNPs from MxNx2 to MxN shape (averaged)
    # For dataset 1...
    length_1, num_dogs_1, num_strands_1 = snps_1.shape
    snps_1 = snps_1.reshape(length_1, num_dogs_1*2).T
    snps_1 = combine_chrom_strands(snps_1)
    # And for dataset 2...
    length_2, num_dogs_2, num_strands_2 = snps_2.shape
    snps_2 = snps_2.reshape(length_2, num_dogs_2*2).T
    snps_2 = combine_chrom_strands(snps_2)
    
    ## Search indexes of common markers
    _, _, indexes_1, indexes_2 = search_and_keep_common_markers_single_chr(vcf_data_1.copy(), vcf_data_2.copy(), track_name)
    
    ## Select SNPs data of common markers
    snps_1 = snps_1[:, indexes_1]
    snps_2 = snps_2[:, indexes_2]
    
    ## Compute the mean of the common SNPs in each dataset
    mean_1 = np.mean(snps_1, axis=0)
    mean_2 = np.mean(snps_2, axis=0)
    
    ## Compute the indexes of the SNPs to be removed in datasets 1 and 2
    idx_1_to_remove = np.array(indexes_1)[abs(mean_1 - mean_2) > threshold]
    idx_2_to_remove = np.array(indexes_2)[abs(mean_1 - mean_2) > threshold]
    
    assert len(idx_1_to_remove) == len(idx_2_to_remove) == sum(abs(mean_1 - mean_2) > threshold), 'The amount of indexes removed in each dataset is not the same. There has been an error.'
    
    track('Found {} SNPs with a mean absolute difference higher than {}'.format(sum(abs(mean_1 - mean_2) > threshold), threshold), track_name)
        
    ## Remove the SNPs that have a mean absolute difference > threshold
    vcf_data_1 = remove_snps(vcf_data_1, idx_1_to_remove)
    vcf_data_2 = remove_snps(vcf_data_2, idx_2_to_remove)
    
    track('All SNPs with a mean absolute difference higher than {} removed correctly'.format(threshold), track_name)
    
    return vcf_data_1, vcf_data_2
