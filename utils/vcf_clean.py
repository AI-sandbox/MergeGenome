################################################################################
# Performs small preprocessing to .vcf data: 
# * Removes ambiguous SNPs in a .vcf file.
# * Corrects SNP flips between SNPs in two .vcf files
# * Remove SNP mismatches between SNPs in two .vcf files.
# * etc.
################################################################################

import numpy as np
import logging
from typing import List
from utils.vcf_utils import combine_chrom_strands


def select_snps(vcf_data: dict, indexes: List) -> dict:
    """
    Selects the SNPs at the specified index positions.
    
    Args:
        vcf_data (dict): allel.read_vcf output.
        indexes (List): list with the indexes of the SNPs to be kept. 
                        The SNPs that are not in this list will be removed.
    
    Returns:
        vcf_data (dict): vcf_data containing the SNPs at the specified indexes.
    
    """

    # Keep all the SNPs at the specified index positions for all the keys different from samples
    for key in vcf_data.keys():
        if key != 'samples':
            vcf_data[key] = vcf_data[key][indexes]

    return vcf_data


def remove_snps(vcf_data: dict, indexes: List) -> dict:
    """
    Removes the SNPs at the specified index positions (in indexes param).
    
    Args:
        vcf_data (dict): allel.read_vcf output.
        indexes (List): list with the indexes to be removed. The SNPs that are not in this list will be kept.
    
    Returns:
        vcf_data (dict): vcf_data containing the SNPs at the specified indexes.
    
    """
    
    # Remove all the SNPs at the specified index positions for all the keys different from samples
    for key in vcf_data.keys():
        if key != 'samples':
            vcf_data[key] = np.delete(vcf_data[key], indexes, axis=0)
            
    return vcf_data


def remove_ambiguous_snps(vcf_data: dict, logger: logging.Logger) -> dict:
    """
    Searches and removes non-ambiguous SNPs from vcf data.
    
    Args:
        vcf_data (dict): allel.read_vcf output. Might contain ambiguous SNPs to be removed.
                         Note: ambiguous SNPs can be A/T, T/A, C/G, or G/C pairs where the 
                         first element is the reference (REF) and the second the alternate (ALT).
        logger (logging.Logger): debug/information tracker.
    
    Returns:
    
        vcf_data (dict): allel.read_vcf output only with non-ambiguous SNPs (ambiguous SNPs are removed).
    
    """

    # Obtain REF and ALT from vcf_data
    REF = vcf_data['variants/REF']
    ALT = vcf_data['variants/ALT']
    
    # Create list that will store the index positions of non-ambiguous SNPs
    non_ambiguous_idx = []
    
    # Create counters of the number of ambiguities of type A/T, T/A, C/G, and G/C in this order
    A_T_count = 0
    T_A_count = 0
    C_G_count = 0
    G_C_count = 0

    # Seach for ambiguous SNPs
    for i in range (0, len(REF)):
        if REF[i] == 'A' and ALT[i][0] == 'T':
            # If the reference is an A and the alternate is a T...
            # Found ambiguity of A/T type
            A_T_count += 1
        elif REF[i] == 'T' and ALT[i][0] == 'A':
            # If the reference is an T and the alternate is a A...
            # Found ambiguity of T/A type
            T_A_count += 1
        elif REF[i] == 'C' and ALT[i][0] == 'G':
            # If the reference is an C and the alternate is a G...
            # Found ambiguity of C/G type
            C_G_count += 1
        elif REF[i] == 'G' and ALT[i][0] == 'C':
            # If the reference is an G and the alternate is a C...
            # Found ambiguity of G/C type
            G_C_count += 1
        else:
            # The SNP at position i is not ambiguous
            non_ambiguous_idx.append(i)
            
    # Write how many ambiguities of each type were found in .txt file with name track_name
    logger.info(f'--> {A_T_count} ambiguities found of A-T type.')
    logger.info(f'--> {T_A_count} ambiguities found of T-A type.')
    logger.info(f'--> {C_G_count} ambiguities found of C-G type.')
    logger.info(f'--> {G_C_count} ambiguities found of G-C type.')
    logger.info(f'--> {A_T_count+T_A_count+C_G_count+G_C_count} ambiguous SNPs in total.')
    
    # Keep the SNPs that are not ambiguous, mbigous SNPs are removed next
    logger.debug('Removing all amgiguities.')
    vcf_data = select_snps(vcf_data, non_ambiguous_idx)
    
    return vcf_data


def correct_flips(vcf_data: dict, indexes: List) -> dict:
    """
    Corrects flips in the dataset by swapping the reference (REF) and the alternate (ALT), 
    and also changing 0's by 1's and 1's by 0's. The index positions of the SNPs with 
    flips to be corrected are found in indexes.
    
    Args:
        vcf_data (dict): allel.read_vcf output.
        indexes (List): list with the index positions of the SNPs that are flipped.
    
    Returns:
        vcf_data (dict): allel.read_vcf output with corrected SNP flips.
    
    """
    
    if len(indexes) > 0:
        # If there is some SNP flip...
        # Swap the reference and the alternate
        ref_aux = vcf_data['variants/REF'][indexes]
        alt_aux = vcf_data['variants/ALT'][indexes,0]
        vcf_data['variants/REF'][indexes] = vcf_data['variants/ALT'][indexes,0]
        vcf_data['variants/ALT'][indexes,0] = ref_aux
        
        # Ensure that the reference and the alternate were correctly swapped
        assert (vcf_data['variants/ALT'][indexes,0] == ref_aux).all(), 'The reference and the alternate were not swapped correctly'
        assert (vcf_data['variants/REF'][indexes] == alt_aux).all(), 'The reference and the alternate were not swapped correctly'
        
        # Change 0's by 1's and 1's by 0's
        snps = vcf_data['calldata/GT'][indexes,:,:]
        vcf_data['calldata/GT'][indexes,:,:] = np.where((snps==0)|(snps==1), snps^1, snps)
        
        # Ensure the 0's and 1's where correctly swapped, by checking that the number of 1's before the swap
        # is the same as the number of 0' after the swap
        number_ones_before = np.sum(snps == 1)
        number_zeros_after = np.sum(vcf_data['calldata/GT'][indexes,:,:] == 0)
        assert number_ones_before == number_zeros_after, 'The zeros and ones where not swapped correctly'
    
    return vcf_data


def correct_flips_by_pos(reference: dict, query: dict, logger: logging.Logger) -> dict:
    """
    Searches for possible SNP flips in the reference and alternate, for the SNPs that 
    are at the same position between the two datasets (reference and query). 
    Note: a SNP is flipped when the reference (REF) in the first dataset (reference)
    is the alternate (ALT) in the second dataset (query), and the alternate in the 
    former is the reference in the latter. 
    
    Corrects SNPs flips by swapping the reference and the alternate of the second dataset.
    Also, the 0's are changed by 1's and viceversa for the those SNPs with a flip.
        
    Args:
        query (dict): query data.
        reference (dict): reference data.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        query (dict): query with corrected SNP flips.
    
    """
    
    # Define counters of SNPs that are at the same position (same POS)
    n_coincidences = 0
    
    # Create list that will store the index positions of the SNPs with a flip (in dataset 2)
    flip_idx2 = []
    
    # Define iterators over the SNPs in dataset 1 and 2
    index1 = 0
    index2 = 0
    
    while index1 < len(reference['variants/POS']) and index2 < len(query['variants/POS']):
        if reference['variants/POS'][index1] < query['variants/POS'][index2]:
            # If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            # Increase the iterator of dataset 1
            index1 += 1
        elif reference['variants/POS'][index1] > query['variants/POS'][index2]:
            # If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            # Increase the iterator of dataset 2
            index2 += 1
        else:   
            # Found a SNP at the same position...
            # Search if there is a SNP flip at this position
            if (reference['variants/REF'][index1] == query['variants/ALT'][index2,0] and 
            reference['variants/ALT'][index1,0] == query['variants/REF'][index2]):
                # Found a SNP flip. The reference in dataset 1 is the alternate in dataset 2, and the alternate in dataset 1 is the reference in dataset 2...
                # Therefore, append the index position of the SNP with a flip in dataset 2 to flip_idx2
                flip_idx2.append(index2)
            
            # Increase the counter of coincidences and the iterator of datasets 1 and 2
            n_coincidences += 1
            index1 += 1
            index2 += 1
    
    logger.info(f'--> {n_coincidences} SNPs found at the same position in the two datasets')
    logger.info(f'--> {len(flip_idx2)} flips found in total')
    
    # Correct SNPs flips at the indexes found in flip_idx2
    logger.debug('Correcting all SNP flips.')
    query = correct_flips(query, flip_idx2)
    
    return query


def remove_mismatches_by_pos(reference: dict, query: dict, logger: logging.Logger) -> (dict, dict):
    """
    Searches and remove mismatches between reference and query. Note: there is a mismatch 
    when there is a difference in the reference (REF) or alternate (ALT) between SNPs at 
    the same position in two datasets (reference and query). Note: it is recommended to
    first correct the SNP flips prior to removing SNP mismatches. You can correct the SNP 
    flips with correct_flips_by_pos.
    
    Args:
        reference (dict): reference without SNPs with mismatches.
        query (dict): query without removed SNPs with mismatches.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        (dict, dict)
    
    """
    
    # Define counters of SNPs that are at the same position (same POS)
    n_coincidences = 0
    
    # Create lists that will store all the indexes of the SNPs that contain a mismatch at the same position
    mismatch_idx1 = []
    mismatch_idx2 = []
    
    # Define iterators over dataset 1 and 2
    index1 = 0
    index2 = 0
    while index1 < len(reference['variants/POS']) and index2 < len(query['variants/POS']):
        if reference['variants/POS'][index1] < query['variants/POS'][index2]:
            # If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            # Increase the iterator of dataset 1
            index1 += 1
        elif reference['variants/POS'][index1] > query['variants/POS'][index2]:
            # If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            # Increase the iterator of dataset 2
            index2 += 1
        else:
            # Found a SNP at the same position...
            # Search if there is a mismatch
            if (reference['variants/REF'][index1] != query['variants/REF'][index2] or 
                reference['variants/ALT'][index1,0] != query['variants/ALT'][index2,0]):
                # If the reference or the alternate between datasets 1 and 2 do not match...
                # Save the index positions of the SNPs that present a mismatch in each dataset
                mismatch_idx1.append(index1)
                mismatch_idx2.append(index2)
            
            # Increase the counter of coincidences and the iterator of datasets 1 and 2
            n_coincidences += 1
            index1 += 1
            index2 += 1
    
    logger.info(f'--> {n_coincidences} SNPs found at the same position in the two datasets')
    logger.info(f'--> {len(mismatch_idx1)} mismatches found in total')
  
    # Remove the SNPs that have a mismatch
    logger.debug('Correcting all SNP mismatches.')
    reference = remove_snps(reference, mismatch_idx1)
    query = remove_snps(query, mismatch_idx2)
    
    logger.info(f'There are {len(reference["variants/ID"])} SNPs and {len(reference["samples"])} samples ' \
                        f'in total in the reference after removing mismatching SNPs.')
    logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples ' \
                f'in total in the query after removing mismatching SNPs.')
    
    return reference, query
    

def keep_common_markers_single_chr(reference: dict, query: dict, logger: logging.Logger) -> (dict, dict, List, List):
    """
    Searches and keeps common markers between two datasets (reference and query).
    Note: common markers are SNPs with same chromosome (CHROM), position (POS), reference (REF), 
    and alernate (ALT). This function assumes reference and query contain data for the same chromosome.
    
    Args:
        reference (dict): dictionary with the content of .vcf file 1.
        query (dict): dictionary with the content of .vcf file 2.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        reference (dict): reference for the common markers.
        query (dict): query for the common markers.
        idxs_reference (List): indexes of the SNPs in the .vcf file 1 that are common markers with the .vcf file 2.
        idxs_query (List): indexes of the SNPs in the .vcf file 2 that are common markers with the .vcf file 1.
    
    """
    
    # Define iterators over dataset 1 and 2
    i = 0
    j = 0
    
    # Create lists that will store all the indexes of the SNPs that are common markers
    idxs_reference = []
    idxs_query = []

    while i < len(reference['variants/POS']) and j < len(query['variants/POS']):
        if reference['variants/POS'][i] < query['variants/POS'][j]:
            # If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            # Increase the iterator of dataset 1
            i += 1
        elif reference['variants/POS'][i] > query['variants/POS'][j]:
            # If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            # Increase the iterator of dataset 2
            j += 1
        else:
            # Found a SNP at the same position...
            # Search if it corresponds to a common marker
            if (reference['variants/REF'][i] == query['variants/REF'][j] and 
                reference['variants/ALT'][i][0] == query['variants/ALT'][j][0]):
                # If the reference and the alternate between datasets 1 and 2 match...
                # Save the index positions of the SNPs that are a common marker in each dataset
                idxs_reference.append(i)
                idxs_query.append(j)
            
            # Increase the iterator of datasets 1 and 2
            i += 1
            j += 1
    
    # Keep the SNPs that are common markers
    # The SNPs that differe in position, reference or alternate are removed
    reference = select_snps(reference, idxs_reference)
    query = select_snps(query, idxs_query)
     
    logger.info(f'--> {len(idxs_reference)} common markers in total')
    
    return reference, query, idxs_reference, idxs_query


def keep_common_markers_several_chr(reference: dict, query: dict, logger: logging.Logger):
    """
    Searches and keeps common markers between two datasets (reference and query).
    Note: common markers are SNPs with same chromosome (CHROM), position (POS), reference (REF), 
    and alernate (ALT). This function accepts reference and query data for multiple chromosomes.
    
    Args:
        reference (dict): dictionary with the content of .vcf file 1.
        query (dict): dictionary with the content of .vcf file 2.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        reference (dict): reference for the common markers.
        query (dict): query for the common markers.
        idxs_reference (List): indexes of the SNPs in the .vcf file 1 that are common markers with the .vcf file 2.
        idxs_query (List): indexes of the SNPs in the .vcf file 2 that are common markers with the .vcf file 1.
    
    """
    
    # Create lists that will store all the indexes of the SNPs that are common markers
    idxs_reference = []
    idxs_query = []

    # Define iterators over dataset 1 and 2
    i = 0
    j = 0

    # Create lists that will store all the indexes of the SNPs that are common markers
    idxs_reference = []
    idxs_query = []

    while i < len(reference['variants/POS']) and j < len(query['variants/POS']):
        if int(reference['variants/CHROM'][i][3:]) < int(query['variants/CHROM'][j][3:]):
            i += 1
        elif int(reference['variants/CHROM'][i][3:]) > int(query['variants/CHROM'][j][3:]):
            j += 1
        elif reference['variants/POS'][i] < query['variants/POS'][j]:
            # If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            # Increase the iterator of dataset 1
            i += 1
        elif reference['variants/POS'][i] > query['variants/POS'][j]:
            # If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            # Increase the iterator of dataset 2
            j += 1
        else:
            # Found a SNP at the same position...
            # Search if it corresponds to a common marker
            if (reference['variants/REF'][i] == query['variants/REF'][j] and 
                reference['variants/ALT'][i][0] == query['variants/ALT'][j][0]):
                # If the reference and the alternate between datasets 1 and 2 match...
                # Save the index positions of the SNPs that are a common marker in each dataset
                idxs_reference.append(i)
                idxs_query.append(j)

                # Increase the iterator of datasets 1 and 2
                i += 1
                j += 1
                
    logger.info(f'--> {len(idxs_reference)} common markers in total')
    
    # Keep the SNPs that are common markers
    # The SNPs that differe in position, reference or alternate are removed
    reference = select_snps(reference, idxs_reference)
    query = select_snps(query, idxs_query)
    
    return reference, query, idxs_reference, idxs_query


def remove_snps_different_means(reference: dict, query: dict, threshold: float, logger: logging.Logger) -> (dict, dict):
    """
    Searches and removes common SNPs between reference and query with a mean 
    absolute difference > threshold. Note: common SNPs are SNPs with same chromosome (CHROM), 
    position (POS), reference (REF), and alernate (ALT). This function assumes reference and 
    query contain data for a single chromosome (the same in both datasets).
    
    Args:
        reference (dict): dictionary with the content of .vcf file 1.
        query (dict): dictionary with the content of .vcf file 2.
        threshold (float): value between 0.0 and 1.0. All SNPs with a mean absolute difference 
                           higher than this value will be removed.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        reference (dict): dictionary with the content of .vcf file 1 without the SNPs with a mean absolute difference > threshold.
        query (dict): dictionary with the content of .vcf file 2 without the SNPs with a mean absolute difference > threshold.
    
    """

    # Obtain SNPs data for datasets 1 and 2
    snps_1 = reference['calldata/GT']
    snps_2 = query['calldata/GT']
        
    # Convert SNPs from MxNx2 to MxN shape (averaged)
    # For dataset 1...
    length_1, num_dogs_1, num_strands_1 = snps_1.shape
    snps_1 = snps_1.reshape(length_1, num_dogs_1*2).T
    snps_1 = combine_chrom_strands(snps_1)
    # And for dataset 2...
    length_2, num_dogs_2, num_strands_2 = snps_2.shape
    snps_2 = snps_2.reshape(length_2, num_dogs_2*2).T
    snps_2 = combine_chrom_strands(snps_2)
    
    # Search indexes of common markers
    _, _, idxs_reference, idxs_query = keep_common_markers_single_chr(reference.copy(), query.copy(), track_name)
    
    # Select SNPs data of common markers
    snps_1 = snps_1[:, idxs_reference]
    snps_2 = snps_2[:, idxs_query]
    
    # Compute the mean of the common SNPs in each dataset
    mean_1 = np.mean(snps_1, axis=0)
    mean_2 = np.mean(snps_2, axis=0)
    
    # Compute the indexes of the SNPs to be removed in datasets 1 and 2
    idx_1_to_remove = np.array(idxs_reference)[abs(mean_1 - mean_2) > threshold]
    idx_2_to_remove = np.array(idxs_query)[abs(mean_1 - mean_2) > threshold]
    
    assert len(idx_1_to_remove) == len(idx_2_to_remove) == sum(abs(mean_1 - mean_2) > threshold), 'The amount of indexes removed in each dataset is not the same. There has been an error.'
    
    track('Found {} SNPs with a mean absolute difference higher than {}'.format(sum(abs(mean_1 - mean_2) > threshold), threshold), track_name)
        
    # Remove the SNPs that have a mean absolute difference > threshold
    reference = remove_snps(reference, idx_1_to_remove)
    query = remove_snps(query, idx_2_to_remove)
    
    track('All SNPs with a mean absolute difference higher than {} removed correctly'.format(threshold), track_name)
    
    return reference, query


def count_snps_same_position(reference, query, dataset_name_1, dataset_name_2, logger: logging.Logger):
    """
    Objective: count the number of SNPs at the same position between the two datasets (reference and query).
    
    Args:
        reference (dict): dictionary with the content of .vcf file 1.
        query (dict): dictionary with the content of .vcf file 2.
        dataset_name_1: name name given to datasets 1.
        dataset_name_2: name given to datasets 2.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        (None)
    
    """
    
    # Define iterators over dataset 1 and 2
    index1 = 0
    index2 = 0
    
    # Define counters of SNPs that are at the same position (same POS)
    n_coincidences = 0
    
    while index1 < len(reference['variants/POS']) and index2 < len(query['variants/POS']):
        if reference['variants/POS'][index1] < query['variants/POS'][index2]:
            # If the position of the SNP in dataset 1 is smaller than the position of the SNP in dataset 2...
            # Increase the iterator of dataset 1
            index1 += 1
        elif reference['variants/POS'][index1] > query['variants/POS'][index2]:
            # If the position of the SNP in dataset 2 is smaller than the position of the SNP in dataset 1...
            # Increase the iterator of dataset 2
            index2 += 1
        else:
            # Found a SNP at the same position
            # Increase the counter of coincidences and the iterator of datasets 1 and 2
            n_coincidences += 1
            index1 += 1
            index2 += 1
    
    track('--> {} SNPs found at the same position between {} and {} datasets'.format(n_coincidences, dataset_name_1, dataset_name_2), track_name)
