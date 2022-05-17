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
import re
from utils.vcf_utils import combine_chrom_strands


def select_snps(query: dict, indexes: List) -> dict:
    """
    Selects the SNPs at the specified index positions.
    
    Args:
        query (dict): query data.
        indexes (List): list with the indexes of the SNPs to be kept.
                        The rest of SNPs are removed from the query.

    Returns:
        query (dict): query with the SNPs at the specified indexes.
    
    """

    # Keep all the SNPs at the specified index positions for all the keys 
    # different from samples
    for key in query.keys():
        if key != 'samples':
            query[key] = query[key][indexes]

    return query


def remove_snps(query: dict, indexes: List) -> dict:
    """
    Removes the SNPs at the specified index positions.
    
    Args:
        query (dict): query data.
        indexes (List): list with the indexes of the SNPs to be removed.
                        The rest of SNPs are kept in the query.
    
    Returns:
        query (dict): query without the SNPs at the specified indexes.
    
    """
    
    # Remove all the SNPs at the specified index positions for all the keys 
    # different from samples
    for key in query.keys():
        if key != 'samples':
            query[key] = np.delete(query[key], indexes, axis=0)
            
    return query


def remove_ambiguous_snps(query: dict, logger: logging.Logger) -> dict:
    """
    Searches and removes non-ambiguous SNPs from vcf data.
    
    Args:
        query (dict): allel.read_vcf output. Might contain ambiguous SNPs to be removed.
        Ambiguous SNPs can be A/T, T/A, C/G, or G/C pairs where the first element is the 
        reference (REF) and the second the alternate (ALT).
        logger (logging.Logger): debug/information tracker.
    
    Returns:
    
        query (dict): query without ambiguous SNPs.
    
    """
    
    non_ambiguous_idx = []  # empty list to store indexes of non-ambiguous SNPs
    
    A_T_count = 0  # counter of A/T ambiguities
    T_A_count = 0  # counter of T/A ambiguities
    C_G_count = 0  # counter of C/G ambiguities
    G_C_count = 0  # counter of G/C ambiguities

    # Obtain REF and ALT from query
    REF = query['variants/REF']
    ALT = query['variants/ALT']
    
    for i in range (0, len(REF)):
        # Check if it is an ambiguous SNP...
        if REF[i] == 'A' and ALT[i][0] == 'T':
            # If the REF is an A and the ALT is a T...
            # Found ambiguity of A/T type
            A_T_count += 1
        elif REF[i] == 'T' and ALT[i][0] == 'A':
            # If the REF is an T and the ALT is a A...
            # Found ambiguity of T/A type
            T_A_count += 1
        elif REF[i] == 'C' and ALT[i][0] == 'G':
            # If the REF is an C and the ALT is a G...
            # Found ambiguity of C/G type
            C_G_count += 1
        elif REF[i] == 'G' and ALT[i][0] == 'C':
            # If the REF is an G and the ALT is a C...
            # Found ambiguity of G/C type
            G_C_count += 1
        else:
            # The SNP at index i is not ambiguous
            non_ambiguous_idx.append(i)
    
    logger.info(f'--> {A_T_count} ambiguities found of A-T type.')
    logger.info(f'--> {T_A_count} ambiguities found of T-A type.')
    logger.info(f'--> {C_G_count} ambiguities found of C-G type.')
    logger.info(f'--> {G_C_count} ambiguities found of G-C type.')
    logger.info(f'--> {A_T_count+T_A_count+C_G_count+G_C_count} ambiguous SNPs in total.')
    
    # Keep the SNPs that are not ambiguous and
    # remove ambigous SNPs
    logger.debug('Removing all amgiguities.')
    query = select_snps(query, non_ambiguous_idx)
    
    return query


def correct_flips(query: dict, indexes: List) -> dict:
    """
    Corrects all flips in the dataset by swapping the REF and ALT fields.
    Consequently, the zeros and ones in the calldata/GT field in the 
    query are also swapped.
    
    Args:
        query (dict): query data.
        indexes (List): indexes of the SNPs with a flip.
    
    Returns:
        query (dict):query data with corrected SNP flips.
    
    """
    
    if len(indexes) > 0:
        # If there is a SNP flip...
        # Swap the REF and ALT fields
        ref_aux = query['variants/REF'][indexes]
        alt_aux = query['variants/ALT'][indexes,0]
        query['variants/REF'][indexes] = query['variants/ALT'][indexes,0]
        query['variants/ALT'][indexes,0] = ref_aux
        
        # Ensure the REF and the ALT were correctly swapped
        assert (query['variants/ALT'][indexes,0] == ref_aux).all(),\
        'The reference and the alternate were not swapped correctly.'
        assert (query['variants/REF'][indexes] == alt_aux).all(),\
        'The reference and the alternate were not swapped correctly.'
        
        # Change 0's by 1's and 1's by 0's
        snps = query['calldata/GT'][indexes,:,:]
        query['calldata/GT'][indexes,:,:] = np.where((snps==0)|(snps==1), snps^1, snps)
        
        # Ensure the 0's and 1's where correctly swapped
        # Check that the number of 1's before the swap is the same as the number of 0's after the swap
        number_ones_before = np.sum(snps == 1)
        number_zeros_after = np.sum(query['calldata/GT'][indexes,:,:] == 0)
        assert number_ones_before == number_zeros_after, 'An error occured while swapping the zeros and '\
        'ones in the variant/CHROM field.'
    
    return query


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
        reference (dict): reference data.
        query (dict): query data.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        query (dict): query with corrected SNP flips.
    
    """
    
    # Define counters of SNPs that at the same POS
    n_coincidences = 0
    
    flip_idx_query = []  # empty list to store the indexes of the SNPs with a SNP flip in the query
    
    i = 0  # reference iterator
    j = 0  # query iterator
    
    while i < len(reference['variants/POS']) and j < len(query['variants/POS']):
        if reference['variants/POS'][i] < query['variants/POS'][j]:
            # If the position of the SNP in reference is smaller than the position of the SNP in query...
            # Increase the iterator of reference
            i += 1
        elif reference['variants/POS'][i] > query['variants/POS'][j]:
            # If the position of the SNP in query is smaller than the position of the SNP in reference...
            # Increase the iterator of query
            j += 1
        else:   
            # Found a SNP at the same position...
            # Search if there is a SNP flip at this position
            if (reference['variants/REF'][i] == query['variants/ALT'][j,0] and 
            reference['variants/ALT'][i,0] == query['variants/REF'][j]):
                #The REF in the reference dataset is the ALT in the query dataset
                # and the ALT in the reference dataset is the REF in query dataset...
                # Found a SNP flip...
                # Append the index position of the SNP with a flip in query to flip_idx_query
                flip_idx_query.append(j)
            
            # Increase the counter of coincidences and the iterator of datasets 1 and 2
            n_coincidences += 1
            i += 1
            j += 1
    
    logger.info(f'--> {n_coincidences} SNPs found at the same position in the two datasets')
    logger.info(f'--> {len(flip_idx_query)} flips found in total')
    
    # Correct SNPs flips
    logger.debug('Correcting all SNP flips.')
    query = correct_flips(query, flip_idx_query)
    
    return query


def remove_mismatches_by_pos(reference: dict, query: dict, logger: logging.Logger) -> (dict, dict):
    """
    Searches and remove mismatches between reference and query. Note: there is a mismatch 
    when there is a difference in the reference (REF) or alternate (ALT) between SNPs at 
    the same position in two datasets (reference and query). Note: it is recommended to
    first correct the SNP flips prior to removing SNP mismatches. You can correct the SNP 
    flips with correct_flips_by_pos.
    
    Args:
        reference (dict): reference data.
        query (dict): query data.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        reference (dict): reference data without SNP mismatches.
        query (dict): query data without SNP mismatches.
    
    """
    
    # Define counters of SNPs that are at the same position (same POS)
    n_coincidences = 0
    
    mismatch_idx_reference = []  # empty list to store the indexes of the SNPs that contain a SNP mismatch in the reference
    mismatch_idx_query = []  # empty list to store the indexes of the SNPs that contain a SNP mismatch in the query
    
    i = 0  # reference iterator
    j = 0  # query iterator
    
    while i < len(reference['variants/POS']) and j < len(query['variants/POS']):
        if reference['variants/POS'][i] < query['variants/POS'][j]:
            # If the position of the SNP in reference is smaller than the position of the SNP in query...
            # Increase the iterator of reference
            i += 1
        elif reference['variants/POS'][i] > query['variants/POS'][j]:
            # If the position of the SNP in query is smaller than the position of the SNP in reference...
            # Increase the iterator of query
            j += 1
        else:
            # Found a SNP at the same position...
            # Search if there is a mismatch
            if (reference['variants/REF'][i] != query['variants/REF'][j] or 
                reference['variants/ALT'][i,0] != query['variants/ALT'][j,0]):
                # If the reference or the alternate between datasets 1 and 2 do not match...
                # Save the index positions of the SNPs that present a mismatch in each dataset
                mismatch_idx_reference.append(i)
                mismatch_idx_query.append(j)
            
            # Increase number of matches
            n_coincidences += 1
            
            # Increase the reference and query counters
            i += 1
            j += 1
    
    logger.info(f'--> {n_coincidences} SNPs found at the same position in the two datasets')
    logger.info(f'--> {len(mismatch_idx_reference)} mismatches found in total')
  
    # Remove mismatching SNPs
    logger.debug('Correcting all SNP mismatches.')
    reference = remove_snps(reference, mismatch_idx_reference)
    query = remove_snps(query, mismatch_idx_query)
    
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
        reference (dict): reference data.
        query (dict): query data.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        reference (dict): reference for the common markers.
        query (dict): query for the common markers.
        idxs_reference (List): indexes of the SNPs in the .vcf file 1 that are common markers with the .vcf file 2.
        idxs_query (List): indexes of the SNPs in the .vcf file 2 that are common markers with the .vcf file 1.
    
    """
    
    idxs_reference = []  # empty list to store indexes of the common markers in the reference
    idxs_query = []  # empty list to store indexes of the common markers in the query
    
    i = 0  # reference iterator
    j = 0  # query iterator
    
    while i < len(reference['variants/POS']) and j < len(query['variants/POS']):
        if reference['variants/POS'][i] < query['variants/POS'][j]:
            # If the position of the SNP in the reference is smaller than the position 
            # of the SNP in the reference...
            # Increase reference counter
            i += 1
        elif reference['variants/POS'][i] > query['variants/POS'][j]:
            # If the position of the SNP in the query is smaller than the position 
            # of the SNP in the query...
            # Increase query counter
            j += 1
        else:
            # Found a SNP at the same position...
            # Search if it is a common marker
            if (reference['variants/REF'][i] == query['variants/REF'][j] and 
                reference['variants/ALT'][i][0] == query['variants/ALT'][j][0]):
                # If the reference and the alternate between both datasets match...
                # A common marker was found
                # Save the indexes of the common marker in both datasets
                idxs_reference.append(i)
                idxs_query.append(j)
            
            # Increase the reference and query counters
            i += 1
            j += 1
    
    # Keep the SNPs that are common markers and remove the rest
    reference = select_snps(reference, idxs_reference)
    query = select_snps(query, idxs_query)
     
    logger.info(f'--> {len(idxs_reference)} common markers in total.')
    
    return reference, query, idxs_reference, idxs_query


def keep_common_markers_several_chr(reference: dict, query: dict, logger: logging.Logger):
    """
    Searches and keeps common markers between two datasets (reference and query).
    Note: common markers are SNPs with same chromosome (CHROM), position (POS), reference (REF), 
    and alernate (ALT). This function accepts reference and query data for multiple chromosomes.
    
    Args:
        reference (dict): reference data.
        query (dict): query data.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        reference (dict): reference for the common markers.
        query (dict): query for the common markers.
        idxs_reference (List): indexes of the SNPs in the .vcf file 1 that are common markers with the .vcf file 2.
        idxs_query (List): indexes of the SNPs in the .vcf file 2 that are common markers with the .vcf file 1.
    
    """
    
    idxs_reference = []  # empty list to store indexes of the common markers in the reference
    idxs_query = []  # empty list to store indexes of the common markers in the query
    
    i = 0  # reference iterator
    j = 0  # query iterator

    while i < len(reference['variants/POS']) and j < len(query['variants/POS']):
        
        # Obtain chromosome number
        if reference['variants/CHROM'][i].isdigit(): chrom_ref = reference['variants/CHROM'][i]
        else: chrom_ref = int(re.search(r'\d+', reference['variants/CHROM'][i]).group())
        if query['variants/CHROM'][j].isdigit(): chrom_query = query['variants/CHROM'][j]
        else: chrom_query = int(re.search(r'\d+', query['variants/CHROM'][j]).group())
        
        if chrom_ref < chrom_query:
            i += 1
        elif chrom_ref > chrom_query:
            j += 1
        elif reference['variants/POS'][i] < query['variants/POS'][j]:
            # If the position of the SNP in the reference is smaller than the position 
            # of the SNP in the reference...
            # Increase reference counter
            i += 1
        elif reference['variants/POS'][i] > query['variants/POS'][j]:
            # If the position of the SNP in the query is smaller than the position 
            # of the SNP in the query...
            # Increase query counter
            j += 1
        else:
            # Found a SNP at the same position...
            # Search if it is a common marker
            if (reference['variants/REF'][i] == query['variants/REF'][j] and 
                reference['variants/ALT'][i][0] == query['variants/ALT'][j][0]):
                # If the reference and the alternate between datasets 1 and 2 match...
                # Save the index of the SNPs that are a common marker in each dataset
                idxs_reference.append(i)
                idxs_query.append(j)

                # Increase the reference and query counters
                i += 1
                j += 1
                
    logger.info(f'--> {len(idxs_reference)} common markers in total.')
    
    # Keep the SNPs that are common markers and remove the rest
    reference = select_snps(reference, idxs_reference)
    query = select_snps(query, idxs_query)
    
    return reference, query, idxs_reference, idxs_query


def remove_snps_different_means(reference: dict, query: dict, threshold: float, 
                                logger: logging.Logger) -> (dict, dict):
    """
    Searches the common markers (i.e., SNPs at the same CHROM, POS, REF, and ALT) 
    between the query and the reference and removes the ones with a mean absolute 
    difference higher than a threshold.
    
    Args:
        reference (dict): reference data.
        query (dict): query data.
        threshold (float): value between 0.0 and 1.0. All SNPs with a mean absolute difference 
                           higher than this value will be removed.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        reference (dict): reference data without the SNPs with a mean absolute 
        difference higher than a threshold.
        query (dict): query data without the SNPs with a mean absolute difference 
        higher than a threshold.
    
    """

    # Obtain SNPs data for datasets 1 and 2
    snps_reference = reference['calldata/GT']
    snps_query = query['calldata/GT']
        
    # Convert SNPs from MxNx2 to MxN shape (averaged)
    # For the reference
    length_reference, num_dogs_reference, num_strands_reference = snps_reference.shape
    snps_reference = snps_reference.reshape(length_reference, num_dogs_reference*2).T
    snps_reference = combine_chrom_strands(snps_reference)
    
    # Convert SNPs from MxNx2 to MxN shape (averaged)
    # For the query
    length_query, num_dogs_query, num_strands_query = snps_query.shape
    snps_query = snps_query.reshape(length_query, num_dogs_query*2).T
    snps_query = combine_chrom_strands(snps_query)
    
    # Search indexes of common markers
    _, _, idxs_reference, idxs_query = keep_common_markers_single_chr(reference.copy(), query.copy(), logger)
    
    # Filter common markers
    snps_reference = snps_reference[:, idxs_reference]
    snps_query = snps_query[:, idxs_query]
    
    # Compute the mean of common markers
    mean_reference = np.mean(snps_reference, axis=0)
    mean_query = np.mean(snps_query, axis=0)
    
    # Compute the indexes of the SNPs to be removed in datasets 1 and 2
    idx_reference_to_remove = np.array(idxs_reference)[abs(mean_reference - mean_query) > threshold]
    idx_query_to_remove = np.array(idxs_query)[abs(mean_reference - mean_query) > threshold]
    
    assert len(idx_reference_to_remove) == len(idx_query_to_remove) ==\
    sum(abs(mean_reference - mean_query) > threshold),\
    'The indexes removed in each dataset is not the same. Terminate.'
    
    logger.info(f'Found {sum(abs(mean_reference - mean_query) > threshold)} SNPs with a mean '\
                f' absolute difference higher than {threshold}.')
    
    # Remove the SNPs that have a mean absolute difference than threshold
    # from the reference and the query
    reference = remove_snps(reference, idx_reference_to_remove)
    query = remove_snps(query, idx_query_to_remove)
    
    return reference, query
