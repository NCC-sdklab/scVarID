# classification.py

import logging
import re
from collections import defaultdict
from joblib import Parallel, delayed
import pysam

# Import necessary functions from classification_functions.py
from classification_functions import (
    parse_cigar,
    classify_variant,
    create_variants_dict,
    prime_read,
    handle_snv,
    handle_insertion,
    handle_deletion
)

def split_list(lst, num_sublists):
    avg_size = len(lst) // num_sublists
    remainder = len(lst) % num_sublists
    result = []
    start = 0
    for i in range(num_sublists):
        end = start + avg_size + (1 if i < remainder else 0)
        result.append(lst[start:end])
        start = end
    return result

def process_read(read, variants_dict, compute_missing_unknown=True):
    ref_classifications_local = set()
    alt_classifications_local = set()
    missing_classifications_local = set()
    unknown_classifications_local = set()

    read_name, cell_barcode, sequences, positions, operations = prime_read(read)
    chrom = read.reference_name

    if cell_barcode is None:
        logging.debug(f"Read {read_name} has no CB tag. Assigning to 'unknown'.")
        # 모든 변이에 대해 'unknown'으로 분류
        for key in variants_dict:
            for variant_data in variants_dict[key]:
                variant = variant_data['variant']
                variant_type = variant_data['type']
                classification, bc, p, rn = classify_variant(
                    sequences, positions, operations,
                    variant[1], variant[2], variant[3],
                    variant_type, read_name, 'unknown'
                )
                classification_data = {
                    'read_name': rn,
                    'cell_barcode': 'unknown',
                    'variant': f"{variant[0]}:{variant[1]}_{variant[2]}->{variant[3]}"
                }
                classification_frozenset = frozenset(classification_data.items())
                if classification == 'ref':
                    ref_classifications_local.add(classification_frozenset)
                elif classification == 'alt':
                    alt_classifications_local.add(classification_frozenset)
                elif classification == 'missing' and compute_missing_unknown:
                    missing_classifications_local.add(classification_frozenset)
                elif classification == 'unknown' and compute_missing_unknown:
                    unknown_classifications_local.add(classification_frozenset)
        return ref_classifications_local, alt_classifications_local, missing_classifications_local, unknown_classifications_local

    visited_positions = set()

    for pos in positions:
        key = (chrom, pos)
        if key in variants_dict and key not in visited_positions:
            visited_positions.add(key)
            for variant_data in variants_dict[key]:
                variant = variant_data['variant']
                variant_type = variant_data['type']
                classification, bc, p, rn = classify_variant(
                    sequences, positions, operations,
                    variant[1], variant[2], variant[3],
                    variant_type, read_name, cell_barcode
                )
                classification_data = {
                    'read_name': rn,
                    'cell_barcode': bc,
                    'variant': f"{variant[0]}:{variant[1]}_{variant[2]}->{variant[3]}"
                }
                classification_frozenset = frozenset(classification_data.items())
                if classification == 'ref':
                    ref_classifications_local.add(classification_frozenset)
                elif classification == 'alt':
                    alt_classifications_local.add(classification_frozenset)
                elif classification == 'missing' and compute_missing_unknown:
                    missing_classifications_local.add(classification_frozenset)
                elif classification == 'unknown' and compute_missing_unknown:
                    unknown_classifications_local.add(classification_frozenset)

    return ref_classifications_local, alt_classifications_local, missing_classifications_local, unknown_classifications_local

def process_read_chunk(
    bam_file_path, 
    read_unique_names_set, 
    variants_dict, 
    read_mapping, 
    compute_missing_unknown=True
):
    ref_classifications = []
    alt_classifications = []
    missing_classifications = []
    unknown_classifications = []

    try:
        with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
            read_name_to_unique = {}
            for read in bam_file.fetch(until_eof=True):
                read_name = read.query_name
                if read_name not in read_name_to_unique:
                    read_name_to_unique[read_name] = 0
                else:
                    read_name_to_unique[read_name] += 1

                unique_id = read_name_to_unique[read_name]
                if read_name in read_mapping and unique_id < len(read_mapping[read_name]):
                    read_unique_name = read_mapping[read_name][unique_id]
                else:
                    read_unique_name = f"{read_name}_{unique_id}"

                if read_unique_name in read_unique_names_set:
                    ref, alt, missing, unknown = process_read(
                        read, variants_dict, compute_missing_unknown
                    )
                    ref_classifications.extend(ref)
                    alt_classifications.extend(alt)
                    missing_classifications.extend(missing)
                    unknown_classifications.extend(unknown)

                    read_unique_names_set.remove(read_unique_name)
                    if not read_unique_names_set:
                        break
    except Exception as e:
        logging.error(f"Error processing read chunk: {e}", exc_info=True)

    return ref_classifications, alt_classifications, missing_classifications, unknown_classifications

def process_bam_data_parallel(
    bam_path, 
    selected_read_unique_names, 
    selected_variants, 
    read_mapping, 
    num_cores=4, 
    compute_missing_unknown=True
):
    """
    Process BAM data in parallel and return classification results.
    """
    variants_dict = create_variants_dict(selected_variants)
    read_unique_names_chunks = split_list(selected_read_unique_names, num_cores)

    # Initialize result lists
    ref_classifications = []
    alt_classifications = []
    missing_classifications = []
    unknown_classifications = []

    # Parallel processing
    results = Parallel(n_jobs=num_cores, backend='multiprocessing')(
        delayed(process_read_chunk)(
            bam_path, 
            set(chunk), 
            variants_dict, 
            read_mapping, 
            compute_missing_unknown
        ) for chunk in read_unique_names_chunks
    )

    for result in results:
        ref_classifications.extend(result[0])
        alt_classifications.extend(result[1])
        missing_classifications.extend(result[2])
        unknown_classifications.extend(result[3])

    return ref_classifications, alt_classifications, missing_classifications, unknown_classifications
