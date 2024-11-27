# utils.py

import multiprocessing
import logging
import h5py
from collections import defaultdict
import joblib
from scipy.sparse import csr_matrix

def check_num_cores(num_cores):
    max_cores = multiprocessing.cpu_count()
    if num_cores > max_cores:
        print(f"Warning: Specified number of cores ({num_cores}) exceeds available cores ({max_cores}). Using {max_cores} cores instead.")
        num_cores = max_cores
    return num_cores

def read_barcodes_from_file(barcode_path):
    """
    Reads barcodes from a file if barcode_path is provided, else returns None.
    """
    if barcode_path:
        with open(barcode_path, 'r') as f:
            barcodes = [line.strip() for line in f if line.strip()]
        print(f"Total barcodes from file: {len(barcodes)}")
        return barcodes
    else:
        return None  # Return None if barcode file is not provided

def save_classification_matrices(
    save_dir, 
    union_variants, 
    barcode_path=None, 
    ref_classifications=None, 
    alt_classifications=None, 
    missing_classifications=None, 
    unknown_classifications=None
):
    """
    Save classified variant data.
    """
    # Extract variant list
    variants = [f"{variant['CHROM']}:{variant['POS']}_{variant['REF']}->{variant['ALT']}" for variant in union_variants]
    
    print(f"Total variants: {len(variants)}")
    
    # Create mapping dictionaries
    variant_to_index = {variant: idx for idx, variant in enumerate(variants)}
    
    # Prepare classifications dictionary
    classifications_dict = {}
    if ref_classifications:
        classifications_dict['ref'] = ref_classifications
    if alt_classifications:
        classifications_dict['alt'] = alt_classifications
    if missing_classifications:
        classifications_dict['missing'] = missing_classifications
    if unknown_classifications:
        classifications_dict['unknown'] = unknown_classifications

    # Read barcodes from file or generate from data
    barcodes = read_barcodes_from_file(barcode_path)
    if barcodes is None:
        barcode_set = set()
        for classifications in classifications_dict.values():
            for classification in classifications:
                classification_dict = dict(classification)
                barcode = classification_dict['cell_barcode']
                barcode_set.add(barcode)
        barcodes = sorted(barcode_set)
        print(f"Total barcodes generated from data: {len(barcodes)}")
    
    barcode_to_index = {barcode: idx for idx, barcode in enumerate(barcodes)}
    
    # Generate and save sparse matrices
    classification_matrices = {}
    
    for classification_name, classifications in classifications_dict.items():
        if not classifications:
            print(f"No data for classification '{classification_name}'. Skipping.")
            continue  # Skip if no data
        
        data = []
        row_indices = []
        col_indices = []
        count_dict = defaultdict(int)
    
        for classification in classifications:
            classification_dict = dict(classification)
            variant = classification_dict['variant']
            barcode = classification_dict['cell_barcode']
            
            # Check if variant and barcode are in mapping dictionaries
            if variant in variant_to_index and barcode in barcode_to_index:
                variant_idx = variant_to_index[variant]
                barcode_idx = barcode_to_index[barcode]
                key = (variant_idx, barcode_idx)
                count_dict[key] += 1  # Increment read count
            else:
                logging.warning(f"Variant '{variant}' or barcode '{barcode}' not found in mapping dictionaries.")
                continue
    
        for (variant_idx, barcode_idx), count in count_dict.items():
            row_indices.append(variant_idx)
            col_indices.append(barcode_idx)
            data.append(count)
    
        # Create sparse matrix
        matrix_shape = (len(variants), len(barcodes))
        matrix = csr_matrix((data, (row_indices, col_indices)), shape=matrix_shape)
        classification_matrices[classification_name] = matrix
    
        # Save sparse matrix to h5 file
        with h5py.File(f'{save_dir}/{classification_name}_matrix.h5', 'w') as h5f:
            h5f.create_dataset('data', data=matrix.data)
            h5f.create_dataset('indices', data=matrix.indices)
            h5f.create_dataset('indptr', data=matrix.indptr)
            h5f.attrs['shape'] = matrix.shape
    
    # Save variant and barcode lists
    joblib.dump((variants, barcodes), f'{save_dir}/variant_barcode_mappings.pkl')
