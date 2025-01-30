import h5py
import pickle
import pandas as pd
from scipy.sparse import csr_matrix

def load_h5_pkl_to_df(h5_file_path, pkl_file_path):
    """
    Load HDF5 and PKL files and construct a DataFrame.

    Parameters:
    h5_file_path (str): Path to the HDF5 file (containing CSR matrix structure: data, indices, indptr, shape).
    pkl_file_path (str): Path to the PKL file (containing variants and barcodes information).

    Returns:
    pd.DataFrame: A DataFrame where rows are variants and columns are barcodes.
    """

    # Load CSR matrix from HDF5
    with h5py.File(h5_file_path, 'r') as f:
        data = f['data'][()]
        indices = f['indices'][()]
        indptr = f['indptr'][()]
        
        # Determine the shape (dataset or attribute)
        if 'shape' in f.keys():
            shape = tuple(f['shape'][()])
        else:
            shape = tuple(f.attrs['shape'])
    
    csr_mat = csr_matrix((data, indices, indptr), shape=shape)

    # Load variants and barcodes from PKL
    with open(pkl_file_path, 'rb') as f:
        mapping = pickle.load(f)

    if isinstance(mapping, dict):
        variants = mapping['variants']
        barcodes = mapping['barcodes']
    elif isinstance(mapping, (list, tuple)):
        variants, barcodes = mapping[0], mapping[1]
    else:
        raise TypeError("pkl file must be a dict or tuple/list containing variants and barcodes.")

    # Check if the CSR matrix matches the variants and barcodes
    if csr_mat.shape[0] != len(variants) or csr_mat.shape[1] != len(barcodes):
        raise ValueError(
            f"CSR matrix shape {csr_mat.shape} does not match "
            f"variants ({len(variants)}) or barcodes ({len(barcodes)})."
        )

    # Create DataFrame from CSR matrix
    df = pd.DataFrame.sparse.from_spmatrix(csr_mat, index=variants, columns=barcodes)
    return df

if __name__ == "__main__":
    # Example usage
    h5_file_path = "/path/to/your/alt_matrix.h5"
    pkl_file_path = "/path/to/your/variant_barcode_mappings.pkl"
    
    try:
        df = load_h5_pkl_to_df(h5_file_path, pkl_file_path)
        print("[INFO] DataFrame created successfully.")
        print(df.head())
    except Exception as e:
        print(f"[ERROR] {e}")
