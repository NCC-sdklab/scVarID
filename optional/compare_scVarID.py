#%%

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import h5py
import pickle
import pandas as pd
from scipy.sparse import csr_matrix

########################################################################
# 1. H5 + PKL => DataFrame
########################################################################

def load_h5_pkl_to_df(h5_file_path, pkl_file_path):
    """
    HDF5의 sparse CSR 정보를 읽고,
    PKL에서 (variants, barcodes) 매핑을 불러와,
    (variants × barcodes) 형태의 DataFrame을 생성.
    """
    # H5에서 CSR read
    with h5py.File(h5_file_path, 'r') as f:
        data = f['data'][()]
        indices = f['indices'][()]
        indptr = f['indptr'][()]

        # shape
        if 'shape' in f.keys():
            shape = tuple(f['shape'][()])
        else:
            shape = tuple(f.attrs['shape'])

    csr_mat = csr_matrix((data, indices, indptr), shape=shape)

    # PKL에서 (variants, barcodes) 불러옴
    with open(pkl_file_path, 'rb') as f:
        mapping = pickle.load(f)

    # mapping이 (variants, barcodes) 형태인지 검사
    if isinstance(mapping, dict):
        variants = mapping['variants']
        barcodes = mapping['barcodes']
    elif isinstance(mapping, (list, tuple)):
        variants, barcodes = mapping[0], mapping[1]
    else:
        raise TypeError("pkl file must be a dict or tuple/list containing variants and barcodes.")

    # shape 검증
    if csr_mat.shape[0] != len(variants) or csr_mat.shape[1] != len(barcodes):
        raise ValueError(
            f"CSR matrix shape {csr_mat.shape} != (#variants={len(variants)}, #barcodes={len(barcodes)})"
        )

    df = pd.DataFrame.sparse.from_spmatrix(csr_mat, index=variants, columns=barcodes)
    return df

########################################################################
# 2. 차이점 확인 함수
########################################################################

def compare_scVarID_matrices(df_old, df_new, label=""):
    """
    1) shape 비교
    2) variant(행) 비교
    3) barcode(열) 비교
    4) 교집합에서 element-wise 비교 -> 최대 10개 정도 출력
    """
    print(f"\n=== Compare {label} matrix ===")

    # 1) shape
    print(f"[OLD] shape = {df_old.shape}, [NEW] shape = {df_new.shape}")

    # 2) 행(variants) 비교
    old_variants = set(df_old.index)
    new_variants = set(df_new.index)
    common_variants = old_variants & new_variants
    only_old_variants = old_variants - new_variants
    only_new_variants = new_variants - old_variants

    print(f" - old variants: {len(old_variants)}, new variants: {len(new_variants)}, common={len(common_variants)}")
    if only_old_variants:
        print(f"   > variants in OLD only: {len(only_old_variants)}")
    if only_new_variants:
        print(f"   > variants in NEW only: {len(only_new_variants)}")

    # 3) 열(barcodes) 비교
    old_barcodes = set(df_old.columns)
    new_barcodes = set(df_new.columns)
    common_barcodes = old_barcodes & new_barcodes
    only_old_barcodes = old_barcodes - new_barcodes
    only_new_barcodes = new_barcodes - old_barcodes

    print(f" - old barcodes: {len(old_barcodes)}, new barcodes: {len(new_barcodes)}, common={len(common_barcodes)}")
    if only_old_barcodes:
        print(f"   > barcodes in OLD only: {len(only_old_barcodes)}")
    if only_new_barcodes:
        print(f"   > barcodes in NEW only: {len(only_new_barcodes)}")

    # 4) 교집합에서 값 비교
    if not common_variants or not common_barcodes:
        print("[WARN] No common variant or barcode => skip value comparison.")
        return

    df_old_common = df_old.loc[list(common_variants), list(common_barcodes)].fillna(0)
    df_new_common = df_new.loc[list(common_variants), list(common_barcodes)].fillna(0)

    diff_matrix = (df_old_common - df_new_common).abs()
    nonzero_count = (diff_matrix != 0).sum().sum()

    if nonzero_count > 0:
        print(f" - There are {nonzero_count} cells with different values in the intersection.")
        # 최대 10개만 표시
        show_diff_in_console(diff_matrix, df_old_common, df_new_common, label)
    else:
        print(f" - All intersecting cells have identical values.")

def show_diff_in_console(diff_matrix, df_old_common, df_new_common, label):
    """
    diff_matrix에서 != 0인 항목 최대 30개만 출력.
    """
    diff_mask = (diff_matrix != 0)
    coords = diff_mask.stack()  # Series of booleans
    diff_positions = coords[coords].index.tolist()  # (variant, barcode)
    
    print(f"[INFO] Listing up to 30 differences in {label} matrix:")
    count_shown = 0
    for (variant, barcode) in diff_positions:
        old_val = df_old_common.loc[variant, barcode]
        new_val = df_new_common.loc[variant, barcode]
        diff_val= abs(old_val - new_val)
        print(f"   variant={variant}, barcode={barcode}, old={old_val}, new={new_val}, diff={diff_val}")
        count_shown +=1
        if count_shown>=30:
            break

########################################################################
# 3. main()
########################################################################

def main():
    """
    OLD: /path/old_dir
    NEW: /path/new_dir
    특정 4개(matrix) => [ref, alt, missing, unknown]
    """
    # 디렉토리 설정
    old_dir = "/mnt/dks_nas2/data/raw/Breast/GSE246613/Processed/SNV/scVarID/h01A/chr1"
    new_dir = "/mnt/dks_nas2/data/raw/Breast/GSE246613/Processed/SNV/scVarID/h01A/chr1_window26"

    # 비교할 행렬들
    matrices = ["ref","alt","missing","unknown"]

    for m in matrices:
        h5_old= os.path.join(old_dir, f"{m}_matrix.h5")
        h5_new= os.path.join(new_dir, f"{m}_matrix.h5")
        pkl_old= os.path.join(old_dir, "variant_barcode_mappings.pkl")
        pkl_new= os.path.join(new_dir, "variant_barcode_mappings.pkl")

        # 파일 존재여부
        if not (os.path.exists(h5_old) and os.path.exists(h5_new)):
            print(f"[SKIP] {m}_matrix.h5 not found in old/new directory => skip.")
            continue

        # DF 로드
        df_old = load_h5_pkl_to_df(h5_old, pkl_old)
        df_new = load_h5_pkl_to_df(h5_new, pkl_new)

        # 비교
        compare_scVarID_matrices(df_old, df_new, label=m)

if __name__=="__main__":
    main()

# %%
