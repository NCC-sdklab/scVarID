# scVarID_window_padding_chunk_parallel_final_integrated_ver6.py

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
scVarID Pipeline - 버전 6
- Chunk-level Parallel Processing
- Region-based VCF/BAM Fetch
- 정확한 Read Parsing 및 Classification Logic
- Temporary Intermediate Files 관리 후 최종 HDF5 저장
"""

import os
import sys
import time
import logging
import argparse
import re
from collections import defaultdict
import pysam
import h5py
import joblib
import pandas as pd
from scipy.sparse import csr_matrix
from joblib import Parallel, delayed

# 임포트한 분류 함수들
from classification_functions import (
    parse_cigar,
    classify_variant,
    create_variants_dict,
    prime_read
)

############################################################################
# 0. Step Logging
############################################################################
def log_step_start(step_name):
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    end_time = time.time()
    elapsed = end_time - start_time
    hh = int(elapsed // 3600)
    mm = int((elapsed % 3600) // 60)
    ss = int(elapsed % 60)
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {hh:02d}h {mm:02d}m {ss:02d}s ===")

############################################################################
# 1. Argument Parsing & Logging Setup
############################################################################
def parse_arguments():
    parser = argparse.ArgumentParser(description='scVarID Integrated Pipeline with Window Processing and Parallelism.')
    parser.add_argument('--variant-files', nargs='+', required=True,
                        help='Pairs: SAMPLE, FILE (e.g. sample1 file1.vcf sample2 file2.vcf)')
    parser.add_argument('--bam-path', required=True, help='Path to the BAM file.')
    parser.add_argument('--save-dir', required=True, help='Directory to save output files.')
    parser.add_argument('--barcode-path', required=False, help='(Optional) Path to the barcode file.')
    parser.add_argument('--num-cores', type=int, default=4, help='Number of parallel processes.')
    parser.add_argument('--ref-alt-only', action='store_true', help='If set, skip missing/unknown classifications.')
    parser.add_argument('--chromosomes', nargs='+', required=False, help='List of chromosomes to process.')
    parser.add_argument('--window-size', type=int, default=10_000_000, help='Size of each window (bp).')
    parser.add_argument('--padding', type=int, default=5000, help='Padding around window boundaries (bp).')
    args = parser.parse_args()

    if len(args.variant_files) % 2 != 0:
        parser.error("Each sample must be paired with a VCF/TXT file in --variant-files.")

    vcf_files_with_origins = []
    for i in range(0, len(args.variant_files), 2):
        sample = args.variant_files[i]
        file_path = args.variant_files[i+1]
        vcf_files_with_origins.append((file_path, sample))

    return args, vcf_files_with_origins

def setup_logging(save_dir):
    os.makedirs(save_dir, exist_ok=True)
    log_file = os.path.join(save_dir, 'processing.log')
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.DEBUG,  # 디버깅을 위해 DEBUG로 설정
        format='%(asctime)s:%(levelname)s:%(message)s'
    )
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)  # 콘솔은 INFO 레벨
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

############################################################################
# 2. Variant Processing
############################################################################
def is_vcf_file(file_path):
    """Check if the file is in VCF or BCF format."""
    return file_path.endswith('.vcf') or file_path.endswith('.vcf.gz') or file_path.endswith('.bcf')

def extract_variants_vcf_region(vcf_file, origin, chrom, start, end):
    """Extract variants from a VCF file within a specific region."""
    variants = []
    try:
        with pysam.VariantFile(vcf_file) as vf:
            for record in vf.fetch(chrom, start, end):
                if "PASS" in record.filter.keys():
                    variant = {
                        "CHROM": record.chrom,
                        "POS": record.pos,
                        "REF": record.ref,
                        "ALT": ",".join(record.alts) if record.alts else "",
                        "QUAL": record.qual if record.qual is not None else ".",
                        "FILTER": ";".join(record.filter.keys()) if record.filter else ".",
                        "ORIGIN": origin
                    }
                    variants.append(variant)
    except Exception as e:
        logging.warning(f"Error extracting variants from VCF {vcf_file}: {e}")
    logging.info(f"Extracted {len(variants)} variants from {vcf_file} in region {chrom}:{start}-{end}")
    return variants

def extract_variants_txt_region(txt_file, origin, chrom, start, end):
    """Extract variants from a pre-filtered TXT file within a specific region."""
    variants = []
    try:
        df = pd.read_csv(txt_file, sep='\t', header=None, names=["CHROM", "POS", "REF", "ALT"],
                         dtype={"CHROM": str, "POS": int, "REF": str, "ALT": str})
        if df.shape[1] != 4:
            raise ValueError(f"TXT file {txt_file} must have exactly 4 columns.")
        subset = df[(df["CHROM"] == chrom) & (df["POS"] >= start) & (df["POS"] <= end)]
        for _, row in subset.iterrows():
            variant = {
                "CHROM": row["CHROM"],
                "POS": row["POS"],
                "REF": row["REF"],
                "ALT": row["ALT"],
                "QUAL": ".",
                "FILTER": ".",
                "ORIGIN": origin
            }
            variants.append(variant)
    except Exception as e:
        logging.warning(f"Error extracting variants from TXT {txt_file}: {e}")
    logging.info(f"Extracted {len(variants)} variants from {txt_file} in region {chrom}:{start}-{end}")
    return variants

def filter_variants(variants):
    """Filter variants based on specific conditions and merge origins."""
    union_map = {}
    for variant in variants:
        ref = variant["REF"]
        alt = variant["ALT"]
        key = (variant["CHROM"], variant["POS"], ref, alt)

        # Exclude variants with multiple REF or ALT alleles or mismatched first characters
        if "," in ref or "," in alt:
            continue
        if (len(ref) > 1 or len(alt) > 1) and ref[0] != alt[0]:
            continue

        if key not in union_map:
            union_map[key] = variant
        else:
            existing_origins = union_map[key]["ORIGIN"].split(", ")
            new_origin = variant["ORIGIN"]
            combined_origins = existing_origins + [new_origin]
            # Maintain desired order without duplicates
            ordered_origins = []
            for origin in ['normal', 'tumor', 'SCREADCOUNT']:
                if origin in combined_origins and origin not in ordered_origins:
                    ordered_origins.append(origin)
            union_map[key]["ORIGIN"] = ", ".join(ordered_origins)

    # Assign VINDEX
    unioned_variants = list(union_map.values())
    for idx, variant in enumerate(unioned_variants, 1):
        variant['VINDEX'] = f'vi_{idx}'

    logging.info(f"Filtered variants count: {len(unioned_variants)}")
    return unioned_variants

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    """Process VCF/TXT files within a specified region and return filtered variants."""
    all_variants = []
    for (file_path, origin) in vcf_files_with_origins:
        if is_vcf_file(file_path):
            logging.info(f"Processing VCF/BCF file: {file_path}")
            variants = extract_variants_vcf_region(file_path, origin, chrom, start, end)
        elif file_path.endswith('.txt'):
            logging.info(f"Processing TXT file: {file_path}")
            variants = extract_variants_txt_region(file_path, origin, chrom, start, end)
        else:
            logging.warning(f"Unsupported variant file format: {file_path}. Skipping.")
            variants = []
        all_variants.extend(variants)
    
    union_variants = filter_variants(all_variants)
    return union_variants

############################################################################
# 3. Read Processing & Classification Logic
############################################################################
def process_read(read, variants_dict, compute_missing_unknown=True):
    """Process a single read and classify variants."""
    ref_classifications_local = set()
    alt_classifications_local = set()
    missing_classifications_local = set()
    unknown_classifications_local = set()

    read_name, cell_barcode, sequences, positions, operations = prime_read(read)
    chrom = read.reference_name

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
    """Process a chunk of reads and classify variants."""
    ref_classifications = []
    alt_classifications = []
    missing_classifications = []
    unknown_classifications = []

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

def split_list(lst, num_sublists):
    """Split a list into nearly equal sublists."""
    avg_size = len(lst) // num_sublists
    remainder = len(lst) % num_sublists
    result = []
    start = 0
    for i in range(num_sublists):
        end = start + avg_size + (1 if i < remainder else 0)
        result.append(lst[start:end])
        start = end
    return result

############################################################################
# 4. Chunk Processing
############################################################################
def extract_reads_with_padding(bam_path, chrom, start, end, padding=5000):
    """Extract reads from BAM within a window with padding."""
    data = []
    read_name_counts = defaultdict(int)
    fs = max(1, start - padding)
    fe = end + padding

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(chrom, fs - 1, fe):
                read_name = read.query_name
                try:
                    bc = read.get_tag("CB")
                except KeyError:
                    bc = None
                unique_id = read_name_counts[read_name]
                read_unique = f"{read_name}_{unique_id}"
                read_name_counts[read_name] += 1

                data.append([
                    read_name, read_unique, bc, chrom,
                    read.reference_start, read.reference_end
                ])
    except Exception as e:
        logging.error(f"Error extracting reads from BAM: {e}", exc_info=True)

    df_reads = pd.DataFrame(data, columns=[
        "Read_Name", "Read_Unique_Name", "Cell_Barcode", "Chromosome",
        "Start_Pos", "End_Pos"
    ])
    logging.info(f"Extracted {len(df_reads)} reads from {chrom}:{start}-{end} in {bam_path}")
    return df_reads

def process_one_chunk_and_dump(
    chrom, win_start, win_end,
    vcf_files_with_origins,
    bam_path,
    padding=5000,
    compute_missing_unknown=True,
    tmp_dir=None
):
    """Process a single window chunk and save intermediate results."""
    chunk_str = f"{chrom}:{win_start}-{win_end}"
    logging.info(f"[CHUNK] Start {chunk_str}")

    # 1. Extract variants within the window
    union_variants = process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(union_variants)} variants")
    if not union_variants:
        return None

    # 2. Extract reads within the window with padding
    df_reads = extract_reads_with_padding(bam_path, chrom, win_start, win_end, padding=padding)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(df_reads)} reads")
    if df_reads.empty:
        return None

    # 3. Build read mapping
    read_mapping = df_reads.groupby("Read_Name")["Read_Unique_Name"].apply(list).to_dict()
    read_uniq_list = df_reads["Read_Unique_Name"].tolist()

    # 4. Filter variants strictly within the window
    final_variants = [v for v in union_variants if (win_start <= v["POS"] <= win_end)]
    variants_dict = create_variants_dict(final_variants)

    # 5. Classify reads
    refC, altC, missC, unkC = process_bam_data_parallel(
        bam_path,
        read_uniq_list,
        final_variants,
        read_mapping,
        num_cores=1,  # 이미 파이프라인 레벨에서 병렬화
        compute_missing_unknown=compute_missing_unknown
    )

    # 6. Log statistics
    m_reads = len(read_uniq_list)
    t_reads = len(df_reads)
    p_reads = (100.0 * m_reads / t_reads) if t_reads > 0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable reads: {m_reads} ({p_reads:.2f}%)")

    m_vars = len(final_variants)
    t_vars = len(union_variants)
    p_vars = (100.0 * m_vars / t_vars) if t_vars > 0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable variants: {m_vars} ({p_vars:.2f}%)")

    # 7. Save intermediate chunk data
    chunk_data = {
        "union_variants": final_variants,
        "ref": refC,
        "alt": altC,
        "missing": missC,
        "unknown": unkC
    }
    if not tmp_dir:
        tmp_dir = "./tmp_chunks"
    os.makedirs(tmp_dir, exist_ok=True)
    out_pkl = os.path.join(tmp_dir, f"chunk_{chrom}_{win_start}_{win_end}.pkl")
    joblib.dump(chunk_data, out_pkl)
    logging.info(f"[CHUNK] {chunk_str} => {out_pkl}")
    return out_pkl

############################################################################
# 5. Save Final Classification Matrices
############################################################################
def read_barcodes_from_file(barcode_path):
    """Read barcodes from a file if provided."""
    if barcode_path:
        try:
            with open(barcode_path, 'r') as f:
                barcodes = [line.strip() for line in f if line.strip()]
            logging.info(f"Total barcodes from file: {len(barcodes)}")
            return barcodes
        except Exception as e:
            logging.warning(f"Error reading barcode file {barcode_path}: {e}")
            return None
    else:
        return None

def save_classification_matrices(
    save_dir,
    union_variants,
    barcode_path=None,
    ref_classifications=None,
    alt_classifications=None,
    missing_classifications=None,
    unknown_classifications=None
):
    """Save classification results as sparse matrices in HDF5 format."""
    # Extract unique variants
    variants = sorted(set(
        f"{variant['CHROM']}:{variant['POS']}_{variant['REF']}->{variant['ALT']}" 
        for variant in union_variants
    ))
    var_to_idx = {variant: idx for idx, variant in enumerate(variants)}

    # Gather barcodes
    barcodes = read_barcodes_from_file(barcode_path)
    if barcodes is None:
        barcode_set = set()
        for classifications in [ref_classifications, alt_classifications, missing_classifications, unknown_classifications]:
            if classifications:
                for classification in classifications:
                    classification_dict = dict(classification)
                    barcode = classification_dict.get('cell_barcode')
                    if barcode:
                        barcode_set.add(barcode)
        barcodes = sorted(barcode_set)
        logging.info(f"# Total barcodes generated from data: {len(barcodes)}")
    bc_to_idx = {barcode: idx for idx, barcode in enumerate(barcodes)}

    # Initialize classification matrices
    classification_matrices = {
        'ref': defaultdict(int),
        'alt': defaultdict(int),
        'missing': defaultdict(int) if missing_classifications else None,
        'unknown': defaultdict(int) if unknown_classifications else None
    }

    # Function to populate classification counts
    def populate_counts(classifications, cls_type):
        if not classifications:
            return
        for classification in classifications:
            cls_dict = dict(classification)
            variant = cls_dict.get('variant')
            barcode = cls_dict.get('cell_barcode')
            if variant in var_to_idx and barcode in bc_to_idx:
                r = var_to_idx[variant]
                c = bc_to_idx[barcode]
                classification_matrices[cls_type][(r, c)] += 1

    # Populate counts for each classification
    populate_counts(ref_classifications, 'ref')
    populate_counts(alt_classifications, 'alt')
    if missing_classifications:
        populate_counts(missing_classifications, 'missing')
    if unknown_classifications:
        populate_counts(unknown_classifications, 'unknown')

    # Function to create and save CSR matrix
    def build_csr_and_save(cls_type):
        if cls_type not in classification_matrices or not classification_matrices[cls_type]:
            logging.info(f"No data for classification '{cls_type}'. Skipping.")
            return
        data = []
        row_indices = []
        col_indices = []
        for (r, c), count in classification_matrices[cls_type].items():
            row_indices.append(r)
            col_indices.append(c)
            data.append(count)
        matrix = csr_matrix((data, (row_indices, col_indices)), shape=(len(variants), len(barcodes)))
        out_h5 = os.path.join(save_dir, f"{cls_type}_matrix.h5")
        with h5py.File(out_h5, "w") as h5f:
            h5f.create_dataset("data", data=matrix.data)
            h5f.create_dataset("indices", data=matrix.indices)
            h5f.create_dataset("indptr", data=matrix.indptr)
            h5f.attrs['shape'] = matrix.shape
        logging.info(f"Saved {cls_type} classification matrix to {out_h5}")

    # Build and save matrices
    for cls_type in ['ref', 'alt', 'missing', 'unknown']:
        if cls_type == 'missing' and not missing_classifications:
            continue
        if cls_type == 'unknown' and not unknown_classifications:
            continue
        build_csr_and_save(cls_type)

    # Save variant and barcode mappings
    joblib.dump((variants, barcodes), os.path.join(save_dir, "variant_barcode_mappings.pkl"))
    logging.info("Saved variant and barcode mappings.")

############################################################################
# 6. Main Function
############################################################################
def main():
    args, vcf_files_with_origins = parse_arguments()
    setup_logging(args.save_dir)
    logging.info("Command Line: " + " ".join(sys.argv))
    logging.info("=== [scVarID] Program Started ===")

    compute_missing_unknown = not args.ref_alt_only

    # Step 1: Retrieve chromosome lengths from BAM
    step_start = log_step_start("Chromosome length retrieval")
    try:
        with pysam.AlignmentFile(args.bam_path, "rb") as bam:
            all_chroms = list(bam.references)
            all_lens = list(bam.lengths)
        chrom_len = dict(zip(all_chroms, all_lens))
        if args.chromosomes:
            filtered = {c: chrom_len[c] for c in args.chromosomes if c in chrom_len}
            chrom_len = filtered
        if not chrom_len:
            logging.error("No valid chromosomes found. Aborting.")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error retrieving chromosome lengths: {e}", exc_info=True)
        sys.exit(1)
    log_step_end("Chromosome length retrieval", step_start)

    window_size = args.window_size
    padding = args.padding
    num_cores = args.num_cores
    logging.info(f"[INFO] Parallel Processing => n_jobs={num_cores}, window_size={window_size}, padding={padding}")

    # Step 2: Generate chunk list
    step_start = log_step_start("Generate chunk list")
    all_chunks = []
    for chrom, length in chrom_len.items():
        st = 1
        while st <= length:
            e = min(st + window_size - 1, length)
            all_chunks.append((chrom, st, e))
            st = e + 1
    logging.info(f"Total chunk count: {len(all_chunks)}")
    log_step_end("Generate chunk list", step_start)

    # Step 3: Parallel Chunk Processing
    tmp_dir = os.path.join(args.save_dir, "tmp_chunks")
    os.makedirs(tmp_dir, exist_ok=True)
    step_start = log_step_start("Parallel chunk processing")
    chunk_files = Parallel(n_jobs=num_cores)(
        delayed(process_one_chunk_and_dump)(
            chrom=c[0],
            win_start=c[1],
            win_end=c[2],
            vcf_files_with_origins=vcf_files_with_origins,
            bam_path=args.bam_path,
            padding=padding,
            compute_missing_unknown=compute_missing_unknown,
            tmp_dir=tmp_dir
        )
        for c in all_chunks
    )
    log_step_end("Parallel chunk processing", step_start)

    # Step 4: Merge Chunk Results
    step_start = log_step_start("Merge chunk results")
    final_variants = []
    ref = []
    alt = []
    missing = []
    unknown = []
    for cf in chunk_files:
        if cf is None:
            continue
        try:
            chunk_data = joblib.load(cf)
            final_variants.extend(chunk_data.get("union_variants", []))
            ref.extend(chunk_data.get("ref", []))
            alt.extend(chunk_data.get("alt", []))
            if compute_missing_unknown:
                missing.extend(chunk_data.get("missing", []))
                unknown.extend(chunk_data.get("unknown", []))
        except Exception as e:
            logging.warning(f"Error loading chunk file {cf}: {e}")
    log_step_end("Merge chunk results", step_start)

    # Step 5: Save Final Classification Matrices
    step_start = log_step_start("Save classification results")
    save_classification_matrices(
        save_dir=args.save_dir,
        union_variants=final_variants,
        barcode_path=args.barcode_path,
        ref_classifications=ref,
        alt_classifications=alt,
        missing_classifications=missing if compute_missing_unknown else None,
        unknown_classifications=unknown if compute_missing_unknown else None
    )
    log_step_end("Save classification results", step_start)

    logging.info("=== All steps completed successfully ===")

if __name__ == "__main__":
    main()
