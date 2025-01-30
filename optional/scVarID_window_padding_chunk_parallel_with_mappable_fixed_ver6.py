#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
scVarID pipeline:
 - Window+padding => chunk-level parallel
 - PyRanges-based overlap => 'mappable reads', 'mappable variants'
 - Accurate Indel handling using original classification logic
 - Logs all relevant info (mappable, classification counts)
 - Final HDF5 output
"""

import os
import sys
import time
import logging
import argparse
import re
from collections import defaultdict

import pysam
import pyranges as pr
import h5py
import joblib
import pandas as pd
from scipy.sparse import csr_matrix
from joblib import Parallel, delayed

############################################################################
# 0. Logging helpers
############################################################################

def log_step_start(step_name):
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    end_ = time.time()
    elapsed = end_ - start_time
    hh = int(elapsed // 3600)
    mm = int((elapsed % 3600) // 60)
    ss = int(elapsed % 60)
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {hh:02d}h {mm:02d}m {ss:02d}s ===\n")

############################################################################
# 1. Parse args, setup logging
############################################################################

def parse_arguments():
    """
    Parses command-line arguments and returns them along with variant file origins.
    """
    parser = argparse.ArgumentParser(description='scVarID chunk parallel + mappable + accurate Indel handling.')

    # Add argument for variant files
    parser.add_argument('--variant-files', nargs='+', required=True, metavar=('ORIGIN', 'FILE_PATH'),
                        help='List of origin and variant file path pairs.')

    # Add arguments for BAM file, barcode file, and save directory
    parser.add_argument('--bam-path', required=True, help='Path to the BAM file.')
    parser.add_argument('--barcode-path', required=False, help='Path to the barcode file (optional).')
    parser.add_argument('--save-dir', required=True, help='Save directory path.')

    # Add argument for number of cores
    parser.add_argument('--num-cores', type=int, default=4, help='Number of cores to use for parallel processing.')

    # Add action
    parser.add_argument('--ref-alt-only', action='store_true',
                        help='Compute only ref and alt classifications, skipping missing and unknown.')

    # Add argument for chromosomes
    parser.add_argument('--chromosomes', nargs='+', required=False,
                        help='List of chromosomes to process (e.g., chr1 chr2). If not specified, process all.')

    # Add argument for window size
    parser.add_argument('--window-size', type=int, default=10_000_000,
                        help='Window size (bp) for chunk processing. Default=10,000,000bp')

    # Add argument for padding
    parser.add_argument('--padding', type=int, default=500,
                        help='Padding around chunk boundary (bp)')

    # Add argument for barcode tag
    parser.add_argument('--barcode-tag', type=str, default='CB', help='BAM tag for cell barcode (default: CB)')

    # Parse arguments
    args = parser.parse_args()

    # Generate list of variant files
    vcf_files_with_origins = []
    if len(args.variant_files) % 2 != 0:
        parser.error('Each origin must be paired with a file path in --variant-files.')

    for i in range(0, len(args.variant_files), 2):
        origin = args.variant_files[i]
        file_path = args.variant_files[i + 1]
        vcf_files_with_origins.append((file_path, origin))

    return args, vcf_files_with_origins

def setup_logging(save_dir):
    """
    Initializes the logging settings.
    Sets up both log file and console output.
    """
    os.makedirs(save_dir, exist_ok=True)
    log_file = os.path.join(save_dir, 'processing.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

############################################################################
# 2. VCF region fetch + filter
############################################################################

def is_vcf_file(fpath):
    return fpath.endswith('.vcf') or fpath.endswith('.vcf.gz') or fpath.endswith('.bcf')

def extract_variants_vcf_region(vcf_file, origin, chrom, start, end):
    vs = []
    try:
        with pysam.VariantFile(vcf_file) as vf:
            for rec in vf.fetch(chrom, start, end):
                if "PASS" in rec.filter:
                    c_ = rec.chrom
                    p_ = rec.pos
                    ref_ = rec.ref
                    alts_ = ",".join(rec.alts) if rec.alts else ""
                    var = {
                        "CHROM": c_,
                        "POS": p_,
                        "REF": ref_,
                        "ALT": alts_,
                        "QUAL": rec.qual if rec.qual else ".",
                        "FILTER": ";".join(rec.filter) if rec.filter else ".",
                        "ORIGIN": origin
                    }
                    vs.append(var)
    except Exception as e:
        logging.warning(f"extract_variants_vcf_region error in {vcf_file}: {e}")
    logging.info(f"Extracted {len(vs)} variants from region {chrom}:{start}-{end} in {vcf_file}")
    return vs

def extract_variants_txt_region(txt_file, origin, chrom, start, end):
    vs = []
    try:
        df = pd.read_csv(txt_file, sep='\t', header=None,
                         names=["CHROM","POS","REF","ALT"],
                         dtype={"CHROM":str,"POS":int,"REF":str,"ALT":str})
        if df.shape[1] != 4:
            raise ValueError(f"TXT file {txt_file} must have 4 columns.")
        sub = df[(df["CHROM"] == chrom) & (df["POS"] >= start) & (df["POS"] <= end)]
        for _, row in sub.iterrows():
            var = {
                "CHROM": row["CHROM"],
                "POS": row["POS"],
                "REF": row["REF"],
                "ALT": row["ALT"],
                "QUAL": ".",
                "FILTER": ".",
                "ORIGIN": origin
            }
            vs.append(var)
    except Exception as e:
        logging.warning(f"extract_variants_txt_region error in {txt_file}: {e}")
    logging.info(f"Extracted {len(vs)} variants from region {chrom}:{start}-{end} in {txt_file}")
    return vs

def filter_variants(variants):
    union_map = {}
    for v_ in variants:
        ref = v_["REF"]
        alt = v_["ALT"]
        key = (v_["CHROM"], v_["POS"], ref, alt)

        # Exclude multi-allelic references or alternates
        if "," in ref or "," in alt:
            continue
        if (len(ref) > 1 or len(alt) > 1) and ref[0] != alt[0]:
            continue

        if key not in union_map:
            union_map[key] = v_
        else:
            # Merge ORIGIN
            exist_o = union_map[key]["ORIGIN"].split(", ")
            new_o = v_["ORIGIN"]
            combo = exist_o + [new_o]
            # Keep order normal, tumor, SCREADCOUNT
            origin_ordered = []
            for tag_ in ["normal","tumor","SCREADCOUNT"]:
                if tag_ in combo and tag_ not in origin_ordered:
                    origin_ordered.append(tag_)
            union_map[key]["ORIGIN"] = ", ".join(origin_ordered)

    unioned = list(union_map.values())
    for i, v_ in enumerate(unioned, 1):
        v_["VINDEX"] = f"vi_{i}"
    logging.info(f"Filtered variants count: {len(unioned)}")
    return unioned

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    all_vs = []
    for (fpath, origin) in vcf_files_with_origins:
        if is_vcf_file(fpath):
            vs = extract_variants_vcf_region(fpath, origin, chrom, start, end)
        elif fpath.endswith('.txt'):
            vs = extract_variants_txt_region(fpath, origin, chrom, start, end)
        else:
            vs = []
            logging.warning(f"Skipping unsupported file: {fpath}")
        all_vs.extend(vs)
    union_ = filter_variants(all_vs)
    return union_

############################################################################
# 3. Read fetch + PyRanges overlap => mappable
############################################################################

def parse_cigar(cigar_string):
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)]

def calculate_end_position_cigar(read):
    ref_pos = read.reference_start
    ops = parse_cigar(read.cigarstring)
    for length, op_ in ops:
        if op_ in "M=XDN":
            ref_pos += length
    return ref_pos

def extract_reads_info_region(bam_path, chrom, start, end, barcode_tag):
    data = []
    read_name_counts = defaultdict(int)
    total_reads = 0
    total_barcodes = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start-1, end):
            try:
                bc = read.get_tag(barcode_tag)
                if not bc:
                    continue
            except KeyError:
                continue
            rname = read.query_name
            c_ = read_name_counts[rname]
            read_unique = f"{rname}_{c_}"
            read_name_counts[rname] += 1

            read_end = calculate_end_position_cigar(read)
            data.append([
                rname, read_unique, bc, chrom,
                read.reference_start, read_end
            ])
            total_reads += 1
            total_barcodes += 1
    df = pd.DataFrame(data, columns=[
        "Read_Name","Read_Unique_Name","Cell_Barcode","Chromosome",
        "Start_Position","End_Position"
    ])
    logging.info(f"Extracted {len(df)} reads from region {chrom}:{start}-{end} in {bam_path}")
    logging.info(f"Total reads with barcode: {total_barcodes}")
    return df

def annotate_reads_with_variants_pyranges(df_reads, union_variants):
    if not len(union_variants):
        df_reads["VINDEX"] = ""
        return df_reads

    df_var = pd.DataFrame(union_variants).copy()
    def calc_end(row):
        ref_len = len(row["REF"])
        alt_len = len(row["ALT"])
        st = row["POS"] - 1
        if ref_len > alt_len:
            return st + ref_len
        elif ref_len < alt_len:
            return st + alt_len
        else:
            return st + ref_len

    df_var["Start"] = df_var["POS"] - 1
    df_var["End"] = df_var.apply(calc_end, axis=1)
    # Rename CHROM -> Chromosome for PyRanges
    df_var.rename(columns={"CHROM": "Chromosome"}, inplace=True)

    var_pr = pr.PyRanges(df_var[["Chromosome","Start","End","VINDEX"]])

    tmp_rd = df_reads[["Read_Unique_Name","Chromosome","Start_Position","End_Position"]].copy()
    tmp_rd.rename(columns={"Start_Position":"Start","End_Position":"End"}, inplace=True)
    rd_pr = pr.PyRanges(tmp_rd)

    overlap = rd_pr.join(var_pr)
    odf = overlap.df
    # Group by read => all VINDEX
    read_vindex = odf.groupby("Read_Unique_Name")["VINDEX"].apply(lambda x: ",".join(x)).reset_index()

    df_reads = df_reads.merge(read_vindex, on="Read_Unique_Name", how="left")
    df_reads["VINDEX"] = df_reads["VINDEX"].fillna("")

    total_ = len(df_reads)
    mappable_ = (df_reads["VINDEX"] != "").sum()
    perc = (100 * mappable_ / total_) if total_ > 0 else 0
    logging.info(f"# mappable reads: {mappable_} ({perc:.2f}%)")

    return df_reads

def annotate_variants_with_reads_pyranges(union_variants, df_reads):
    if not union_variants:
        return union_variants

    df_var = pd.DataFrame(union_variants).copy()
    def calc_end(row):
        ref_len = len(row["REF"])
        alt_len = len(row["ALT"])
        st = row["POS"] - 1
        if ref_len > alt_len:
            return st + ref_len
        elif ref_len < alt_len:
            return st + alt_len
        else:
            return st + ref_len

    df_var["Start"] = df_var["POS"] - 1
    df_var["End"] = df_var.apply(calc_end, axis=1)
    df_var.rename(columns={"CHROM": "Chromosome"}, inplace=True)

    var_pr = pr.PyRanges(df_var[["Chromosome","Start","End","VINDEX"]])

    tmp_rd = df_reads[["Chromosome","Start_Position","End_Position","Read_Unique_Name"]].copy()
    tmp_rd.rename(columns={"Start_Position":"Start","End_Position":"End"}, inplace=True)
    rd_pr = pr.PyRanges(tmp_rd)

    overlap = var_pr.join(rd_pr)
    odf = overlap.df
    # Group by variant => read list
    var_reads = odf.groupby("VINDEX")["Read_Unique_Name"].apply(lambda x: ",".join(x)).reset_index()

    df_var = df_var.merge(var_reads, on="VINDEX", how="left")
    df_var["Read_Names"] = df_var["Read_Unique_Name"].fillna("")
    df_var.drop(columns=["Start","End","Read_Unique_Name"], inplace=True)

    # Restore 'CHROM' from 'Chromosome'
    df_var.rename(columns={"Chromosome": "CHROM"}, inplace=True)

    updated_vars = df_var.to_dict(orient="records")

    total_var = len(updated_vars)
    mappable_var = sum(1 for vv in updated_vars if vv["Read_Names"] != "")
    p_ = (100 * mappable_var / total_var) if total_var > 0 else 0
    logging.info(f"# mappable variants: {mappable_var} ({p_:.2f}%)")

    return updated_vars

def process_reads_and_variants_pyranges(df_reads, union_variants):
    df_reads2 = annotate_reads_with_variants_pyranges(df_reads, union_variants)
    updated_vars = annotate_variants_with_reads_pyranges(union_variants, df_reads2)
    return df_reads2, updated_vars

############################################################################
# 4. Classification handle_insertion / handle_deletion
############################################################################

def prime_read(read, barcode_tag):
    read_name = read.query_name
    try:
        cell_barcode = read.get_tag(barcode_tag)
    except KeyError:
        cell_barcode = None

    sequences, positions, operations = [], [], []
    current_position = read.reference_start + 1
    current_read_pos = 0

    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', read.cigarstring)

    for length, op in cigar_ops:
        length = int(length)
        if op in "M=X":
            for _ in range(length):
                sequences.append(read.query_sequence[current_read_pos])
                positions.append(current_position)
                operations.append(op)
                current_position += 1
                current_read_pos += 1
        elif op == "I":
            for _ in range(length):
                sequences.append(read.query_sequence[current_read_pos])
                positions.append(-1)  # Use -1 as the standard representation for insertion positions
                operations.append(op)
                current_read_pos += 1
        elif op in "DNP":
            for _ in range(length):
                sequences.append('-')
                positions.append(current_position)
                operations.append(op)
                current_position += 1
        elif op == "S":
            current_read_pos += length
        elif op == "H":
            pass

    return read_name, cell_barcode, sequences, positions, operations

def handle_snv(sequences, positions, pos, ref, alt, read_name, cell_barcode):
    for p, seq in zip(positions, sequences):
        if p == pos:
            if seq == ref:
                return 'ref', cell_barcode, pos, read_name
            elif seq == alt:
                return 'alt', cell_barcode, pos, read_name
            elif seq == '-':
                return 'missing', cell_barcode, pos, read_name
            else:
                return 'unknown', cell_barcode, pos, read_name
    return 'unknown', cell_barcode, pos, read_name

def handle_insertion(sequences, positions, operations, pos_, ref_, alt_, read_name, cell_barcode):
    try:
        pos_index = positions.index(pos_)

        # Check for missing status first
        if sequences[pos_index] == '-':
            return 'missing', cell_barcode, pos_, read_name

        # Check for ref status
        # All operations in the reference region should be '='
        if all(op == '=' for op in operations[pos_index : pos_index + len(ref_)]) and \
           (pos_index + len(ref_) < len(operations) and operations[pos_index + len(ref_)] == '='):
            return 'ref', cell_barcode, pos_, read_name

        # Check for alt status
        insert_length = len(alt_) - len(ref_)
        if insert_length <= 0:
            return 'unknown', cell_barcode, pos_, read_name

        inserted_bases = ''.join(sequences[pos_index + len(ref_) : pos_index + len(ref_) + insert_length])
        expected_insert = alt_[len(ref_):]

        # Using equality: check if inserted_bases match expected_insert
        if inserted_bases == expected_insert:
            # Also, ensure that the operations are 'I'
            if all(op == 'I' for op in operations[pos_index + len(ref_) : pos_index + len(ref_) + insert_length]):
                return 'alt', cell_barcode, pos_, read_name
            else:
                return 'missing', cell_barcode, pos_, read_name
        else:
            return 'missing', cell_barcode, pos_, read_name

    except (ValueError, IndexError):
        return 'unknown', cell_barcode, pos_, read_name

    return 'unknown', cell_barcode, pos_, read_name

def handle_deletion(sequences, positions, operations, pos_, ref_, alt_, read_name, cell_barcode):
    try:
        pos_index = positions.index(pos_)

        # Check for missing status first
        if sequences[pos_index] == '-':
            return 'missing', cell_barcode, pos_, read_name

        # Check for ref status
        # All operations in the reference region should be '='
        if all(op == '=' for op in operations[pos_index : pos_index + len(ref_)]) and \
           (pos_index + len(ref_) < len(operations) and operations[pos_index + len(ref_)] == '='):
            return 'ref', cell_barcode, pos_, read_name

        # Check for alt status
        del_length = len(ref_) - len(alt_)
        if del_length <= 0:
            return 'unknown', cell_barcode, pos_, read_name

        deleted_bases = ''.join(sequences[pos_index + len(alt_) : pos_index + len(alt_) + del_length])
        expected_delete = ref_[:del_length]

        # Using equality: check if deleted_bases match expected_delete
        if deleted_bases == expected_delete:
            # Also, ensure that the operations are 'D'
            if all(op == 'D' for op in operations[pos_index + len(alt_) : pos_index + len(alt_) + del_length]):
                return 'alt', cell_barcode, pos_, read_name
            else:
                return 'missing', cell_barcode, pos_, read_name
        else:
            return 'missing', cell_barcode, pos_, read_name

    except (ValueError, IndexError):
        return 'unknown', cell_barcode, pos_, read_name

    return 'unknown', cell_barcode, pos_, read_name

def classify_variant(seqs, poss, ops, pos_, ref_, alt_, variant_type, op_, rname, bc):
    if variant_type == 'deletion':
        return handle_deletion(seqs, poss, ops, pos_, ref_, alt_, rname, bc)
    elif variant_type == 'insertion':
        return handle_insertion(seqs, poss, ops, pos_, ref_, alt_, rname, bc)
    else:
        return handle_snv(seqs, poss, pos_, ref_, alt_, rname, bc)

############################################################################
# 5. Classification Functions
############################################################################

def create_variants_dict(selected_variants):
    variants_dict = {}
    for variant in selected_variants:
        chrom = variant['CHROM']
        pos = int(variant['POS'])
        ref = variant['REF']
        alt_list = variant['ALT'].split(',')  # Handle cases where ALT values are multiple

        for alt in alt_list:
            if len(ref) > len(alt):
                variant_type = 'deletion'
            elif len(ref) < len(alt):
                variant_type = 'insertion'
            else:
                variant_type = 'snv'

            key = (chrom, pos)
            if key not in variants_dict:
                variants_dict[key] = []
            variants_dict[key].append({'variant': (chrom, pos, ref, alt), 'type': variant_type})
    return variants_dict

def process_read(read, variants_dict, compute_missing_unknown=True, barcode_tag='CB'):
    ref_classifications_local = set()
    alt_classifications_local = set()
    missing_classifications_local = set()
    unknown_classifications_local = set()

    read_name, cell_barcode, sequences, positions, operations = prime_read(read, barcode_tag)
    chrom = read.reference_name

    read_start = read.reference_start + 1
    read_end = read.reference_end

    processed_positions = set()

    for pos in positions:
        key = (chrom, pos)

        if key in variants_dict and key not in processed_positions:
            processed_positions.add(key)

            for variant_data in variants_dict[key]:
                variant = variant_data['variant']
                variant_type = variant_data['type']

                ref_len = len(variant[2])
                alt_len = len(variant[3])

                if pos < read_start or pos > read_end - (alt_len - 1) or pos > read_end - (ref_len - 1):
                    continue

                classification, barcode, variant_id, read_id = classify_variant(
                    sequences, positions, operations, pos, variant[2], variant[3], variant_type,
                    operations[positions.index(pos)], read_name, cell_barcode
                )

                classification_data = {
                    'read_name': read_id,
                    'cell_barcode': barcode,
                    'variant': f"{variant[0]}:{variant[1]}_{variant[2]}->{variant[3]}",
                }

                if classification == 'ref':
                    ref_classifications_local.add(frozenset(classification_data.items()))
                elif classification == 'alt':
                    alt_classifications_local.add(frozenset(classification_data.items()))
                elif classification == 'missing' and compute_missing_unknown:
                    missing_classifications_local.add(frozenset(classification_data.items()))
                elif classification == 'unknown' and compute_missing_unknown:
                    unknown_classifications_local.add(frozenset(classification_data.items()))

    return (
        ref_classifications_local,
        alt_classifications_local,
        missing_classifications_local,
        unknown_classifications_local,
    )

def process_read_chunk(
    bam_file_path, 
    read_unique_names_set, 
    variants_dict, 
    read_mapping, 
    compute_missing_unknown=True,
    barcode_tag='CB'
):
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
                    read, variants_dict, compute_missing_unknown, barcode_tag
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
    compute_missing_unknown=True,
    barcode_tag='CB'
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
    results = Parallel(n_jobs=num_cores, backend='loky')(
        delayed(process_read_chunk)(
            bam_path, 
            set(chunk), 
            variants_dict, 
            read_mapping, 
            compute_missing_unknown,
            barcode_tag
        ) for chunk in read_unique_names_chunks
    )

    for result in results:
        ref_classifications.extend(result[0])
        alt_classifications.extend(result[1])
        missing_classifications.extend(result[2])
        unknown_classifications.extend(result[3])

    return ref_classifications, alt_classifications, missing_classifications, unknown_classifications

############################################################################
# 6. Save final matrix
############################################################################

def read_barcodes_from_file(barcode_path):
    """
    Reads barcodes from a file if barcode_path is provided, else returns None.
    Assumes one barcode per line.
    """
    if barcode_path:
        with open(barcode_path, 'r') as f:
            barcodes = [line.strip() for line in f if line.strip()]
        logging.info(f"Total barcodes from file: {len(barcodes)}")
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
    logging.info(f"=== Save Classification Matrices ===")
    
    # Extract variant list
    variants = [f"{variant['CHROM']}:{variant['POS']}_{variant['REF']}->{variant['ALT']}" for variant in union_variants]
    variants = list(set(variants))
    variants.sort()
    var_to_idx = {v:i for i,v in enumerate(variants)}
    logging.info(f"Total unique variants => {len(variants)}")

    # Read barcodes from file or generate from data
    if barcode_path:
        barcodes = read_barcodes_from_file(barcode_path)
    else:
        # Generate barcodes from classification data
        all_barcodes = set()
        def gather_barcodes(cls_):
            if cls_:
                for c_ in cls_:
                    classification_dict = dict(c_)
                    bc_ = classification_dict.get("cell_barcode", "")
                    if bc_:
                        all_barcodes.add(bc_)
        
        gather_barcodes(ref_classifications)
        gather_barcodes(alt_classifications)
        if missing_classifications:
            gather_barcodes(missing_classifications)
        if unknown_classifications:
            gather_barcodes(unknown_classifications)
        
        if all_barcodes:
            barcodes = sorted(all_barcodes)
            logging.info(f"Total unique barcodes generated from data: {len(barcodes)}")
        else:
            # Assign barcodes based on sorted order if no barcodes are found
            barcodes = []
            logging.info("No barcodes found in data. Assigning default barcodes.")
        
    if not barcodes:
        # If no barcodes are found or provided, assign default barcodes
        barcodes = [f"barcode_{i}" for i in range(1, len(variants)+1)]
        logging.info(f"Assigned {len(barcodes)} default barcodes.")

    bc_to_idx = {b:i for i,b in enumerate(barcodes)}
    logging.info(f"Total unique barcodes => {len(barcodes)}")

    from collections import defaultdict

    def build_csr_and_save(cls_data, name):
        if not cls_data:
            logging.info(f"No data for {name} matrix.")
            return
        dmap = defaultdict(int)
        for c_ in cls_data:
            classification_dict = dict(c_)
            var_ = classification_dict["variant"]
            bc_ = classification_dict["cell_barcode"]
            if var_ not in var_to_idx or bc_ not in bc_to_idx:
                continue
            r = var_to_idx[var_]
            co = bc_to_idx[bc_]
            dmap[(r, co)] += 1

        data = []
        row_inds = []
        col_inds = []
        for (r, co), val in dmap.items():
            row_inds.append(r)
            col_inds.append(co)
            data.append(val)

        mat = csr_matrix((data, (row_inds, col_inds)), shape=(len(variants), len(barcodes)))
        out_h5 = os.path.join(save_dir, f"{name}_matrix.h5")
        with h5py.File(out_h5, "w") as hf:
            hf.create_dataset("data", data=mat.data)
            hf.create_dataset("indices", data=mat.indices)
            hf.create_dataset("indptr", data=mat.indptr)
            hf.attrs['shape'] = mat.shape
        logging.info(f"Saved => {out_h5} | shape={mat.shape} | nnz={mat.nnz}")

    build_csr_and_save(ref_classifications, "ref")
    build_csr_and_save(alt_classifications, "alt")
    if missing_classifications is not None:
        build_csr_and_save(missing_classifications, "missing")
    if unknown_classifications is not None:
        build_csr_and_save(unknown_classifications, "unknown")

    # Save variant and barcode mappings
    joblib.dump((variants, barcodes), os.path.join(save_dir, "variant_barcode_mappings.pkl"))
    logging.info(f"Saved => variant_barcode_mappings.pkl")

############################################################################
# 7. Helper Functions
############################################################################

def split_list(lst, num_sublists):
    """Splits a list into num_sublists approximately equal parts."""
    if num_sublists <= 0:
        raise ValueError("num_sublists must be positive")
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
# 8. main => chunk-level parallel
############################################################################

def process_one_chunk_and_dump(
    chrom,
    win_start,
    win_end,
    vcf_files_with_origins,
    bam_path,
    padding=500,
    compute_missing_unknown=True,
    tmp_dir=None,
    barcode_tag='CB'
):
    chunk_str = f"{chrom}:{win_start}-{win_end}"
    logging.info(f"[CHUNK] Start {chunk_str}")

    # 1) VCF => union_variants
    union_variants = process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(union_variants)} variants")
    if not union_variants:
        return None

    # 2) Reads => region-based
    df_reads = extract_reads_info_region(bam_path, chrom, win_start, win_end, barcode_tag)
    if df_reads.empty:
        logging.info(f"[CHUNK] {chunk_str} no reads with barcode => skip.")
        return None
    logging.info(f"[CHUNK] {chunk_str} extracted {len(df_reads)} reads")

    # 3) PyRanges overlap => mappable reads => mappable variants
    df_reads2, updated_union_variants = process_reads_and_variants_pyranges(df_reads, union_variants)

    # 4) classification => accurate Indel handling => store
    read_mapping = df_reads2.set_index("Read_Name")["Read_Unique_Name"].to_dict()
    refC, altC, missC, unkC = process_bam_data_parallel(
        bam_path,
        [r["Read_Unique_Name"] for r in df_reads2.to_dict(orient="records")],
        updated_union_variants,
        read_mapping,
        num_cores=1,  # Already parallelized at higher level
        compute_missing_unknown=compute_missing_unknown,
        barcode_tag=barcode_tag
    )

    # Log classification counts
    logging.info(f"[CHUNK] {chunk_str} classification => ref:{len(refC)}, alt:{len(altC)}, missing:{len(missC)}, unknown:{len(unkC)}")

    chunk_data = {
        "union_variants": updated_union_variants,
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

def main():
    args, vcf_files_with_origins = parse_arguments()
    setup_logging(args.save_dir)
    log_command_line()
    logging.info("=== [scVarID] Program Started ===")

    compute_missing_unknown = not args.ref_alt_only
    num_cores = check_num_cores(args.num_cores)
    wsize = args.window_size
    pad  = args.padding
    barcode_tag = args.barcode_tag

    # 1) get chrom lengths
    step_t = log_step_start("Chromosome length retrieval")
    chrom_len = get_chromosome_lengths(args.bam_path, target_chromosomes=args.chromosomes)
    if not chrom_len:
        logging.error("No valid chromosome => abort.")
        sys.exit(1)
    log_step_end("Chromosome length retrieval", step_t)

    logging.info(f"[INFO] chunk-level parallel => n_jobs={num_cores}, window_size={wsize}, padding={pad}")

    # 2) build chunk list
    step_t = log_step_start("Generate chunk list")
    all_chunks = []
    for c_, length_ in chrom_len.items():
        st = 1
        while st <= length_:
            e_ = min(st + wsize - 1, length_)
            all_chunks.append((c_, st, e_))
            st = e_ + 1
    logging.info(f"Total chunk count: {len(all_chunks)}")
    log_step_end("Generate chunk list", step_t)

    tmp_dir = os.path.join(args.save_dir, "tmp_chunks")
    os.makedirs(tmp_dir, exist_ok=True)

    # 3) parallel chunk
    step_t = log_step_start("Parallel chunk processing")
    chunk_files = Parallel(n_jobs=num_cores)(
        delayed(process_one_chunk_and_dump)(
            chrom=c[0],
            win_start=c[1],
            win_end=c[2],
            vcf_files_with_origins=vcf_files_with_origins,
            bam_path=args.bam_path,
            padding=pad,
            compute_missing_unknown=compute_missing_unknown,
            tmp_dir=tmp_dir,
            barcode_tag=barcode_tag
        )
        for c in all_chunks
    )
    log_step_end("Parallel chunk processing", step_t)

    # 4) merge
    step_t = log_step_start("Merge chunk results")
    final_variants = []
    ref_ = []
    alt_ = []
    missing_ = []
    unknown_ = []
    for cf_ in chunk_files:
        if cf_ is None:
            continue
        cdata = joblib.load(cf_)
        final_variants.extend(cdata["union_variants"])
        ref_.extend(cdata["ref"])
        alt_.extend(cdata["alt"])
        if compute_missing_unknown:
            missing_.extend(cdata["missing"])
            unknown_.extend(cdata["unknown"])
        del cdata
    log_step_end("Merge chunk results", step_t)

    # 5) save
    step_t = log_step_start("Save classification results")
    save_classification_matrices(
        save_dir=args.save_dir,
        union_variants=final_variants,
        barcode_path=args.barcode_path,
        ref_classifications=ref_,
        alt_classifications=alt_,
        missing_classifications=missing_ if compute_missing_unknown else None,
        unknown_classifications=unknown_ if compute_missing_unknown else None
    )
    log_step_end("Save classification results", step_t)

    logging.info("=== All steps completed successfully ===")

############################################################################
# 9. Additional Helper Functions
############################################################################

def log_command_line():
    """
    Logs the command-line arguments used to run the script.
    """
    command_line = " ".join(sys.argv)  # Get the full command-line arguments
    logging.info(f"Command Line: {command_line}")

def get_chromosome_lengths(bam_path, target_chromosomes=None):
    """
    Extract chromosome lengths from BAM file.
    If target_chromosomes is provided, filter chromosomes accordingly.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        all_chroms = list(bam_file.references)
        all_lengths = list(bam_file.lengths)
    chrom_len = dict(zip(all_chroms, all_lengths))

    # Filter if needed
    if target_chromosomes:
        filtered = {}
        for c_ in target_chromosomes:
            if c_ in chrom_len:
                filtered[c_] = chrom_len[c_]
        return filtered
    return chrom_len

############################################################################
# 10. Run main
############################################################################

if __name__ == "__main__":
    main()
