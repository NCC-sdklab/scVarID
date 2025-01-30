#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse
import pysam
import h5py
import joblib
import pandas as pd
from collections import defaultdict
from scipy.sparse import csr_matrix

from joblib import Parallel, delayed   # ★★★ 중요: chunk 병렬화

#####################################
# 1) 유틸 함수들
#####################################

def parse_arguments():
    """
    --window-size, --padding, --num-cores, --chromosomes 등
    """
    parser = argparse.ArgumentParser(description='scVarID-like window+padding + parallel chunk example.')

    # 기존 인자들
    parser.add_argument('--variant-files', nargs='+', required=True, metavar=('ORIGIN', 'FILE_PATH'),
                        help='List of origin and variant file path pairs.')
    parser.add_argument('--bam-path', required=True, help='Path to the BAM file.')
    parser.add_argument('--barcode-path', required=False, help='Path to the barcode file (optional).')
    parser.add_argument('--save-dir', required=True, help='Save directory path.')
    parser.add_argument('--num-cores', type=int, default=4, help='Number of cores for parallel processing.')
    parser.add_argument('--ref-alt-only', action='store_true',
                        help='Compute only ref and alt classifications, skipping missing and unknown.')
    parser.add_argument('--chromosomes', nargs='+', required=False,
                        help='List of chromosomes. If not specified, process all from BAM header.')

    # 추가: window-size, padding
    parser.add_argument('--window-size', type=int, default=10_000_000,
                        help='Window size (bp) for chunk processing. Default=10,000,000.')
    parser.add_argument('--padding', type=int, default=5000,
                        help='Padding size (bp) to fetch more reads around chunk boundary. Default=5000.')

    args = parser.parse_args()

    # scVarID처럼 (origin, filepath) 형태로 묶기
    vcf_files_with_origins = []
    if len(args.variant_files) % 2 != 0:
        parser.error('Each origin must be paired with a file path in --variant-files.')

    for i in range(0, len(args.variant_files), 2):
        origin = args.variant_files[i]
        file_path = args.variant_files[i+1]
        vcf_files_with_origins.append((file_path, origin))

    return args, vcf_files_with_origins


def setup_logging(save_dir):
    """
    로그 설정
    """
    os.makedirs(save_dir, exist_ok=True)
    log_file = os.path.join(save_dir, 'processing.log')

    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s'
    )
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


def check_num_cores(num_cores):
    import multiprocessing
    max_cores = multiprocessing.cpu_count()
    if num_cores > max_cores:
        logging.info(f"Warning: requested {num_cores} cores > available {max_cores}. Using {max_cores} instead.")
        return max_cores
    return num_cores


def get_chromosome_lengths(bam_path, target_chromosomes=None):
    """
    BAM 헤더에서 (chrom -> length) 추출
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        all_chroms = list(bam.references)
        all_lengths = list(bam.lengths)
    chrom_lengths = dict(zip(all_chroms, all_lengths))

    if target_chromosomes:
        filtered = {}
        for c in target_chromosomes:
            if c in chrom_lengths:
                filtered[c] = chrom_lengths[c]
        return filtered
    else:
        return chrom_lengths


def get_windows_for_chromosome(chrom, chrom_length, window_size):
    """
    chr 길이를 window_size씩 나누어 (start, end) yield
    """
    start = 1
    while start <= chrom_length:
        end = min(start + window_size - 1, chrom_length)
        yield (start, end)
        start = end + 1

#####################################
# 2) Variant 추출 (region 기반)
#####################################

def is_vcf_file(file_path):
    return (file_path.endswith('.vcf') or file_path.endswith('.vcf.gz') or file_path.endswith('.bcf'))

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    all_variants = []
    for file_path, origin in vcf_files_with_origins:
        if is_vcf_file(file_path):
            vs = extract_variants_vcf_region(file_path, origin, chrom, start, end)
        elif file_path.endswith('.txt'):
            vs = extract_variants_txt_region(file_path, origin, chrom, start, end)
        else:
            logging.warning(f"Unsupported file format: {file_path}")
            continue
        all_variants.extend(vs)
    union_variants = filter_variants(all_variants)
    return union_variants

def extract_variants_vcf_region(vcf_file, origin, chrom, start, end):
    formatted_variants = []
    try:
        with pysam.VariantFile(vcf_file) as vf:
            for record in vf.fetch(chrom, start, end):
                if "PASS" in record.filter:
                    v = {
                        "CHROM": record.chrom,
                        "POS": record.pos,
                        "REF": record.ref,
                        "ALT": ",".join(record.alts) if record.alts else "",
                        "QUAL": record.qual if record.qual else ".",
                        "FILTER": ";".join(record.filter),
                        "ORIGIN": origin
                    }
                    formatted_variants.append(v)
    except Exception as e:
        logging.warning(f"Error reading {vcf_file}: {e}")

    logging.info(f"Extracted {len(formatted_variants)} variants in region {chrom}:{start}-{end} from {vcf_file}")
    return formatted_variants

def extract_variants_txt_region(txt_file, origin, chrom, start, end):
    formatted_variants = []
    try:
        df = pd.read_csv(txt_file, sep='\t', header=None, names=["CHROM", "POS", "REF", "ALT"])
        region_df = df[(df["CHROM"] == chrom) & (df["POS"] >= start) & (df["POS"] <= end)]
        for _, row in region_df.iterrows():
            v = {
                "CHROM": row["CHROM"],
                "POS": row["POS"],
                "REF": row["REF"],
                "ALT": row["ALT"],
                "QUAL": ".",
                "FILTER": ".",
                "ORIGIN": origin
            }
            formatted_variants.append(v)
    except Exception as e:
        logging.warning(f"Error reading {txt_file}: {e}")

    logging.info(f"Extracted {len(formatted_variants)} variants in region {chrom}:{start}-{end} from {txt_file}")
    return formatted_variants

def filter_variants(variants):
    union_dict = {}
    for v in variants:
        key = (v["CHROM"], v["POS"], v["REF"], v["ALT"])
        if key not in union_dict:
            union_dict[key] = v
        else:
            old_origins = union_dict[key]["ORIGIN"].split(",")
            new_origins = v["ORIGIN"].split(",")
            union_dict[key]["ORIGIN"] = ",".join(set(old_origins + new_origins))
    union_list = list(union_dict.values())
    for i, variant in enumerate(union_list, 1):
        variant["VINDEX"] = f"vi_{i}"
    logging.info(f"Filtered variants count: {len(union_list)}")
    return union_list

#####################################
# 3) Reads 추출 (패딩 적용)
#####################################

def extract_reads_info_region_with_padding(bam_path, chrom, start, end, padding=3000):
    data = []
    read_name_counts = {}
    fetch_start = max(1, start - padding)
    fetch_end = end + padding

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # 0-based exclusive
        for read in bam.fetch(chrom, fetch_start-1, fetch_end):
            try:
                cb = read.get_tag("CB")
            except KeyError:
                continue

            if read.query_name not in read_name_counts:
                read_name_counts[read.query_name] = 0
            else:
                read_name_counts[read.query_name] += 1

            unique_id = read_name_counts[read.query_name]
            read_unique_name = f"{read.query_name}_{unique_id}"

            data.append([
                read.query_name,
                read_unique_name,
                cb,
                chrom,
                read.reference_start,  # 0-based
                read.reference_end     # 0-based
            ])

    df = pd.DataFrame(data, columns=[
        "Read_Name", "Read_Unique_Name", "Cell_Barcode", "Chromosome",
        "Start_Position", "End_Position"
    ])
    logging.info(f"[WithPadding] {chrom}:{start}-{end} (padding={padding}) => {len(df)} reads")
    return df

#####################################
# 4) 간단한 Overlap / Classification (데모용)
#####################################

def create_read_mapping(df_reads):
    g = df_reads.groupby("Read_Name")["Read_Unique_Name"].apply(list)
    return g.to_dict()

def process_reads_and_variants(df_reads, union_variants, chunk_start, chunk_end):
    selected_read_unique_names = df_reads["Read_Unique_Name"].tolist()
    final_variants = []
    for v in union_variants:
        pos = v["POS"]
        if pos < chunk_start or pos > chunk_end:
            continue
        final_variants.append(v)
    return df_reads, selected_read_unique_names, final_variants

def classify_in_parallel(bam_path, selected_read_unique_names, selected_variants, read_mapping,
                         compute_missing_unknown=True):
    """
    예시로, random classification을 만듦. 실제 scVarID로 치환 가능.
    """
    import random
    ref_classifications = []
    alt_classifications = []
    missing_classifications = []
    unknown_classifications = []

    for rn in selected_read_unique_names:
        for var in selected_variants:
            classification = random.choice(["ref", "alt", "missing", "unknown"])
            d = {
                "read_name": rn,
                "cell_barcode": "BC_TEST",
                "variant": f"{var['CHROM']}:{var['POS']}_{var['REF']}->{var['ALT']}"
            }
            if classification == "ref":
                ref_classifications.append(d)
            elif classification == "alt":
                alt_classifications.append(d)
            elif classification == "missing":
                if compute_missing_unknown:
                    missing_classifications.append(d)
            else:
                if compute_missing_unknown:
                    unknown_classifications.append(d)

    return (ref_classifications,
            alt_classifications,
            missing_classifications,
            unknown_classifications)

def save_classification_matrices(
    save_dir,
    union_variants,
    barcode_path=None,
    ref_classifications=None,
    alt_classifications=None,
    missing_classifications=None,
    unknown_classifications=None
):
    variants = [f"{v['CHROM']}:{v['POS']}_{v['REF']}->{v['ALT']}" for v in union_variants]
    variant_to_index = {v: i for i, v in enumerate(variants)}

    all_barcodes = set()
    def gather_barcodes(classifs):
        if not classifs: return
        for c in classifs:
            all_barcodes.add(c["cell_barcode"])

    gather_barcodes(ref_classifications)
    gather_barcodes(alt_classifications)
    gather_barcodes(missing_classifications)
    gather_barcodes(unknown_classifications)
    barcodes = sorted(all_barcodes)
    barcode_to_index = {bc: i for i, bc in enumerate(barcodes)}

    def build_csr_and_save(classifs, name):
        from collections import defaultdict
        if not classifs: return
        data = []
        row_inds = []
        col_inds = []
        count_dict = defaultdict(int)
        for c in classifs:
            var = c["variant"]
            bc  = c["cell_barcode"]
            if var not in variant_to_index or bc not in barcode_to_index:
                continue
            r = variant_to_index[var]
            co = barcode_to_index[bc]
            count_dict[(r, co)] += 1

        for (r, co), val in count_dict.items():
            row_inds.append(r)
            col_inds.append(co)
            data.append(val)

        mat = csr_matrix((data, (row_inds, col_inds)), shape=(len(variants), len(barcodes)))
        with h5py.File(os.path.join(save_dir, f"{name}_matrix.h5"), "w") as hf:
            hf.create_dataset("data", data=mat.data)
            hf.create_dataset("indices", data=mat.indices)
            hf.create_dataset("indptr", data=mat.indptr)
            hf.attrs['shape'] = mat.shape

    build_csr_and_save(ref_classifications, "ref")
    build_csr_and_save(alt_classifications, "alt")
    build_csr_and_save(missing_classifications, "missing")
    build_csr_and_save(unknown_classifications, "unknown")

    joblib.dump((variants, barcodes), os.path.join(save_dir, "variant_barcode_mappings.pkl"))

#####################################
# 5) 실제 chunk 처리 함수 (병렬용)
#####################################

def process_one_chunk(
    chrom, win_start, win_end,
    vcf_files_with_origins, bam_path,
    padding=3000,
    compute_missing_unknown=True,
    save_dir=None
):
    """
    한 개 chunk (chrom, win_start, win_end)에 대해
    1) variant fetch
    2) reads fetch (with padding)
    3) overlap & classification
    4) 임시결과 dict return
       (원한다면 여기서 바로 joblib.dump() 할 수도 있음)
    """
    logging.info(f"[Chunk] {chrom}:{win_start}-{win_end} start")

    # 1) Variants
    union_variants = process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    if not union_variants:
        logging.info(f"[Chunk] {chrom}:{win_start}-{win_end}: No variants, skip.")
        return None

    # 2) Reads
    df_reads = extract_reads_info_region_with_padding(bam_path, chrom, win_start, win_end, padding=padding)
    if df_reads.empty:
        logging.info(f"[Chunk] {chrom}:{win_start}-{win_end}: No reads, skip.")
        return None

    # 3) overlap & classification
    read_mapping = create_read_mapping(df_reads)
    df_reads, selected_read_unique_names, final_variants = process_reads_and_variants(
        df_reads, union_variants, win_start, win_end
    )
    refc, altc, missingc, unknownc = classify_in_parallel(
        bam_path,
        selected_read_unique_names,
        final_variants,
        read_mapping,
        compute_missing_unknown=compute_missing_unknown
    )

    # 결과를 dict 형태로 반환
    chunk_data = {
        "union_variants": final_variants,
        "ref": refc,
        "alt": altc,
        "missing": missingc,
        "unknown": unknownc
    }
    logging.info(f"[Chunk] {chrom}:{win_start}-{win_end} done => {len(final_variants)} variants, #reads={len(df_reads)}")

    return chunk_data

#####################################
# 6) MAIN: 윈도우+패딩 + chunk 병렬
#####################################

def main():
    args, vcf_files_with_origins = parse_arguments()

    setup_logging(args.save_dir)
    num_cores = check_num_cores(args.num_cores)
    compute_missing_unknown = not args.ref_alt_only

    # 1) chrom_lengths
    chrom_lengths = get_chromosome_lengths(args.bam_path, args.chromosomes)
    if not chrom_lengths:
        logging.error("No chromosomes found to process.")
        sys.exit(1)

    window_size = args.window_size
    padding = args.padding
    logging.info(f"Using window_size={window_size}, padding={padding}, num_cores={num_cores}")

    os.makedirs(args.save_dir, exist_ok=True)

    # 2) 모든 (chrom, win_start, win_end) 구간을 리스트업
    all_chunks = []
    for chrom, length in chrom_lengths.items():
        for (win_start, win_end) in get_windows_for_chromosome(chrom, length, window_size):
            all_chunks.append((chrom, win_start, win_end))

    logging.info(f"Total chunk count: {len(all_chunks)}")

    # 3) 병렬로 각 chunk 처리
    chunk_results = Parallel(n_jobs=num_cores)(
        delayed(process_one_chunk)(
            chrom=c[0], win_start=c[1], win_end=c[2],
            vcf_files_with_origins=vcf_files_with_origins,
            bam_path=args.bam_path,
            padding=padding,
            compute_missing_unknown=compute_missing_unknown
        ) for c in all_chunks
    )

    # 4) 결과 머지
    final_union_variants = []
    final_ref = []
    final_alt = []
    final_missing = []
    final_unknown = []

    for res in chunk_results:
        if res is None:
            continue
        final_union_variants.extend(res["union_variants"])
        final_ref.extend(res["ref"])
        final_alt.extend(res["alt"])
        if compute_missing_unknown:
            final_missing.extend(res["missing"])
            final_unknown.extend(res["unknown"])

    # 5) save
    save_classification_matrices(
        save_dir=args.save_dir,
        union_variants=final_union_variants,
        barcode_path=args.barcode_path,
        ref_classifications=final_ref,
        alt_classifications=final_alt,
        missing_classifications=final_missing if compute_missing_unknown else None,
        unknown_classifications=final_unknown if compute_missing_unknown else None
    )

    logging.info("=== All done ===")


if __name__ == '__main__':
    main()
