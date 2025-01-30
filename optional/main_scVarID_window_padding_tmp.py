#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import time
import logging
import argparse
import pysam
import joblib
import pandas as pd
from collections import defaultdict
from scipy.sparse import csr_matrix

from joblib import Parallel, delayed

############################################################################
# 0. 유틸 함수: Step 로그(Start/End), 시간 경과 포맷
############################################################################

def log_step_start(step_name):
    """
    스텝 시작 로그를 찍고, 시작 시간을 리턴한다.
    """
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    """
    스텝 종료 로그를 찍고, 걸린 시간(hh:mm:ss)을 로깅한다.
    """
    end_time = time.time()
    elapsed = end_time - start_time
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {h:02d}h {m:02d}m {s:02d}s ===")

############################################################################
# 1. 인자 파싱 & 로깅 설정
############################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description='scVarID-like pipeline with window+padding+parallel + scVarID-style logs, dumping chunk results to disk.')
    # 기존 인자
    parser.add_argument('--variant-files', nargs='+', required=True, metavar=('ORIGIN', 'FILE_PATH'),
                        help='List of origin and variant file path pairs.')
    parser.add_argument('--bam-path', required=True, help='Path to the BAM file.')
    parser.add_argument('--barcode-path', required=False, help='Path to the barcode file (optional).')
    parser.add_argument('--save-dir', required=True, help='Directory where results will be saved.')
    parser.add_argument('--num-cores', type=int, default=4, help='Number of CPU cores for parallel processing.')
    parser.add_argument('--ref-alt-only', action='store_true',
                        help='Compute only ref and alt classifications, skipping missing and unknown.')
    parser.add_argument('--chromosomes', nargs='+', required=False,
                        help='List of chromosomes to process (e.g. chr1 chr2). If not specified, process all from BAM.')

    # 윈도우+패딩 추가
    parser.add_argument('--window-size', type=int, default=10_000_000,
                        help='Window size (bp) for chunk processing (default=10Mb).')
    parser.add_argument('--padding', type=int, default=500,
                        help='Padding size (bp) around chunk boundary (default=500).')

    args = parser.parse_args()

    # variant-files를 (file_path, origin) 형태로 묶기
    if len(args.variant_files) % 2 != 0:
        parser.error('Each origin must be paired with a file path in --variant-files.')

    vcf_files_with_origins = []
    for i in range(0, len(args.variant_files), 2):
        origin = args.variant_files[i]
        file_path = args.variant_files[i+1]
        vcf_files_with_origins.append((file_path, origin))

    return args, vcf_files_with_origins

def setup_logging(save_dir):
    """
    scVarID-style 로깅 초기화.
    """
    os.makedirs(save_dir, exist_ok=True)
    log_file = os.path.join(save_dir, 'processing.log')
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s'
    )
    # 콘솔에도 로그 출력
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

############################################################################
# 2. Chromosome lengths, windows
############################################################################

def get_chromosome_lengths(bam_path, target_chromosomes=None):
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
    start = 1
    while start <= chrom_length:
        end = min(start + window_size - 1, chrom_length)
        yield (start, end)
        start = end + 1

############################################################################
# 3. Variant 추출 & 필터
############################################################################

def is_vcf_file(file_path):
    return file_path.endswith('.vcf') or file_path.endswith('.vcf.gz') or file_path.endswith('.bcf')

def extract_variants_vcf_region(vcf_file, origin, chrom, start, end):
    import pysam
    vs = []
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
                    vs.append(v)
    except Exception as e:
        logging.warning(f"Error reading {vcf_file} region: {e}")
    return vs

def extract_variants_txt_region(txt_file, origin, chrom, start, end):
    vs = []
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
            vs.append(v)
    except Exception as e:
        logging.warning(f"Error reading {txt_file} region: {e}")
    return vs

def filter_variants(variants_list):
    union_dict = {}
    for v in variants_list:
        key = (v["CHROM"], v["POS"], v["REF"], v["ALT"])
        if key not in union_dict:
            union_dict[key] = v
        else:
            old_origins = union_dict[key]["ORIGIN"].split(",")
            new_origins = v["ORIGIN"].split(",")
            union_dict[key]["ORIGIN"] = ",".join(set(old_origins + new_origins))

    union_variants = list(union_dict.values())
    for i, var in enumerate(union_variants, 1):
        var["VINDEX"] = f"vi_{i}"
    return union_variants

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    all_variants = []
    for file_path, origin in vcf_files_with_origins:
        if is_vcf_file(file_path):
            chunk_vs = extract_variants_vcf_region(file_path, origin, chrom, start, end)
        elif file_path.endswith('.txt'):
            chunk_vs = extract_variants_txt_region(file_path, origin, chrom, start, end)
        else:
            logging.warning(f"Unsupported file format: {file_path}")
            continue
        all_variants.extend(chunk_vs)

    union_variants = filter_variants(all_variants)
    return union_variants

############################################################################
# 4. Reads 추출 (with padding), Overlap, Classification
############################################################################

def extract_reads_info_region_with_padding(bam_path, chrom, start, end, padding=500):
    data = []
    read_name_counts = {}
    fetch_start = max(1, start - padding)
    fetch_end   = end + padding

    with pysam.AlignmentFile(bam_path, "rb") as bam:
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
                read.reference_start,
                read.reference_end
            ])

    df = pd.DataFrame(data, columns=[
        "Read_Name", "Read_Unique_Name", "Cell_Barcode", "Chromosome",
        "Start_Position", "End_Position"
    ])
    return df

def create_read_mapping(df_reads):
    g = df_reads.groupby("Read_Name")["Read_Unique_Name"].apply(list)
    return g.to_dict()

def process_reads_and_variants(df_reads, union_variants, chunk_start, chunk_end):
    selected_read_unique_names = df_reads["Read_Unique_Name"].tolist()
    final_variants = []
    for v in union_variants:
        pos = v["POS"]
        if chunk_start <= pos <= chunk_end:
            final_variants.append(v)
    return df_reads, selected_read_unique_names, final_variants

def classify_in_parallel(selected_read_unique_names, selected_variants, read_mapping, compute_missing_unknown=True):
    """
    예시로 random 분류.
    """
    import random
    ref_list = []
    alt_list = []
    missing_list = []
    unknown_list = []

    for rn in selected_read_unique_names:
        for var in selected_variants:
            classification = random.choice(["ref", "alt", "missing", "unknown"])
            rec = {
                "read_name": rn,
                "cell_barcode": "BC_TEST",
                "variant": f"{var['CHROM']}:{var['POS']}_{var['REF']}->{var['ALT']}"
            }
            if classification == "ref":
                ref_list.append(rec)
            elif classification == "alt":
                alt_list.append(rec)
            elif classification == "missing":
                if compute_missing_unknown:
                    missing_list.append(rec)
            else:
                if compute_missing_unknown:
                    unknown_list.append(rec)

    return (ref_list, alt_list, missing_list, unknown_list)

############################################################################
# 5. 결과 저장 (ref_matrix.h5, alt_matrix.h5 등)
############################################################################

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
        if not classifs:
            return
        for c in classifs:
            all_barcodes.add(c["cell_barcode"])

    gather_barcodes(ref_classifications)
    gather_barcodes(alt_classifications)
    gather_barcodes(missing_classifications)
    gather_barcodes(unknown_classifications)

    barcodes = sorted(all_barcodes)
    barcode_to_index = {bc: i for i, bc in enumerate(barcodes)}

    from scipy.sparse import csr_matrix
    import h5py

    def build_csr_and_save(classifs, name):
        if not classifs:
            return
        from collections import defaultdict
        data = []
        row_inds = []
        col_inds = []
        count_map = defaultdict(int)
        for c in classifs:
            var = c["variant"]
            bc  = c["cell_barcode"]
            if var not in variant_to_index or bc not in barcode_to_index:
                continue
            r = variant_to_index[var]
            co = barcode_to_index[bc]
            count_map[(r, co)] += 1

        for (r, co), val in count_map.items():
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

    build_csr_and_save(ref_classifications, "ref")
    build_csr_and_save(alt_classifications, "alt")
    build_csr_and_save(missing_classifications, "missing")
    build_csr_and_save(unknown_classifications, "unknown")

    joblib.dump((variants, barcodes), os.path.join(save_dir, "variant_barcode_mappings.pkl"))

############################################################################
# 6. Chunk 처리 함수: chunk별 .pkl 파일로 저장
############################################################################

def process_one_chunk_and_dump(
    chrom, win_start, win_end,
    vcf_files_with_origins,
    bam_path,
    padding=500,
    compute_missing_unknown=True,
    tmp_dir=None
):
    """
    (chrom, win_start, win_end) -> variant fetch -> read fetch + padding ->
    overlap/classification -> chunk 결과를 .pkl 파일로 저장 + 반환 (pkl 경로)

    메모리에서 chunk_data 제거 -> RAM 절약
    """
    chunk_str = f"{chrom}:{win_start}-{win_end}"
    logging.info(f"Processing chunk {chunk_str} ...")

    # 1) Variants
    chunk_variants = process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    logging.info(f"Extracted {len(chunk_variants)} variants from region {chunk_str}")

    if not chunk_variants:
        return None  # skip

    # 2) Reads
    df_reads = extract_reads_info_region_with_padding(bam_path, chrom, win_start, win_end, padding=padding)
    logging.info(f"Extracted {len(df_reads)} reads from region {chunk_str}")

    if df_reads.empty:
        return None

    # 3) Overlap
    read_mapping = create_read_mapping(df_reads)
    df_reads, selected_read_unique_names, final_variants = process_reads_and_variants(
        df_reads, chunk_variants, win_start, win_end
    )

    # 4) Classification
    refc, altc, missingc, unknownc = classify_in_parallel(
        selected_read_unique_names,
        final_variants,
        read_mapping,
        compute_missing_unknown=compute_missing_unknown
    )

    # mappable(placeholder)
    mappable_reads = len(selected_read_unique_names)
    total_reads = len(df_reads)
    perc_reads = 100.0 * mappable_reads / total_reads if total_reads > 0 else 0
    logging.info(f"# mappable reads: {mappable_reads} ({perc_reads:.2f}%) in chunk {chunk_str}")

    mappable_variants = len(final_variants)
    total_variants = len(chunk_variants)
    perc_vars = 100.0 * mappable_variants / total_variants if total_variants > 0 else 0
    logging.info(f"# mappable variants: {mappable_variants} ({perc_vars:.2f}%) in chunk {chunk_str}")

    # 5) chunk 결과 dict
    chunk_data = {
        "union_variants": final_variants,
        "ref": refc,
        "alt": altc,
        "missing": missingc,
        "unknown": unknownc
    }

    # 6) 임시 파일로 저장
    if tmp_dir is None:
        tmp_dir = "./tmp_chunks"
    os.makedirs(tmp_dir, exist_ok=True)

    out_pkl = os.path.join(tmp_dir, f"chunk_{chrom}_{win_start}_{win_end}.pkl")
    joblib.dump(chunk_data, out_pkl)
    logging.info(f"Chunk {chunk_str} saved to {out_pkl}")

    # 7) 메모리에서 chunk_data 내려놓기
    del chunk_data
    return out_pkl  # 임시파일 경로만 반환


############################################################################
# 7. MAIN: 병렬 + tmp 파일 merge
############################################################################

def main():
    args, vcf_files_with_origins = parse_arguments()

    setup_logging(args.save_dir)
    logging.info("Command Line: " + " ".join(sys.argv))
    logging.info("=== [scVarID] Program Started ===")

    num_cores = check_num_cores(args.num_cores)
    compute_missing_unknown = not args.ref_alt_only

    # Step: Chromosome length retrieval
    step_time = log_step_start("Chromosome length retrieval")
    chrom_lengths = get_chromosome_lengths(args.bam_path, args.chromosomes)
    if not chrom_lengths:
        logging.error("No valid chromosomes found; abort.")
        sys.exit(1)
    log_step_end("Chromosome length retrieval", step_time)

    window_size = args.window_size
    padding = args.padding
    logging.info(f"[INFO] Using window_size={window_size}, padding={padding}, num_cores={num_cores}")

    tmp_dir = os.path.join(args.save_dir, "tmp_chunks")  # 임시 파일 저장 폴더

    # Step: Generate chunk list
    step_time = log_step_start("Generate chunk list")
    all_chunks = []
    for chrom, length in chrom_lengths.items():
        for (win_start, win_end) in get_windows_for_chromosome(chrom, length, window_size):
            all_chunks.append((chrom, win_start, win_end))
    logging.info(f"Total chunk count: {len(all_chunks)}")
    log_step_end("Generate chunk list", step_time)

    # Step: Parallel chunk processing (.pkl dump)
    step_time = log_step_start("Parallel chunk processing")
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
    log_step_end("Parallel chunk processing", step_time)

    # Step: Merge chunk results
    step_time = log_step_start("Merge results")

    final_union_variants = []
    final_ref = []
    final_alt = []
    final_missing = []
    final_unknown = []

    # chunk_files에는 pkl 경로들이 들어있음
    for pkl_path in chunk_files:
        if pkl_path is None:
            continue
        chunk_data = joblib.load(pkl_path)
        final_union_variants.extend(chunk_data["union_variants"])
        final_ref.extend(chunk_data["ref"])
        final_alt.extend(chunk_data["alt"])
        if compute_missing_unknown:
            final_missing.extend(chunk_data["missing"])
            final_unknown.extend(chunk_data["unknown"])
        # 메모리에서 내려놓기
        del chunk_data

    log_step_end("Merge results", step_time)

    # Step: Save classification results
    step_time = log_step_start("Save classification results")
    save_classification_matrices(
        save_dir=args.save_dir,
        union_variants=final_union_variants,
        barcode_path=args.barcode_path,
        ref_classifications=final_ref,
        alt_classifications=final_alt,
        missing_classifications=final_missing if compute_missing_unknown else None,
        unknown_classifications=final_unknown if compute_missing_unknown else None
    )
    log_step_end("Save classification results", step_time)

    logging.info("=== All steps completed successfully ===")


if __name__ == '__main__':
    main()
