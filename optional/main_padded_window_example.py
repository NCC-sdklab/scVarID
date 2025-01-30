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

#####################################
# 1) 유틸 함수들
#####################################

def parse_arguments():
    """
    기존 utils.py의 parse_arguments()와 유사하되
    --window-size, --padding 같은 인자를 추가.
    """
    parser = argparse.ArgumentParser(description='scVarID-like window+padding example.')

    # 기존
    parser.add_argument('--variant-files', nargs='+', required=True, metavar=('ORIGIN', 'FILE_PATH'),
                        help='List of origin and variant file path pairs.')
    parser.add_argument('--bam-path', required=True, help='Path to the BAM file.')
    parser.add_argument('--barcode-path', required=False, help='Path to the barcode file (optional).')
    parser.add_argument('--save-dir', required=True, help='Save directory path.')
    parser.add_argument('--num-cores', type=int, default=4, help='Number of cores to use for parallel processing.')
    parser.add_argument('--ref-alt-only', action='store_true',
                        help='Compute only ref and alt classifications, skipping missing and unknown.')
    parser.add_argument('--chromosomes', nargs='+', required=False,
                        help='List of chromosomes to process. If not specified, process all from BAM header.')

    # **새로 추가**: window-size, padding
    parser.add_argument('--window-size', type=int, default=10_000_000,
                        help='Window size (bp) for chunk processing. Default=10,000,000.')
    parser.add_argument('--padding', type=int, default=500,
                        help='Padding size (bp) to fetch more reads around chunk boundary. Default=500.')

    args = parser.parse_args()

    # 처리
    # (scVarID에서는 variant-files를 2개씩 묶어 (origin, filepath)로 파싱했었음)
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
    로그 세팅.
    """
    os.makedirs(save_dir, exist_ok=True)
    log_file = os.path.join(save_dir, 'processing.log')

    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s'
    )
    # 콘솔에도 찍으려면 추가 핸들러
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


def check_num_cores(num_cores):
    import multiprocessing
    max_cores = multiprocessing.cpu_count()
    if num_cores > max_cores:
        logging.info(f"Warning: requested {num_cores} cores > available {max_cores}. Use {max_cores} instead.")
        return max_cores
    return num_cores


def get_chromosome_lengths(bam_path, target_chromosomes=None):
    """
    BAM 파일 헤더에서 (chrom -> length) 얻기.
    target_chromosomes 지정되면 그중 유효한 것만 반환.
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
    chrom_length 범위를 window_size 씩 나누어 (start, end) yield.
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
    """Check if the file is in VCF or BCF format."""
    return file_path.endswith('.vcf') or file_path.endswith('.vcf.gz') or file_path.endswith('.bcf')


def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    """
    (chrom, start, end)에 해당하는 variant만 추출해서 union_variants 반환.
    """
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
    """
    pysam.VariantFile.fetch(chrom, start, end)로 (chrom, pos in [start, end])만 추출
    """
    import pysam
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
    """
    TXT 파일은 fetch가 안 되므로 전체 읽고, (CHROM==chrom, start<=POS<=end) 필터
    """
    import pandas as pd
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
    """
    scVarID식 중복 제거 & VINDEX 부여. (간단 버전)
    """
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

def extract_reads_info_region_with_padding(bam_path, chrom, start, end, padding=500):
    """
    (chrom, start, end)에 대해서, BAM을 fetch할 때 패딩을 주어
    경계에 걸친 read도 온전히 읽는다.
    - fetch는 0-based exclusive. start/end가 1-based이면 start-1 해야 함
    """
    import pysam
    data = []
    read_name_counts = {}

    fetch_start = max(1, start - padding)
    fetch_end   = end + padding

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # fetch는 0-based (start inclusive, end exclusive)
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
    logging.info(f"[WithPadding] Extracted {len(df)} reads in region {chrom}:{start}-{end} (padding={padding}).")

    return df


#####################################
# 4) 간단한 Overlap / Classification (데모용)
#####################################

def create_read_mapping(df_reads):
    """
    read_name -> [read_unique_name1, read_unique_name2, ...]
    """
    g = df_reads.groupby("Read_Name")["Read_Unique_Name"].apply(list)
    return g.to_dict()


def process_reads_and_variants(df_reads, union_variants, chunk_start, chunk_end):
    """
    기존 scVarID 스타일로 read-variant 오버랩 정보 추가.
    여기서는 'missing' 판단 시 chunk 범위를 벗어나는 variant는 skip.
    """
    # (여기서는 간단히 'selected_read_unique_names' = df_reads의 전체
    #  실제론 variant-pos와 read range가 겹치는지 체크 가능)
    selected_read_unique_names = df_reads["Read_Unique_Name"].tolist()

    # **변이 중에서 pos가 [chunk_start, chunk_end] 범위 벗어나는건 skip 가능**
    # 근데 이미 region fetch에서 variant pos가 [start, end]였으니 크게 문제 없을 수도.
    # 만약 패딩으로 variant도 fetch했다면 여기서 필터 필요:
    final_variants = []
    for v in union_variants:
        pos = v["POS"]
        if pos < chunk_start or pos > chunk_end:
            # skip
            continue
        final_variants.append(v)

    return df_reads, selected_read_unique_names, final_variants


def classify_in_parallel(bam_path, selected_read_unique_names, selected_variants, read_mapping,
                         compute_missing_unknown=True):
    """
    여기서는 실제 scVarID classification보다 간단히 가짜 데이터 생성.
    (실전에서는 classification.process_bam_data_parallel() 사용)
    """
    # 예: 변이개수 * read개수 * "missing" 임의 카운트..
    # 여기서는 'missing' 한 번만 보여주기 위해 더미로 만든다.
    import random
    ref_classifications = []
    alt_classifications = []
    missing_classifications = []
    unknown_classifications = []

    for rn in selected_read_unique_names:
        # 바코드 추출
        # read_mapping: read_name -> [unique_ids]
        # 여기서는 그냥 dummy
        for var in selected_variants:
            classification = random.choice(["ref", "alt", "missing", "unknown"])
            dat = {
                "read_name": rn,
                "cell_barcode": "BC_TEST",  # dummy
                "variant": f"{var['CHROM']}:{var['POS']}_{var['REF']}->{var['ALT']}"
            }
            if classification == "ref":
                ref_classifications.append(dat)
            elif classification == "alt":
                alt_classifications.append(dat)
            elif classification == "missing":
                if compute_missing_unknown:
                    missing_classifications.append(dat)
            else: # unknown
                if compute_missing_unknown:
                    unknown_classifications.append(dat)

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
    """
    scVarID의 save_classification_matrices() 유사.
    (여기서는 간단화 버전)
    """
    import joblib
    from collections import defaultdict
    from scipy.sparse import csr_matrix

    # variants
    variants = [f"{v['CHROM']}:{v['POS']}_{v['REF']}->{v['ALT']}" for v in union_variants]
    variant_to_index = {v: i for i, v in enumerate(variants)}

    # 바코드
    # 실제론 barcode_path가 있으면 읽고, 없으면 classification에서 나온 barcode들의 union
    # 여기서는 단순히 classification에서 barcode 모으기
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

    # build & save each classification
    def build_csr_and_save(classifs, name):
        if not classifs or len(classifs) == 0:
            return
        data = []
        row_inds = []
        col_inds = []
        from collections import defaultdict
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

    # 마지막에 variant, barcode 저장
    joblib.dump((variants, barcodes), os.path.join(save_dir, "variant_barcode_mappings.pkl"))


#####################################
# 5) MAIN: 윈도우+패딩
#####################################

def main():
    args, vcf_files_with_origins = parse_arguments()

    setup_logging(args.save_dir)
    num_cores = check_num_cores(args.num_cores)
    compute_missing_unknown = not args.ref_alt_only

    # 1) BAM 헤더에서 (chrom -> length)
    chrom_lengths = get_chromosome_lengths(args.bam_path, args.chromosomes)
    if not chrom_lengths:
        logging.error("No chromosomes found to process.")
        sys.exit(1)

    window_size = args.window_size
    padding = args.padding
    logging.info(f"Using window_size={window_size}, padding={padding}")

    os.makedirs(args.save_dir, exist_ok=True)
    tmp_dir = os.path.join(args.save_dir, "tmp_chunks")
    os.makedirs(tmp_dir, exist_ok=True)

    # 2) Chromosome → window loop
    for chrom, length in chrom_lengths.items():
        logging.info(f"Processing chrom={chrom}, length={length}")
        for (win_start, win_end) in get_windows_for_chromosome(chrom, length, window_size):
            logging.info(f"Chunk region: {chrom}:{win_start}-{win_end}")

            # 2-1. Variant fetch (POS in [win_start, win_end])
            union_variants = process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
            if not union_variants:
                logging.info("No variants in this window, skip.")
                continue

            # 2-2. Reads with padding
            df_reads = extract_reads_info_region_with_padding(
                args.bam_path, chrom, win_start, win_end, padding=padding
            )
            if df_reads.empty:
                logging.info("No reads found in this window, skip.")
                continue

            # 2-3. overlap & classification
            read_mapping = create_read_mapping(df_reads)
            df_reads, selected_read_unique_names, final_variants = process_reads_and_variants(
                df_reads, union_variants, win_start, win_end
            )

            refc, altc, missingc, unknownc = classify_in_parallel(
                args.bam_path,
                selected_read_unique_names,
                final_variants,
                read_mapping,
                compute_missing_unknown=compute_missing_unknown
            )

            # 2-4. 임시 저장
            out_path = os.path.join(tmp_dir, f"{chrom}_{win_start}_{win_end}.pkl")
            joblib.dump({
                "union_variants": final_variants,
                "ref": refc,
                "alt": altc,
                "missing": missingc,
                "unknown": unknownc
            }, out_path)

    # 3) 임시 파일 merge
    final_union_variants = []
    final_ref = []
    final_alt = []
    final_missing = []
    final_unknown = []

    for fname in os.listdir(tmp_dir):
        if fname.endswith(".pkl"):
            cdata = joblib.load(os.path.join(tmp_dir, fname))
            final_union_variants.extend(cdata["union_variants"])
            final_ref.extend(cdata["ref"])
            final_alt.extend(cdata["alt"])
            if compute_missing_unknown:
                final_missing.extend(cdata["missing"])
                final_unknown.extend(cdata["unknown"])

    # 4) 최종 저장
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


if __name__ == "__main__":
    main()
