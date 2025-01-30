# scVarID_windowed_pipeline.py

import os
import sys
import time
import logging
import argparse
import pysam
import joblib
import pandas as pd
from collections import defaultdict
from joblib import Parallel, delayed
from classification import (
    process_bam_data_parallel,
    create_variants_dict,
    split_list
)

# 윈도우 생성 함수
def generate_windows(chrom_lengths, window_size, padding=5000):
    windows = []
    for chrom, length in chrom_lengths.items():
        start = 1
        while start <= length:
            end = min(start + window_size - 1, length)
            windows.append((chrom, start, end))
            start = end + 1
    return windows

# 변이 데이터 추출 함수 (VCF/TXT 파일 처리)
def extract_variants(vcf_files_with_origins, chrom, start, end):
    all_variants = []
    for file_path, origin in vcf_files_with_origins:
        if file_path.endswith('.vcf') or file_path.endswith('.vcf.gz') or file_path.endswith('.bcf'):
            try:
                with pysam.VariantFile(file_path) as vf:
                    for record in vf.fetch(chrom, start, end):
                        if "PASS" in record.filter.keys():
                            variant = {
                                "CHROM": record.chrom,
                                "POS": record.pos,
                                "REF": record.ref,
                                "ALT": ",".join(record.alts) if record.alts else "",
                                "ORIGIN": origin
                            }
                            all_variants.append(variant)
            except Exception as e:
                logging.warning(f"VCF 파일 {file_path}에서 변이 추출 중 에러 발생: {e}")
        elif file_path.endswith('.txt'):
            try:
                df = pd.read_csv(file_path, sep='\t', header=None, names=["CHROM", "POS", "REF", "ALT"],
                                 dtype={"CHROM": str, "POS": int, "REF": str, "ALT": str})
                subset = df[(df["CHROM"] == chrom) & (df["POS"] >= start) & (df["POS"] <= end)]
                for _, row in subset.iterrows():
                    variant = {
                        "CHROM": row["CHROM"],
                        "POS": row["POS"],
                        "REF": row["REF"],
                        "ALT": row["ALT"],
                        "ORIGIN": origin
                    }
                    all_variants.append(variant)
            except Exception as e:
                logging.warning(f"TXT 파일 {file_path}에서 변이 추출 중 에러 발생: {e}")
        else:
            logging.warning(f"지원되지 않는 파일 형식: {file_path}")
    return all_variants

# 윈도우 단위로 처리하는 함수
def process_window(chrom, start, end, vcf_files_with_origins, bam_path, save_dir, padding=5000):
    logging.info(f"윈도우 처리 시작: {chrom}:{start}-{end}")
    
    # 1. 변이 데이터 추출
    variants = extract_variants(vcf_files_with_origins, chrom, start, end)
    if not variants:
        logging.info(f"윈도우 {chrom}:{start}-{end}에 변이 없음. 스킵합니다.")
        return
    
    # 2. 변이 딕셔너리 생성
    variants_dict = create_variants_dict(variants)
    
    # 3. BAM 파일에서 리드 추출 (패딩 포함)
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            reads = bam_file.fetch(chrom, max(1, start - padding) - 1, end + padding)
            read_data = []
            read_name_counts = defaultdict(int)
            for read in reads:
                read_name = read.query_name
                read_name_counts[read_name] += 1
                unique_id = read_name_counts[read_name] - 1
                read_unique_name = f"{read_name}_{unique_id}"
                read_data.append(read_unique_name)
    except Exception as e:
        logging.error(f"BAM 파일 {bam_path}에서 리드 추출 중 에러 발생: {e}")
        return
    
    if not read_data:
        logging.info(f"윈도우 {chrom}:{start}-{end}에 리드 없음. 스킵합니다.")
        return
    
    # 4. 리드 병렬 처리
    ref, alt, missing, unknown = process_bam_data_parallel(
        bam_path=bam_path,
        selected_read_unique_names=read_data,
        selected_variants=variants,
        read_mapping=defaultdict(list),  # 필요에 따라 리드 매핑 조정
        num_cores=4,
        compute_missing_unknown=True
    )
    
    # 5. 결과 저장
    window_save_dir = os.path.join(save_dir, f"{chrom}_{start}_{end}")
    os.makedirs(window_save_dir, exist_ok=True)
    joblib.dump({
        "variants": variants,
        "ref": ref,
        "alt": alt,
        "missing": missing,
        "unknown": unknown
    }, os.path.join(window_save_dir, "results.pkl"))
    
    logging.info(f"윈도우 처리 완료: {chrom}:{start}-{end}")

def main():
    # 인자 파싱
    parser = argparse.ArgumentParser(description='scVarID 윈도우 기반 파이프라인.')
    parser.add_argument('--variant-files', nargs='+', required=True,
                        help='SAMPLE, FILE의 페어 목록 (예: sample1 file1.vcf sample2 file2.vcf)')
    parser.add_argument('--bam-path', required=True, help='BAM 파일 경로.')
    parser.add_argument('--save-dir', required=True, help='결과를 저장할 디렉토리.')
    parser.add_argument('--num-cores', type=int, default=4, help='병렬 처리 코어 수.')
    parser.add_argument('--window-size', type=int, default=1000000, help='윈도우 크기 (bp).')
    parser.add_argument('--padding', type=int, default=5000, help='윈도우 경계 주변 패딩 (bp).')
    args = parser.parse_args()
    
    # 로깅 설정
    os.makedirs(args.save_dir, exist_ok=True)
    log_file = os.path.join(args.save_dir, 'processing.log')
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.DEBUG,  # 상세 로그를 위해 DEBUG로 설정
        format='%(asctime)s:%(levelname)s:%(message)s'
    )
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)  # 콘솔에는 INFO 이상만 표시
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)
    
    logging.info("=== scVarID 윈도우 기반 파이프라인 시작 ===")
    
    # 변이 파일 페어 파싱
    if len(args.variant_files) % 2 != 0:
        logging.error("각 SAMPLE은 FILE과 페어로 제공되어야 합니다.")
        sys.exit(1)
    
    vcf_files_with_origins = []
    for i in range(0, len(args.variant_files), 2):
        sample = args.variant_files[i]
        file_path = args.variant_files[i+1]
        vcf_files_with_origins.append((file_path, sample))
    
    # 염색체 길이 추출 (BAM 파일에서)
    try:
        with pysam.AlignmentFile(args.bam_path, "rb") as bam_file:
            chrom_lengths = dict(zip(bam_file.references, bam_file.lengths))
    except Exception as e:
        logging.error(f"BAM 파일 {args.bam_path} 열기 실패: {e}")
        sys.exit(1)
    
    # 윈도우 생성
    windows = generate_windows(chrom_lengths, args.window_size, args.padding)
    logging.info(f"총 윈도우 개수: {len(windows)}")
    
    # 병렬로 윈도우 처리
    Parallel(n_jobs=args.num_cores)(
        delayed(process_window)(
            chrom=win[0],
            start=win[1],
            end=win[2],
            vcf_files_with_origins=vcf_files_with_origins,
            bam_path=args.bam_path,
            save_dir=args.save_dir,
            padding=args.padding
        ) for win in windows
    )
    
    logging.info("=== scVarID 윈도우 기반 파이프라인 완료 ===")

if __name__ == "__main__":
    main()
