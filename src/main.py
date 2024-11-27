# main.py

import os
import logging
import argparse

from variant_processing import process_vcf_files
from read_processing import extract_reads_info, create_read_mapping, process_reads_and_variants
from read_processing import process_variants_and_reads
from classification import process_bam_data_parallel
from utils import check_num_cores, save_classification_matrices

def main():
    # ArgumentParser 객체 생성
    parser = argparse.ArgumentParser(description='Process variant files and BAM files.')
    
    # 변이 파일에 대한 인자 추가
    parser.add_argument('--variant-files', nargs='+', required=True, metavar=('ORIGIN', 'FILE_PATH'),
                        help='List of origin and variant file path pairs.')
    
    # BAM 파일, 바코드 파일, 저장 디렉토리에 대한 인자 추가
    parser.add_argument('--bam-path', required=True, help='Path to the BAM file.')
    parser.add_argument('--barcode-path', required=False, help='Path to the barcode file (optional).')
    parser.add_argument('--save-dir', required=True, help='Save directory path.')
    
    # 코어 수 인자 추가
    parser.add_argument('--num-cores', type=int, default=4, help='Number of cores to use for parallel processing.')
    
    # action 추가
    parser.add_argument('--ref-alt-only', action='store_true',
                        help='Compute only ref and alt classifications, skipping missing and unknown.')
    
    # 인자 파싱
    args = parser.parse_args()
    
    # 변이 파일 목록 생성
    vcf_files_with_origins = []
    if len(args.variant_files) % 2 != 0:
        parser.error('Each origin must be paired with a file path in --variant-files.')
    
    for i in range(0, len(args.variant_files), 2):
        origin = args.variant_files[i]
        file_path = args.variant_files[i + 1]
        vcf_files_with_origins.append((file_path, origin))
    
    # 다른 인자들을 변수에 할당
    bam_path = args.bam_path
    barcode_path = args.barcode_path
    save_dir = args.save_dir
    num_cores = args.num_cores
    
    # compute_missing_unknown 값을 설정
    compute_missing_unknown = not args.ref_alt_only

    # 저장 디렉토리 생성
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # 로그 설정
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s',
        handlers=[
            logging.FileHandler(os.path.join(save_dir, 'processing.log')),
            logging.StreamHandler()
        ]
    )
    
    # 0. num_cores를 확인하고 조정
    num_cores = check_num_cores(num_cores)
    
    # 1. 변이 파일 처리
    union_variants = process_vcf_files(vcf_files_with_origins)
    
    # 2. BAM 파일에서 리드 정보 추출
    df_reads = extract_reads_info(bam_path)
    
    # 3. Read_Name과 Read_Unique_Name의 매핑 딕셔너리 생성
    read_mapping = create_read_mapping(df_reads)
    
    # 4. 리드와 변이 간의 오버랩 정보 추가
    df_reads, selected_read_unique_names = process_reads_and_variants(df_reads, union_variants)
    
    # 5. 변이와 리드 간의 오버랩 정보 추가
    updated_union_variants, selected_variants = process_variants_and_reads(union_variants, df_reads)
    
    # 6. 병렬 처리 함수 호출하여 분류 작업 수행
    ref_classifications, alt_classifications, missing_classifications, unknown_classifications = process_bam_data_parallel(
        bam_path,
        selected_read_unique_names,
        selected_variants,
        read_mapping,
        num_cores=num_cores,
        compute_missing_unknown=compute_missing_unknown
    )
    
    # 7. 분류 결과 저장
    save_classification_matrices(
        save_dir=save_dir,
        union_variants=union_variants,
        barcode_path=barcode_path,
        ref_classifications=ref_classifications,
        alt_classifications=alt_classifications,
        missing_classifications=missing_classifications if compute_missing_unknown else None,
        unknown_classifications=unknown_classifications if compute_missing_unknown else None
    )

if __name__ == '__main__':
    main()
