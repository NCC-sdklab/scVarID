# main.py

import os
import logging
from utils import (
    check_num_cores, 
    save_classification_matrices, 
    setup_logging, 
    log_step_start, 
    log_step_end,
    parse_arguments,
    get_all_chromosomes  # 염색체 목록 추출 함수 추가
)
from variant_processing import process_vcf_files, extract_variants_for_chromosome  # 염색체별 변이 처리 함수 추가
from read_processing import (
    extract_reads_for_chromosome,  # 염색체별 리드 처리 함수 추가
    create_read_mapping, 
    process_reads_and_variants, 
    process_variants_and_reads
)
from classification import process_bam_data_parallel  # 병렬 처리는 기존 유지
from joblib import Parallel, delayed  # 병렬 처리를 위한 라이브러리

def process_chromosome(chromosome, args, vcf_files_with_origins):
    """
    Process a single chromosome and save results.
    """
    logging.info(f"=== Start processing chromosome {chromosome} ===")

    # 1. Process variants for the current chromosome
    union_variants = []
    for file_path, origin in vcf_files_with_origins:
        variants = extract_variants_for_chromosome(file_path, origin, chromosome)
        union_variants.extend(variants)
    
    # 2. Process BAM file for the current chromosome
    df_reads = extract_reads_for_chromosome(args.bam_path, chromosome)
    read_mapping = create_read_mapping(df_reads)

    # 3. Classification process
    df_reads, selected_read_unique_names = process_reads_and_variants(df_reads, union_variants)
    updated_union_variants, selected_variants = process_variants_and_reads(union_variants, df_reads)

    ref_classifications, alt_classifications, missing_classifications, unknown_classifications = process_bam_data_parallel(
        args.bam_path,
        selected_read_unique_names,
        selected_variants,
        read_mapping,
        num_cores=args.num_cores,
        compute_missing_unknown=not args.ref_alt_only
    )

    # 4. Save results for the current chromosome
    save_classification_matrices(
        save_dir=os.path.join(args.save_dir, f"chromosome_{chromosome}"),
        union_variants=updated_union_variants,
        barcode_path=args.barcode_path,
        ref_classifications=ref_classifications,
        alt_classifications=alt_classifications,
        missing_classifications=missing_classifications if not args.ref_alt_only else None,
        unknown_classifications=unknown_classifications if not args.ref_alt_only else None
    )

    logging.info(f"=== Finished processing chromosome {chromosome} ===\n")


def main():
    # Parse arguments
    args, vcf_files_with_origins = parse_arguments()

    # Assign other arguments to variables
    bam_path = args.bam_path
    save_dir = args.save_dir
    num_cores = args.num_cores

    # Create save directory
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # Configure logging
    setup_logging(save_dir)

    # Verify and adjust num_cores
    num_cores = check_num_cores(num_cores)
    
    # Log program start
    logging.info("=== [scVarID] Program Started ===\n")
    
    try:
        # Extract list of chromosomes
        chromosomes = args.chromosomes or get_all_chromosomes(bam_path, [file_path for file_path, _ in vcf_files_with_origins])

        # Process each chromosome in parallel
        Parallel(n_jobs=num_cores)(
            delayed(process_chromosome)(chromosome, args, vcf_files_with_origins) for chromosome in chromosomes
        )
    
    except Exception as e:
        logging.error(f"An error occurred: {e}", exc_info=True)
        exit(1)
    
    # Log program completion
    logging.info("=== All steps completed successfully ===")

if __name__ == '__main__':
    main()
