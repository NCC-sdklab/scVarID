# main.py

import os
import logging
from utils import (
    check_num_cores, 
    save_classification_matrices, 
    setup_logging, 
    log_step_start, 
    log_step_end,
    parse_arguments
)
from variant_processing import process_vcf_files
from read_processing import (
    extract_reads_info, 
    create_read_mapping, 
    process_reads_and_variants, 
    process_variants_and_reads
)
from classification import process_bam_data_parallel

def main():
    # Parse arguments
    args, vcf_files_with_origins = parse_arguments()
    
    # Assign other arguments to variables
    bam_path = args.bam_path
    barcode_path = args.barcode_path
    save_dir = args.save_dir
    num_cores = args.num_cores
    chromosomes = args.chromosomes  # New argument for chromosomes

    # Set compute_missing_unknown value
    compute_missing_unknown = not args.ref_alt_only

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
        # 1. Process variant files
        step_name = "Process variant files"
        start_time = log_step_start(step_name)
        union_variants = process_vcf_files(vcf_files_with_origins, chromosomes=chromosomes)
        log_step_end(step_name, start_time)
        
        # 2. Process BAM file
        step_name = "Process BAM file"
        start_time = log_step_start(step_name)
        # 2-1. Extract read information from the BAM file
        df_reads = extract_reads_info(bam_path, chromosomes=chromosomes)
        # 2-2. Create a mapping dictionary for Read_Name and Read_Unique_Name
        read_mapping = create_read_mapping(df_reads)
        log_step_end(step_name, start_time)

        # 3. Classification process
        step_name = "Classification process"    
        start_time = log_step_start(step_name)
        # 3-1. Add overlap information between reads and variants
        df_reads, selected_read_unique_names = process_reads_and_variants(df_reads, union_variants)
        # 3-2. Add overlap information between variants and reads
        updated_union_variants, selected_variants = process_variants_and_reads(union_variants, df_reads)
        # 3-3. Call parallel processing function for classification
        ref_classifications, alt_classifications, missing_classifications, unknown_classifications = process_bam_data_parallel(
            bam_path,
            selected_read_unique_names,
            selected_variants,
            read_mapping,
            num_cores=num_cores,
            compute_missing_unknown=compute_missing_unknown
        )
        log_step_end(step_name, start_time)
        
        # 4. Save classification results
        step_name = "Save classification results"
        start_time = log_step_start(step_name)
        save_classification_matrices(
            save_dir=save_dir,
            union_variants=union_variants,
            barcode_path=barcode_path,
            ref_classifications=ref_classifications,
            alt_classifications=alt_classifications,
            missing_classifications=missing_classifications if compute_missing_unknown else None,
            unknown_classifications=unknown_classifications if compute_missing_unknown else None
        )
        log_step_end(step_name, start_time)
    
    except Exception as e:
        logging.error(f"An error occurred: {e}", exc_info=True)
        exit(1)
    
    # Log program completion
    logging.info("=== All steps completed successfully ===")

if __name__ == '__main__':
    main()