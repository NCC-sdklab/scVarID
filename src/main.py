# main.py

import os
import logging
import argparse

from variant_processing import process_vcf_files
from read_processing import (
    extract_reads_info, 
    create_read_mapping, 
    process_reads_and_variants, 
    process_variants_and_reads
)
from classification import process_bam_data_parallel
from utils import (
    check_num_cores, 
    save_classification_matrices, 
    setup_logging, 
    log_step_start, 
    log_step_end
)

def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Process variant files and BAM files.')
    
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
    
    # Assign other arguments to variables
    bam_path = args.bam_path
    barcode_path = args.barcode_path
    save_dir = args.save_dir
    num_cores = args.num_cores
    
    # Set compute_missing_unknown value
    compute_missing_unknown = not args.ref_alt_only

    # Create save directory
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # Configure logging
    setup_logging(save_dir)
    
    # 0. Verify and adjust num_cores
    num_cores = check_num_cores(num_cores)
    
    # 1. Process variant files
    step_name = "Process variant files"
    start_time = log_step_start(step_name)
    union_variants = process_vcf_files(vcf_files_with_origins)
    log_step_end(step_name, start_time)
    
    # 2. Extract read information from the BAM file
    step_name = "Extract read information from the BAM file"
    start_time = log_step_start(step_name)
    df_reads = extract_reads_info(bam_path)
    log_step_end(step_name, start_time)
    
    # 3. Create a mapping dictionary for Read_Name and Read_Unique_Name
    step_name = "Process classification"
    start_time = log_step_start(step_name)
    read_mapping = create_read_mapping(df_reads)
    
    # 4. Add overlap information between reads and variants
    df_reads, selected_read_unique_names = process_reads_and_variants(df_reads, union_variants)
    
    # 5. Add overlap information between variants and reads
    updated_union_variants, selected_variants = process_variants_and_reads(union_variants, df_reads)
    
    # 6. Call parallel processing function for classification
    ref_classifications, alt_classifications, missing_classifications, unknown_classifications = process_bam_data_parallel(
        bam_path,
        selected_read_unique_names,
        selected_variants,
        read_mapping,
        num_cores=num_cores,
        compute_missing_unknown=compute_missing_unknown
    )
    log_step_end(step_name, start_time)
    
    # 7. Save classification results
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

if __name__ == '__main__':
    main()
