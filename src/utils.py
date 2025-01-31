# utils.py

import os
import time
import multiprocessing
import logging
from logging.handlers import RotatingFileHandler
import h5py
from collections import defaultdict
import joblib
from scipy.sparse import csr_matrix
import argparse
import pysam

# def parse_arguments():
#     """
#     Parses command-line arguments and returns them along with variant file origins.
#     """
#     parser = argparse.ArgumentParser(description='Process variant files and BAM files.')
    
#     # Add argument for variant files
#     parser.add_argument('--variant-files', nargs='+', required=True, metavar=('ORIGIN', 'FILE_PATH'),
#                         help='List of origin and variant file path pairs.')
    
#     # Add arguments for BAM file, barcode file, and save directory
#     parser.add_argument('--bam-path', required=True, help='Path to the BAM file.')
#     parser.add_argument('--barcode-path', required=False, help='Path to the barcode file (optional).')
#     parser.add_argument('--save-dir', required=True, help='Save directory path.')
    
#     # Add argument for number of cores
#     parser.add_argument('--num-cores', type=int, default=4, help='Number of cores to use for parallel processing.')
    
#     # Add action
#     parser.add_argument('--ref-alt-only', action='store_true',
#                         help='Compute only ref and alt classifications, skipping missing and unknown.')
    
#     # Parse arguments
#     args = parser.parse_args()
    
#     # Generate list of variant files
#     vcf_files_with_origins = []
#     if len(args.variant_files) % 2 != 0:
#         parser.error('Each origin must be paired with a file path in --variant-files.')
    
#     for i in range(0, len(args.variant_files), 2):
#         origin = args.variant_files[i]
#         file_path = args.variant_files[i + 1]
#         vcf_files_with_origins.append((file_path, origin))
    
#     return args, vcf_files_with_origins
def parse_arguments():
    """
    Parses command-line arguments and returns them along with variant file origins.
    """
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
    
    # Add argument for chromosomes
    parser.add_argument('--chromosomes', nargs='+', required=False,
                        help='List of chromosomes to process (e.g., chr1 chr2). If not specified, process all.')

    # Add argument for window
    parser.add_argument('--window-size', type=int, default=10000000,
                       help='Processing window size in base pairs (default: 10,000,000)')
    parser.add_argument('--padding', type=int, default=5000000,
                       help='Padding size for window overlaps (default: 5,000,000)')

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
    log_file = os.path.join(save_dir, 'processing.log')
    
    # Use RotatingFileHandler to limit log file size and create backups
    handler = RotatingFileHandler(log_file, maxBytes=10**7, backupCount=5)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s',
        handlers=[
            handler,
            logging.StreamHandler()
        ]
    )

def log_step_start(step_name):
    """
    Logs the start of a step and returns the start time.
    """
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    """
    Logs the end of a step, calculates the elapsed time, and logs it.
    Formats the elapsed time as "00h 00m 00s".
    """
    end_time = time.time()
    duration = end_time - start_time
    formatted_duration = time.strftime("%Hh %Mm %Ss", time.gmtime(duration))
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {formatted_duration} ===\n")

def check_num_cores(num_cores):
    max_cores = multiprocessing.cpu_count()
    if num_cores > max_cores:
        logging.info(f"Warning: Specified number of cores ({num_cores}) exceeds available cores ({max_cores}). Using {max_cores} cores instead.")
        num_cores = max_cores
    return num_cores

def read_barcodes_from_file(barcode_path):
    """
    Reads barcodes from a file if barcode_path is provided, else returns None.
    """
    if barcode_path:
        with open(barcode_path, 'r') as f:
            barcodes = [line.strip() for line in f if line.strip()]
        logging.info(f"Total barcodes from file: {len(barcodes)}")
        return barcodes
    else:
        return None  # Return None if barcode file is not provided

# def save_classification_matrices(
#     save_dir, 
#     union_variants, 
#     barcode_path=None, 
#     ref_classifications=None, 
#     alt_classifications=None, 
#     missing_classifications=None, 
#     unknown_classifications=None
# ):
#     """
#     Save classified variant data.
#     """
#     # Extract variant list
#     variants = [f"{variant['CHROM']}:{variant['POS']}_{variant['REF']}->{variant['ALT']}" for variant in union_variants]
    
#     #logging.info(f"Total variants: {len(variants)}")
    
#     # Create mapping dictionaries
#     variant_to_index = {variant: idx for idx, variant in enumerate(variants)}
    
#     # Prepare classifications dictionary
#     classifications_dict = {}
#     if ref_classifications:
#         classifications_dict['ref'] = ref_classifications
#     if alt_classifications:
#         classifications_dict['alt'] = alt_classifications
#     if missing_classifications:
#         classifications_dict['missing'] = missing_classifications
#     if unknown_classifications:
#         classifications_dict['unknown'] = unknown_classifications

#     # Read barcodes from file or generate from data
#     barcodes = read_barcodes_from_file(barcode_path)
#     if barcodes is None:
#         barcode_set = set()
#         for classifications in classifications_dict.values():
#             for classification in classifications:
#                 classification_dict = dict(classification)
#                 barcode = classification_dict['cell_barcode']
#                 barcode_set.add(barcode)
#         barcodes = sorted(barcode_set)
#         logging.info(f"# Total barcodes generated from data: {len(barcodes)}")
    
#     barcode_to_index = {barcode: idx for idx, barcode in enumerate(barcodes)}
    
#     # Generate and save sparse matrices
#     classification_matrices = {}
    
#     for classification_name, classifications in classifications_dict.items():
#         if not classifications:
#             logging.info(f"No data for classification '{classification_name}'. Skipping.")
#             continue  # Skip if no data
        
#         data = []
#         row_indices = []
#         col_indices = []
#         count_dict = defaultdict(int)
    
#         for classification in classifications:
#             classification_dict = dict(classification)
#             variant = classification_dict['variant']
#             barcode = classification_dict['cell_barcode']
            
#             # Check if variant and barcode are in mapping dictionaries
#             if variant in variant_to_index and barcode in barcode_to_index:
#                 variant_idx = variant_to_index[variant]
#                 barcode_idx = barcode_to_index[barcode]
#                 key = (variant_idx, barcode_idx)
#                 count_dict[key] += 1  # Increment read count
#             else:
#                 #logging.warning(f"Variant '{variant}' or barcode '{barcode}' not found in mapping dictionaries.")
#                 continue
    
#         for (variant_idx, barcode_idx), count in count_dict.items():
#             row_indices.append(variant_idx)
#             col_indices.append(barcode_idx)
#             data.append(count)
    
#         # Create sparse matrix
#         matrix_shape = (len(variants), len(barcodes))
#         matrix = csr_matrix((data, (row_indices, col_indices)), shape=matrix_shape)
#         classification_matrices[classification_name] = matrix
    
#         # Save sparse matrix to h5 file
#         with h5py.File(f'{save_dir}/{classification_name}_matrix.h5', 'w') as h5f:
#             h5f.create_dataset('data', data=matrix.data)
#             h5f.create_dataset('indices', data=matrix.indices)
#             h5f.create_dataset('indptr', data=matrix.indptr)
#             h5f.attrs['shape'] = matrix.shape
    
#     # Save variant and barcode lists
#     joblib.dump((variants, barcodes), f'{save_dir}/variant_barcode_mappings.pkl')

# def save_classification_matrices(
#     save_dir, 
#     union_variants, 
#     barcode_path=None, 
#     ref_classifications=None, 
#     alt_classifications=None, 
#     missing_classifications=None, 
#     unknown_classifications=None,
#     chunk_size=100000  # 청크 크기 파라미터 추가
# ):
#     # 변수 초기화
#     variants = [f"{v['CHROM']}:{v['POS']}_{v['REF']}->{v['ALT']}" for v in union_variants]
#     variant_to_index = {v: i for i, v in enumerate(variants)}
    
#     # 바코드 처리
#     barcodes = read_barcodes_from_file(barcode_path) or sorted({
#         c['cell_barcode'] 
#         for cls in [ref_classifications, alt_classifications, 
#                    missing_classifications, unknown_classifications] 
#         if cls for c in cls
#     })
#     barcode_to_index = {b: i for i, b in enumerate(barcodes)}

#     # 분류별 처리
#     classifications = {
#         'ref': ref_classifications,
#         'alt': alt_classifications,
#         'missing': missing_classifications,
#         'unknown': unknown_classifications
#     }

#     for cls_name, cls_data in classifications.items():
#         if not cls_data:
#             continue

#         h5_path = f'{save_dir}/{cls_name}_matrix.h5'
#         first_chunk = True

#         # 청크 단위 처리
#         for i in range(0, len(cls_data), chunk_size):
#             chunk = cls_data[i:i+chunk_size]
            
#             data = []
#             rows = []
#             cols = []
#             count_dict = defaultdict(int)

#             # 현재 청크 처리
#             for entry in chunk:
#                 var = entry['variant']
#                 bcd = entry['cell_barcode']
#                 if var in variant_to_index and bcd in barcode_to_index:
#                     key = (variant_to_index[var], barcode_to_index[bcd])
#                     count_dict[key] += 1

#             # CSR 형식으로 변환
#             for (r, c), cnt in count_dict.items():
#                 data.append(cnt)
#                 rows.append(r)
#                 cols.append(c)

#             # HDF5 저장
#             with h5py.File(h5_path, 'a') as h5f:
#                 if first_chunk:
#                     # 첫 청크: 데이터셋 생성
#                     h5f.create_dataset('data', data=data, maxshape=(None,))
#                     h5f.create_dataset('indices', data=rows, maxshape=(None,))
#                     h5f.create_dataset('indptr', data=[0, len(cols)], maxshape=(None,))
#                     h5f.attrs['shape'] = (len(variants), len(barcodes))
#                     first_chunk = False
#                 else:
#                     # 추가 청크: 데이터 확장
#                     for ds_name, new_data in [('data', data), ('indices', rows)]:
#                         h5f[ds_name].resize((h5f[ds_name].shape[0] + len(new_data),))
#                         h5f[ds_name][-len(new_data):] = new_data
                    
#                     # indptr 업데이트
#                     h5f['indptr'].resize((h5f['indptr'].shape[0] + 1,))
#                     h5f['indptr'][-1] = h5f['indptr'][-2] + len(cols)

#     # 매핑 정보 저장
#     joblib.dump((variants, barcodes), f'{save_dir}/variant_barcode_mappings.pkl')

def save_classification_matrices(
    save_dir, 
    union_variants, 
    barcode_path=None, 
    ref_classifications=None, 
    alt_classifications=None, 
    missing_classifications=None, 
    unknown_classifications=None
):
    # Extract variant list
    variants = [f"{variant['CHROM']}:{variant['POS']}_{variant['REF']}->{variant['ALT']}" for variant in union_variants]
    
    # Create mapping dictionaries
    variant_to_index = {variant: idx for idx, variant in enumerate(variants)}
    
    # Prepare classifications dictionary
    classifications_dict = {}
    if ref_classifications:
        classifications_dict['ref'] = ref_classifications
    if alt_classifications:
        classifications_dict['alt'] = alt_classifications
    if missing_classifications:
        classifications_dict['missing'] = missing_classifications
    if unknown_classifications:
        classifications_dict['unknown'] = unknown_classifications

    # Read barcodes from file or generate from data
    barcodes = read_barcodes_from_file(barcode_path)
    if barcodes is None:
        barcode_set = set()
        for classifications in classifications_dict.values():
            for classification in classifications:
                classification_dict = dict(classification)
                barcode = classification_dict['cell_barcode']
                barcode_set.add(barcode)
        barcodes = sorted(barcode_set)
    logging.info(f"# Total barcodes: {len(barcodes)}")
    
    barcode_to_index = {barcode: idx for idx, barcode in enumerate(barcodes)}
    
    # Generate and save sparse matrices
    for classification_name, classifications in classifications_dict.items():
        if not classifications:
            logging.info(f"No data for classification '{classification_name}'. Skipping.")
            continue
        
        data = []
        row_indices = []
        col_indices = []
        count_dict = defaultdict(int)
    
        for classification in classifications:
            classification_dict = dict(classification)
            variant = classification_dict['variant']
            barcode = classification_dict['cell_barcode']
            
            if variant in variant_to_index and barcode in barcode_to_index:
                variant_idx = variant_to_index[variant]
                barcode_idx = barcode_to_index[barcode]
                key = (variant_idx, barcode_idx)
                count_dict[key] += 1
    
        for (variant_idx, barcode_idx), count in count_dict.items():
            row_indices.append(variant_idx)
            col_indices.append(barcode_idx)
            data.append(count)
    
        matrix_shape = (len(variants), len(barcodes))
        matrix = csr_matrix((data, (row_indices, col_indices)), shape=matrix_shape)
    
        with h5py.File(f'{save_dir}/{classification_name}_matrix.h5', 'w') as h5f:
            h5f.create_dataset('data', data=matrix.data)
            h5f.create_dataset('indices', data=matrix.indices)
            h5f.create_dataset('indptr', data=matrix.indptr)
            h5f.attrs['shape'] = matrix.shape
    
    # Save variant and barcode lists
    joblib.dump((variants, barcodes), f'{save_dir}/variant_barcode_mappings.pkl')


def get_all_chromosomes(bam_path, variant_paths=None):
    """
    Automatically extract a list of chromosomes from the BAM file or variant files.
    """
    # Extract chromosomes from BAM file
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        chromosomes = list(bam_file.references)

    # Validate chromosomes with variant files if provided
    if variant_paths:
        for variant_path in variant_paths:
            with pysam.VariantFile(variant_path) as vcf:
                file_chromosomes = list(vcf.header.contigs.keys())
                chromosomes = [chrom for chrom in chromosomes if chrom in file_chromosomes]

    return chromosomes
