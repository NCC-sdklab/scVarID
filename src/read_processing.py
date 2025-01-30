# read_processing.py

import logging
import pandas as pd
import pysam
import pyranges as pr
import re

def parse_cigar(cigar):
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

def calculate_end_position_cigar(read):
    ref_pos = read.reference_start
    cigar_ops = parse_cigar(read.cigarstring)
    for length, op in cigar_ops:
        if op in "M=XDN":
            ref_pos += length
    return ref_pos

# def extract_reads_info(bam_file, chromosomes=None):
#     """
#     Extract read information from a BAM file.
#     """
#     data = []
#     read_name_counts = {}
    
#     with pysam.AlignmentFile(bam_file, "rb") as bam:
#         for read in bam.fetch():
#             try:
#                 barcode = read.get_tag("CB")  # Extract the "CB" tag for the cell barcode
#             except KeyError:
#                 continue  # Skip reads without a "CB" tag
            
#             chrom = bam.get_reference_name(read.reference_id)
            
#             # Skip reads that are not in the specified chromosomes
#             if chromosomes and chrom not in chromosomes:
#                 continue
            
#             start_pos = read.reference_start
#             end_pos = calculate_end_position_cigar(read)  # Helper function to compute end position
#             read_name = read.query_name
            
#             # Count occurrences of each read name
#             if read_name in read_name_counts:
#                 read_name_counts[read_name] += 1
#             else:
#                 read_name_counts[read_name] = 1
#             unique_id = read_name_counts[read_name] - 1  # Start numbering from 0
            
#             read_unique_name = f"{read_name}_{unique_id}"
#             data.append([read.query_name, read_unique_name, barcode, chrom, start_pos, end_pos])
    
#     # Convert to a DataFrame
#     df_reads = pd.DataFrame(data, columns=["Read_Name", "Read_Unique_Name", "Cell_Barcode", "Chromosome", "Start_Position", "End_Position"])
#     logging.info(f"Extracted {len(df_reads)} reads from BAM file: {bam_file}")
#     return df_reads

def extract_reads_info(bam_file, chromosomes=None, window_size=10_000_000, padding=5_000):
    data = []
    read_name_counts = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        target_chroms = chromosomes if chromosomes else bam.references
        
        for chrom in target_chroms:
            chrom_length = bam.get_reference_length(chrom)
            
            # Window 단위 처리
            for win_start in range(0, chrom_length, window_size):
                win_end = min(win_start + window_size + padding, chrom_length)
                
                # Read fetch (0-based)
                for read in bam.fetch(chrom, win_start, win_end):
                    try:
                        barcode = read.get_tag("CB")
                    except KeyError:
                        continue
                    
                    # Read 좌표 계산
                    read_start = read.reference_start
                    read_end = calculate_end_position_cigar(read)
                    
                    # 실제 window 범위 필터링 (padding 영역 제외)
                    #if read_end < win_start or read_start >= (win_start + window_size):
                    #    continue
                    
                    # Read 정보 저장
                    read_name = read.query_name
                    if read_name in read_name_counts:
                        read_name_counts[read_name] += 1
                    else:
                        read_name_counts[read_name] = 1
                    unique_id = read_name_counts[read_name] - 1
                    read_unique_name = f"{read_name}_{unique_id}"
                    
                    data.append([read_name, read_unique_name, barcode, chrom, read_start, read_end])
    
    df_reads = pd.DataFrame(data, columns=["Read_Name", "Read_Unique_Name", "Cell_Barcode", "Chromosome", "Start_Position", "End_Position"])
    logging.info(f"Extracted {len(df_reads)} reads from BAM file: {bam_file}")
    return df_reads    

def create_read_mapping(df_reads):
    """
    Create a mapping dictionary for Read_Name to Read_Unique_Name.
    """
    read_mapping = df_reads.groupby('Read_Name')['Read_Unique_Name'].apply(list).to_dict()
    return read_mapping

def annotate_reads_with_variants_pyranges(df_reads, union_variants):
    """
    Add a VINDEX column to df_reads using variant information from union_variants.
    """
    # Convert union_variants to a DataFrame
    variants_df = pd.DataFrame(union_variants)

    # Define a function to calculate the End position based on variant type
    def calculate_end(row):
        ref_len = len(row['REF'])
        alt_len = len(row['ALT'])
        start = row['POS'] - 1  # Convert to 0-based coordinate
        if ref_len > alt_len:
            # Deletion
            return start + ref_len
        elif ref_len < alt_len:
            # Insertion
            return start + alt_len  # Insertions occupy a single position
        else:
            # SNP or equal-length variants
            return start + ref_len

    # Calculate Start and End positions for variants
    variants_df['Start'] = variants_df['POS'] - 1
    variants_df['End'] = variants_df.apply(calculate_end, axis=1)

    # Rename columns to match PyRanges requirements
    variants_df.rename(columns={'CHROM': 'Chromosome'}, inplace=True)

    # Create PyRanges object for variants
    variants_pr = pr.PyRanges(variants_df[['Chromosome', 'Start', 'End', 'VINDEX']])

    # Select necessary columns from df_reads and rename for PyRanges
    reads_df = df_reads[['Read_Unique_Name', 'Chromosome', 'Start_Position', 'End_Position']].copy()
    reads_df.rename(columns={'Start_Position': 'Start', 'End_Position': 'End'}, inplace=True)

    # Create PyRanges object for reads
    reads_pr = pr.PyRanges(reads_df)

    # Find overlaps between reads and variants
    overlaps = reads_pr.join(variants_pr)

    # Convert overlaps to a DataFrame
    overlaps_df = overlaps.df

    # Group overlapping VINDEX values by Read_Unique_Name and concatenate them
    read_vindex = overlaps_df.groupby('Read_Unique_Name')['VINDEX'].apply(lambda x: ','.join(x)).reset_index()

    # Merge VINDEX data back into df_reads
    df_reads = df_reads.merge(read_vindex, on='Read_Unique_Name', how='left')

    # Replace NaN in VINDEX column with an empty string
    df_reads['VINDEX'] = df_reads['VINDEX'].fillna('')

    return df_reads

def process_reads_and_variants(df_reads, union_variants):
    """
    Process reads and variants, outputting statistics and mappable reads.
    """
    df_reads = annotate_reads_with_variants_pyranges(df_reads, union_variants)

    # Total number of reads
    total_rows = len(df_reads)

    # Calculate the number of rows where the VINDEX column is either an empty string '' or NaN
    num_empty_vindex = df_reads['VINDEX'].isnull().sum() + (df_reads['VINDEX'] == '').sum()

    # Calculate the number of rows where the VINDEX column is not empty
    num_non_empty_vindex = total_rows - num_empty_vindex

    # Output mappable read statistics
    percentage_mappable = (100 * num_non_empty_vindex) / total_rows if total_rows > 0 else 0
    logging.info(f"# mappable reads: {num_non_empty_vindex} ({percentage_mappable:.2f}%)")

    # Extract Read_Unique_Name with non-empty VINDEX
    selected_read_unique_names = df_reads[df_reads['VINDEX'] != '']['Read_Unique_Name'].tolist()

    return df_reads, selected_read_unique_names

def annotate_variants_with_reads(union_variants, df_reads):
    """
    Find overlapping Read_Unique_Name for each variant in union_variants and annotate them.
    """
    # Convert union_variants to a DataFrame
    variants_df = pd.DataFrame(union_variants)

    # Define a function to calculate the End position based on variant type
    def calculate_end(row):
        ref_len = len(row['REF'])
        alt_len = len(row['ALT'])
        start = row['POS'] - 1  # Convert to 0-based coordinate
        if ref_len > alt_len:
            # Deletion
            return start + ref_len
        elif ref_len < alt_len:
            # Insertion
            return start + alt_len  # Insertions occupy a single position
        else:
            # SNP or equal-length variants
            return start + ref_len

    # Calculate Start and End positions for variants
    variants_df['Start'] = variants_df['POS'] - 1
    variants_df['End'] = variants_df.apply(calculate_end, axis=1)

    # Rename columns to match PyRanges requirements
    variants_df.rename(columns={'CHROM': 'Chromosome'}, inplace=True)

    # Create PyRanges object for variants
    variants_pr = pr.PyRanges(variants_df[['Chromosome', 'Start', 'End', 'VINDEX']])

    # Select necessary columns from df_reads and rename for PyRanges
    reads_df = df_reads[['Chromosome', 'Start_Position', 'End_Position', 'Read_Unique_Name']].copy()
    reads_df.rename(columns={'Start_Position': 'Start', 'End_Position': 'End'}, inplace=True)

    # Create PyRanges object for reads
    reads_pr = pr.PyRanges(reads_df)

    # Find overlaps between variants and reads
    overlaps = variants_pr.join(reads_pr)

    # Convert overlaps to a DataFrame
    overlaps_df = overlaps.df

    # Group overlapping Read_Unique_Name values by VINDEX and concatenate them
    variant_readnames = overlaps_df.groupby('VINDEX')['Read_Unique_Name'].apply(lambda x: ','.join(x)).reset_index()

    # Merge the overlapping read names into the variants DataFrame
    variants_df = variants_df.merge(variant_readnames, on='VINDEX', how='left')

    # Replace NaN in Read_Unique_Name column with an empty string for variants with no overlapping reads
    variants_df['Read_Names'] = variants_df['Read_Unique_Name'].fillna('')

    # Drop unnecessary columns
    variants_df.drop(columns=['Start', 'End', 'Read_Unique_Name'], inplace=True)

    # Convert the resulting DataFrame back to a list of dictionaries
    updated_variants = variants_df.to_dict(orient='records')

    return updated_variants

def process_variants_and_reads(union_variants, df_reads):
    """
    Process variants and reads, outputting statistics and mappable variants.
    """
    updated_union_variants = annotate_variants_with_reads(union_variants, df_reads)
    
    # Total number of variants
    total_variants = len(updated_union_variants)
    
    # Calculate number of variants without Read_Names
    num_empty_readnames = sum(1 for variant in updated_union_variants if variant['Read_Names'] == '')
    
    # Calculate number of variants with Read_Names
    num_non_empty_readnames = total_variants - num_empty_readnames
    
    # Output mappable variant statistics
    percentage_mappable = (100 * num_non_empty_readnames) / total_variants if total_variants > 0 else 0
    logging.info(f"# mappable variants: {num_non_empty_readnames} ({percentage_mappable:.2f}%)")
    
    # Select variants with non-empty Read_Names
    selected_variants = [variant for variant in updated_union_variants if variant.get('Read_Names')]
    
    return updated_union_variants, selected_variants
