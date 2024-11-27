# variant_processing.py

import logging
import pandas as pd
import pysam

def is_vcf_file(file_path):
    """Check if the file is in VCF or BCF format"""
    return file_path.endswith('.vcf') or file_path.endswith('.vcf.gz') or file_path.endswith('.bcf')

def extract_variants_vcf(vcf_file, origin):
    """Extract and format variants from a VCF file"""
    formatted_variants = []

    try:
        with pysam.VariantFile(vcf_file) as vcf:
            for record in vcf:
                # Extract only 'PASS' filters
                if "PASS" in record.filter:
                    chrom = record.chrom
                    pos = record.pos
                    ref = record.ref
                    alts = ",".join(record.alts) if record.alts else []
                    qual = record.qual if record.qual is not None else "."
                    filter_ = ";".join(record.filter) if record.filter is not None else "."

                    variant = {
                        "CHROM": chrom,
                        "POS": pos,
                        "REF": ref,
                        "ALT": alts,
                        "QUAL": qual,
                        "FILTER": filter_,
                        "ORIGIN": origin  # Source of the variant
                    }
                    formatted_variants.append(variant)
    except ValueError as e:
        logging.warning(f"Skipping file {vcf_file} due to ValueError: {e}")
    except Exception as e:
        logging.warning(f"An unexpected error occurred while processing {vcf_file}: {e}")

    return formatted_variants

def extract_variants_txt(txt_file, origin):
    """
    Extract and format variants from a pre-filtered TXT file.
    """
    formatted_variants = []
    try:
        df = pd.read_csv(txt_file, sep='\t', header=None, names=["CHROM", "POS", "REF", "ALT"],
                         dtype={"CHROM": str, "POS": int, "REF": str, "ALT": str})
        
        # Check for the presence of required columns
        if df.shape[1] != 4:
            raise ValueError(f"The TXT file {txt_file} must have 4 columns.")
        
        for _, row in df.iterrows():
            chrom = row["CHROM"]
            pos = row["POS"]
            ref = row["REF"]
            alt = row["ALT"]
            variant = {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "QUAL": ".",        # Set QUAL to '.'
                "FILTER": ".",      # Set FILTER to '.'
                "ORIGIN": origin
            }
            formatted_variants.append(variant)
    except Exception as e:
        print(f"Warning: Skipping file {txt_file} due to error: {e}")
        logging.warning(f"Skipping file {txt_file} due to error: {e}")
    
    return formatted_variants

def filter_variants(variants):
    """Filter variants based on specific conditions."""
    union_variants_dict = {}
    
    for variant in variants:
        ref = variant["REF"]
        alt = variant["ALT"]
        key = (variant["CHROM"], variant["POS"], variant["REF"], variant["ALT"])

        # Exclude variants with commas or mismatched first characters
        if "," in ref or "," in alt:
            continue
        elif (len(ref) > 1 or len(alt) > 1) and ref[0] != alt[0]:
            continue
        else:
            # Update ORIGIN in the specified order
            if key not in union_variants_dict:
                union_variants_dict[key] = variant
            else:
                existing_origin = union_variants_dict[key]["ORIGIN"].split(", ")
                new_origin = variant["ORIGIN"]
                combined_origin = existing_origin + [new_origin]
                
                # Remove duplicates and maintain desired order
                origins_in_order = []
                for origin in ['normal', 'tumor', 'SCREADCOUNT']:
                    if origin in combined_origin and origin not in origins_in_order:
                        origins_in_order.append(origin)
                
                union_variants_dict[key]["ORIGIN"] = ", ".join(origins_in_order)
    
    # Convert to list and add VINDEX
    union_variants = list(union_variants_dict.values())

    for idx, variant in enumerate(union_variants, 1):
        variant['VINDEX'] = f'vi_{idx}'
    
    return union_variants

def process_vcf_files(vcf_files_with_origins):
    """
    Process multiple VCF/BCF or pre-filtered TXT files and return union_variants.
    """
    all_variants = []
    
    for file_path, origin in vcf_files_with_origins:
        if is_vcf_file(file_path):
            print(f"Processing VCF/BCF file: {file_path}")
            variants = extract_variants_vcf(file_path, origin)
        elif file_path.endswith('.txt'):
            print(f"Processing TXT file: {file_path}")
            variants = extract_variants_txt(file_path, origin)
        else:
            print(f"Warning: Skipping unsupported file format: {file_path}")
            logging.warning(f"Skipping unsupported file format: {file_path}")
            continue  # Skip unsupported file formats
        
        all_variants.extend(variants)
    
    union_variants = filter_variants(all_variants)
    
    return union_variants
