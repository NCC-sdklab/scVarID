#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
scVarID-like single-file code:
 - Window+padding chunk approach
 - classification.py (SNV/indel) logic integrated
 - Temporary pkl for each chunk
 - Final HDF5 output

Usage example:
 python scVarID_window_padding_all_in_one.py \
   --variant-files sample1 sample1.vcf.gz sample2 sample2.vcf.gz \
   --bam-path sample.bam \
   --save-dir results_out \
   --chromosomes chr1 \
   --window-size 10000000 \
   --padding 3000 \
   --num-cores 4
"""

import os
import sys
import time
import logging
import argparse
import re
from collections import defaultdict
import pysam
import h5py
import joblib
import pandas as pd
from scipy.sparse import csr_matrix
from joblib import Parallel, delayed

############################################################################
# 0. 유틸: Step 로그, 시간 포맷
############################################################################
def log_step_start(step_name):
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    end_time = time.time()
    elapsed = end_time - start_time
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {h:02d}h {m:02d}m {s:02d}s ===")

############################################################################
# 1. 인자 & 로깅
############################################################################
def parse_arguments():
    parser = argparse.ArgumentParser(description='scVarID all-in-one: window+padding + classification + tmp pkl.')
    parser.add_argument('--variant-files', nargs='+', required=True,
                        help='List of sample_name and variant file path pairs.')
    parser.add_argument('--bam-path', required=True, help='Path to the BAM file.')
    parser.add_argument('--save-dir', required=True, help='Output directory.')
    parser.add_argument('--barcode-path', required=False, help='(Optional) Barcode file.')
    parser.add_argument('--num-cores', type=int, default=4, help='Number of CPU cores.')
    parser.add_argument('--ref-alt-only', action='store_true',
                        help='Only compute ref/alt, skip missing/unknown.')
    parser.add_argument('--chromosomes', nargs='+', required=False,
                        help='Chromosomes to process.')
    parser.add_argument('--window-size', type=int, default=10_000_000,
                        help='Window size in bp (default=10Mb).')
    parser.add_argument('--padding', type=int, default=3000,
                        help='Padding around chunk boundary (default=3000).')

    args = parser.parse_args()

    if len(args.variant_files) % 2 != 0:
        parser.error('Please provide sample_name and variant_file pairs for --variant-files.')

    vcf_files_with_origins = []
    for i in range(0, len(args.variant_files), 2):
        sample = args.variant_files[i]
        fpath  = args.variant_files[i+1]
        vcf_files_with_origins.append((fpath, sample))

    return args, vcf_files_with_origins

def setup_logging(save_dir):
    os.makedirs(save_dir, exist_ok=True)
    log_file = os.path.join(save_dir, 'processing.log')
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s'
    )
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    fmt = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    console.setFormatter(fmt)
    logging.getLogger().addHandler(console)

def check_num_cores(n):
    import multiprocessing
    max_ = multiprocessing.cpu_count()
    if n>max_:
        logging.info(f"Requested {n} cores > available {max_}, use {max_}")
        return max_
    return n

############################################################################
# 2. Chrom/Windows
############################################################################
def get_chromosome_lengths(bam_path, target_chroms=None):
    with pysam.AlignmentFile(bam_path,"rb") as bam:
        all_chroms = list(bam.references)
        all_lengths = list(bam.lengths)
    chrom_len = dict(zip(all_chroms, all_lengths))
    if target_chroms:
        filtered = {}
        for c in target_chroms:
            if c in chrom_len:
                filtered[c] = chrom_len[c]
        return filtered
    else:
        return chrom_len

def get_windows_for_chromosome(chrom, chrom_length, window_size):
    start = 1
    while start <= chrom_length:
        end = min(start + window_size - 1, chrom_length)
        yield (start, end)
        start = end + 1

############################################################################
# 3. VCF fetch + filter (keys=Chromosome,POS,REF,ALT)
############################################################################
def is_vcf_file(fpath):
    return fpath.endswith('.vcf') or fpath.endswith('.vcf.gz') or fpath.endswith('.bcf')

def extract_variants_vcf_region(vcf_file, origin, chrom, start, end):
    vs = []
    try:
        with pysam.VariantFile(vcf_file) as vf:
            for rec in vf.fetch(chrom, start, end):
                if "PASS" in rec.filter:
                    v = {
                       "Chromosome": rec.chrom,
                       "POS": rec.pos,
                       "REF": rec.ref,
                       "ALT": ",".join(rec.alts) if rec.alts else "",
                       "QUAL": rec.qual if rec.qual else ".",
                       "FILTER": ";".join(rec.filter),
                       "ORIGIN": origin
                    }
                    vs.append(v)
    except Exception as e:
        logging.warning(f"Error reading {vcf_file} region={chrom}:{start}-{end}: {e}")
    return vs

def extract_variants_txt_region(txt_file, origin, chrom, start, end):
    import pandas as pd
    vs = []
    try:
        df = pd.read_csv(txt_file, sep='\t', header=None, names=["Chromosome","POS","REF","ALT"])
        sub = df[(df["Chromosome"]==chrom) & (df["POS"]>=start) & (df["POS"]<=end)]
        for _, row in sub.iterrows():
            v = {
               "Chromosome": row["Chromosome"],
               "POS": row["POS"],
               "REF": row["REF"],
               "ALT": row["ALT"],
               "QUAL": ".",
               "FILTER": ".",
               "ORIGIN": origin
            }
            vs.append(v)
    except Exception as e:
        logging.warning(f"Error reading {txt_file}: {e}")
    return vs

def filter_variants(variants_list):
    union_map = {}
    for v in variants_list:
        key = (v["Chromosome"], v["POS"], v["REF"], v["ALT"])
        if key not in union_map:
            union_map[key] = v
        else:
            oldo = union_map[key]["ORIGIN"].split(",")
            newo = v["ORIGIN"].split(",")
            union_map[key]["ORIGIN"] = ",".join(set(oldo+newo))
    unioned = list(union_map.values())
    for i,var in enumerate(unioned,1):
        var["VINDEX"] = f"vi_{i}"
    return unioned

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    all_vs = []
    for (fpath,origin) in vcf_files_with_origins:
        if is_vcf_file(fpath):
            chunk_vs = extract_variants_vcf_region(fpath, origin, chrom, start, end)
        elif fpath.endswith('.txt'):
            chunk_vs = extract_variants_txt_region(fpath, origin, chrom, start, end)
        else:
            logging.warning(f"Unsupported variant file: {fpath}")
            continue
        all_vs.extend(chunk_vs)
    union_vs = filter_variants(all_vs)
    return union_vs

############################################################################
# 4. Read fetch (padding) + Overlap => Classification
############################################################################
def extract_reads_with_padding(bam_path, chrom, start, end, padding=3000):
    data = []
    read_name_counts = {}
    fetch_start = max(1, start - padding)
    fetch_end   = end + padding

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, fetch_start-1, fetch_end):
            try:
                cb = read.get_tag("CB")
            except KeyError:
                cb = None
            if read.query_name not in read_name_counts:
                read_name_counts[read.query_name] = 0
            else:
                read_name_counts[read.query_name] += 1
            uid = read_name_counts[read.query_name]
            read_unique = f"{read.query_name}_{uid}"

            data.append([
                read.query_name,
                read_unique,
                cb,
                chrom,
                read.reference_start,
                read.reference_end
            ])
    df = pd.DataFrame(data, columns=[
        "Read_Name","Read_Unique_Name","Cell_Barcode","Chromosome",
        "Start_Position","End_Position"
    ])
    return df

def create_read_mapping(df_reads):
    """
    read_mapping: read_name -> list of read_unique_name
    """
    g = df_reads.groupby("Read_Name")["Read_Unique_Name"].apply(list)
    return g.to_dict()

def overlap_reads_variants(df_reads, union_variants, chunk_start, chunk_end):
    """
    Very simple: we just keep all reads; keep variants in [chunk_start, chunk_end]
    """
    final_vars = [v for v in union_variants if (v["POS"]>=chunk_start and v["POS"]<=chunk_end)]
    selected_read_uniques = df_reads["Read_Unique_Name"].tolist()
    return selected_read_uniques, final_vars

############################################################################
# 5. Classification.py 통합
############################################################################
def split_list(lst, num_sublists):
    avg_size = len(lst) // num_sublists
    remainder = len(lst) % num_sublists
    result = []
    start = 0
    for i in range(num_sublists):
        end = start + avg_size + (1 if i < remainder else 0)
        result.append(lst[start:end])
        start = end
    return result

def process_read_chunk(bam_file_path, read_unique_names_set, variants_dict, read_mapping, compute_missing_unknown=True):
    """
    chunk function from classification.py
    """
    import pysam
    ref_classifications = []
    alt_classifications = []
    missing_classifications = []
    unknown_classifications = []

    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        read_name_to_unique = {}
        for read in bam_file.fetch(until_eof=True):
            read_name = read.query_name
            if read_name not in read_name_to_unique:
                read_name_to_unique[read_name] = 0
            else:
                read_name_to_unique[read_name] += 1

            unique_id = read_name_to_unique[read_name]
            if read_name in read_mapping and unique_id < len(read_mapping[read_name]):
                read_unique_name = read_mapping[read_name][unique_id]
            else:
                read_unique_name = f"{read_name}_{unique_id}"

            if read_unique_name in read_unique_names_set:
                r_ref, r_alt, r_missing, r_unknown = process_read(
                    read, variants_dict, compute_missing_unknown
                )
                ref_classifications.extend(r_ref)
                alt_classifications.extend(r_alt)
                missing_classifications.extend(r_missing)
                unknown_classifications.extend(r_unknown)

                read_unique_names_set.remove(read_unique_name)
                if not read_unique_names_set:
                    break

    return ref_classifications, alt_classifications, missing_classifications, unknown_classifications

def process_bam_data_parallel(bam_path, selected_read_unique_names, selected_variants, read_mapping,
                              num_cores=4, compute_missing_unknown=True):
    variants_dict = create_variants_dict(selected_variants)
    read_unique_names_chunks = split_list(selected_read_unique_names, num_cores)

    ref_cls, alt_cls, missing_cls, unknown_cls = [], [], [], []

    results = Parallel(n_jobs=num_cores, backend='multiprocessing')(
        delayed(process_read_chunk)(
            bam_path,
            set(chunk),
            variants_dict,
            read_mapping,
            compute_missing_unknown
        )
        for chunk in read_unique_names_chunks
    )

    for r in results:
        ref_cls.extend(r[0])
        alt_cls.extend(r[1])
        missing_cls.extend(r[2])
        unknown_cls.extend(r[3])

    return ref_cls, alt_cls, missing_cls, unknown_cls

def process_read(read, variants_dict, compute_missing_unknown=True):
    """
    from classification.py: parse read => classify each variant position
    """
    read_name, cell_barcode, sequences, positions, operations = prime_read(read)
    chrom = read.reference_name

    read_start = read.reference_start + 1
    read_end   = read.reference_end

    ref_c = set()
    alt_c = set()
    missing_c = set()
    unknown_c = set()

    processed_positions = set()

    for pos in positions:
        key = (chrom, pos)
        if key in variants_dict and key not in processed_positions:
            processed_positions.add(key)
            # variants_dict[key] -> list of { 'variant':(chrom,pos,ref,alt), 'type':...}
            for vdata in variants_dict[key]:
                variant = vdata['variant']
                vtype   = vdata['type']
                # variant=(chrom,pos,ref,alt)
                if pos < read_start or pos>read_end:
                    continue
                classification, bc, var_id, rd_nm = classify_variant(
                    sequences, positions, operations, pos,
                    variant[2], variant[3], vtype,
                    operations[positions.index(pos)],
                    read_name, cell_barcode
                )
                cdict = {
                   'read_name': rd_nm,
                   'cell_barcode': bc,
                   'variant': f"{variant[0]}:{variant[1]}_{variant[2]}->{variant[3]}"
                }
                froz = frozenset(cdict.items())
                if classification=='ref':
                    ref_c.add(froz)
                elif classification=='alt':
                    alt_c.add(froz)
                elif classification=='missing' and compute_missing_unknown:
                    missing_c.add(froz)
                elif classification=='unknown' and compute_missing_unknown:
                    unknown_c.add(froz)

    return ref_c, alt_c, missing_c, unknown_c

def create_variants_dict(selected_variants):
    """
    from classification.py, expects each variant like:
     { 'Chromosome':..., 'POS':..., 'REF':..., 'ALT':...}
    """
    variants_dict = {}
    for var in selected_variants:
        chrom = var['Chromosome']
        pos   = int(var['POS'])
        ref   = var['REF']
        alt_list = var['ALT'].split(',')
        for alt in alt_list:
            if len(ref)>len(alt):
                vtype = 'deletion'
            elif len(ref)<len(alt):
                vtype = 'insertion'
            else:
                vtype = 'snv'
            key = (chrom, pos)
            if key not in variants_dict:
                variants_dict[key] = []
            variants_dict[key].append({'variant':(chrom,pos,ref,alt), 'type':vtype})
    return variants_dict

def classify_variant(sequences, positions, operations, pos, ref, alt, variant_type, operation, read_name, cell_barcode):
    if variant_type=='deletion':
        return handle_deletion(sequences, positions, operations, pos, ref, alt, read_name, cell_barcode)
    elif variant_type=='insertion':
        return handle_insertion(sequences, positions, operations, pos, ref, alt, read_name, cell_barcode)
    else:
        return handle_snv(sequences, positions, pos, ref, alt, read_name, cell_barcode)

def prime_read(read):
    read_name = read.query_name
    try:
        cell_barcode = read.get_tag("CB")
    except KeyError:
        cell_barcode = None

    sequences, positions, ops = [], [], []
    cur_pos = read.reference_start+1
    cur_readpos = 0
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', read.cigarstring)
    for length, op in cigar_ops:
        length = int(length)
        if op in "M=X":
            for _ in range(length):
                sequences.append(read.query_sequence[cur_readpos])
                positions.append(cur_pos)
                ops.append(op)
                cur_readpos += 1
                cur_pos += 1
        elif op=='I':
            for _ in range(length):
                sequences.append(read.query_sequence[cur_readpos])
                positions.append(-1)
                ops.append(op)
                cur_readpos += 1
        elif op in "DNP":
            for _ in range(length):
                sequences.append('-')
                positions.append(cur_pos)
                ops.append(op)
                cur_pos += 1
        elif op=='S':
            cur_readpos += length
        elif op=='H':
            pass
    return read_name, cell_barcode, sequences, positions, ops

def handle_snv(sequences, positions, pos, ref, alt, read_name, cell_barcode):
    for p, seq in zip(positions, sequences):
        if p==pos:
            if seq==ref:
                return 'ref', cell_barcode, pos, read_name
            elif seq==alt:
                return 'alt', cell_barcode, pos, read_name
            elif seq=='-':
                return 'missing', cell_barcode, pos, read_name
            else:
                return 'unknown', cell_barcode, pos, read_name
    return 'unknown', cell_barcode, pos, read_name

def handle_insertion(sequences, positions, ops, pos, ref, alt, read_name, cell_barcode):
    try:
        idx = positions.index(pos)
        if sequences[idx]=='-':
            return 'missing', cell_barcode, pos, read_name
        # 이 로직은 기본 예시. 실제 인서션 처리는 복잡함
        return 'unknown', cell_barcode, pos, read_name
    except:
        return 'unknown', cell_barcode, pos, read_name

def handle_deletion(sequences, positions, ops, pos, ref, alt, read_name, cell_barcode):
    try:
        idx = positions.index(pos)
        if sequences[idx]=='-':
            return 'missing', cell_barcode, pos, read_name
        return 'unknown', cell_barcode, pos, read_name
    except:
        return 'unknown', cell_barcode, pos, read_name

############################################################################
# 6. process_one_chunk_and_dump
############################################################################
def process_one_chunk_and_dump(
    chrom, win_start, win_end,
    vcf_files_with_origins,
    bam_path,
    padding=3000,
    compute_missing_unknown=True,
    num_cores=4,
    tmp_dir=None
):
    chunk_str = f"{chrom}:{win_start}-{win_end}"
    logging.info(f"[CHUNK] Start {chunk_str}")

    # variant fetch
    union_variants = process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(union_variants)} variants")
    if not union_variants:
        return None

    # read fetch
    df_reads = extract_reads_with_padding(bam_path, chrom, win_start, win_end, padding=padding)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(df_reads)} reads")
    if df_reads.empty:
        return None

    # overlap
    read_mapping = create_read_mapping(df_reads)
    selected_read_unique_names, final_variants = overlap_reads_variants(df_reads, union_variants, win_start, win_end)

    # classification.py 로직: process_bam_data_parallel
    ref_cls, alt_cls, miss_cls, unk_cls = process_bam_data_parallel(
        bam_path,
        selected_read_unique_names,
        final_variants,
        read_mapping,
        num_cores=num_cores,
        compute_missing_unknown=compute_missing_unknown
    )

    m_reads = len(selected_read_unique_names)
    t_reads = len(df_reads)
    p_reads = 100.0*m_reads/t_reads if t_reads>0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable reads: {m_reads} ({p_reads:.2f}%)")

    m_vars = len(final_variants)
    t_vars = len(union_variants)
    p_vars = 100.0*m_vars/t_vars if t_vars>0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable variants: {m_vars} ({p_vars:.2f}%)")

    chunk_data = {
        "union_variants": final_variants,
        "ref": ref_cls,
        "alt": alt_cls,
        "missing": miss_cls,
        "unknown": unk_cls
    }

    if not tmp_dir:
        tmp_dir = "./tmp_chunks"
    os.makedirs(tmp_dir, exist_ok=True)
    out_pkl = os.path.join(tmp_dir, f"chunk_{chrom}_{win_start}_{win_end}.pkl")
    joblib.dump(chunk_data, out_pkl)
    logging.info(f"[CHUNK] {chunk_str} => saved {out_pkl}")
    del chunk_data
    return out_pkl

############################################################################
# 7. Save final matrix
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
    # gather variants
    variants = []
    for v in union_variants:
        # (Chromosome, POS, REF, ALT)
        variants.append(f"{v['Chromosome']}:{v['POS']}_{v['REF']}->{v['ALT']}")
    variants = list(set(variants))  # 중복 제거 가능
    variants.sort()
    var_to_idx = {v:i for i,v in enumerate(variants)}

    # gather barcodes
    all_barcodes = set()
    def gather_barcodes(cls_):
        if not cls_:
            return
        for c in cls_:
            # c is frozenset => dict
            d = dict(c)
            all_barcodes.add(d["cell_barcode"])

    gather_barcodes(ref_classifications)
    gather_barcodes(alt_classifications)
    gather_barcodes(missing_classifications)
    gather_barcodes(unknown_classifications)
    if None in all_barcodes:
        all_barcodes.remove(None)  # remove no CB
    barcodes = sorted(all_barcodes)
    bc_to_idx = {b:i for i,b in enumerate(barcodes)}

    # build each classification
    def build_csr_and_save(cls_data, name):
        if not cls_data or len(cls_data)==0:
            return
        from collections import defaultdict
        data, row_inds, col_inds = [],[],[]
        dmap = defaultdict(int)
        for c_ in cls_data:
            dd = dict(c_)
            var_ = dd["variant"]
            bc_  = dd["cell_barcode"]
            if (bc_ is None) or (var_ not in var_to_idx) or (bc_ not in bc_to_idx):
                continue
            r = var_to_idx[var_]
            co = bc_to_idx[bc_]
            dmap[(r,co)] += 1

        for (r,co), val in dmap.items():
            row_inds.append(r)
            col_inds.append(co)
            data.append(val)

        mat = csr_matrix((data,(row_inds,col_inds)), shape=(len(variants), len(barcodes)))
        out_h5 = os.path.join(save_dir,f"{name}_matrix.h5")
        with h5py.File(out_h5,"w") as hf:
            hf.create_dataset("data", data=mat.data)
            hf.create_dataset("indices", data=mat.indices)
            hf.create_dataset("indptr", data=mat.indptr)
            hf.attrs['shape'] = mat.shape

    build_csr_and_save(ref_classifications,"ref")
    build_csr_and_save(alt_classifications,"alt")
    if missing_classifications:
        build_csr_and_save(missing_classifications,"missing")
    if unknown_classifications:
        build_csr_and_save(unknown_classifications,"unknown")

    # dump pkl
    joblib.dump((variants, barcodes), os.path.join(save_dir,"variant_barcode_mappings.pkl"))

############################################################################
# 8. main
############################################################################
def main():
    args, vcf_files_with_origins = parse_arguments()

    setup_logging(args.save_dir)
    logging.info("Command Line: " + " ".join(sys.argv))
    logging.info("=== [scVarID] Program Started ===")

    num_cores = check_num_cores(args.num_cores)
    compute_missing_unknown = not args.ref_alt_only

    step_t = log_step_start("Chromosome length retrieval")
    chrom_len = get_chromosome_lengths(args.bam_path, args.chromosomes)
    if not chrom_len:
        logging.error("No valid chrom found, abort.")
        sys.exit(1)
    log_step_end("Chromosome length retrieval", step_t)

    wsize = args.window_size
    pad   = args.padding
    logging.info(f"[INFO] window_size={wsize}, padding={pad}, num_cores={num_cores}")

    tmp_dir = os.path.join(args.save_dir, "tmp_chunks")
    os.makedirs(tmp_dir, exist_ok=True)

    # 1) generate chunk list
    step_t = log_step_start("Generate chunk list")
    all_chunks = []
    for c_, length_ in chrom_len.items():
        for (start_, end_) in get_windows_for_chromosome(c_, length_, wsize):
            all_chunks.append((c_,start_,end_))
    logging.info(f"Total chunk count: {len(all_chunks)}")
    log_step_end("Generate chunk list", step_t)

    # 2) parallel chunk => pkl
    step_t = log_step_start("Parallel chunk processing")
    chunk_files = Parallel(n_jobs=num_cores)(
      delayed(process_one_chunk_and_dump)(
        chrom=c[0],
        win_start=c[1],
        win_end=c[2],
        vcf_files_with_origins=vcf_files_with_origins,
        bam_path=args.bam_path,
        padding=pad,
        compute_missing_unknown=compute_missing_unknown,
        num_cores=num_cores,
        tmp_dir=tmp_dir
      )
      for c in all_chunks
    )
    log_step_end("Parallel chunk processing", step_t)

    # 3) merge
    step_t = log_step_start("Merge chunk results")
    final_variants = []
    ref_ = []
    alt_ = []
    missing_ = []
    unknown_ = []

    for f in chunk_files:
        if f is None:
            continue
        cdata = joblib.load(f)
        final_variants.extend(cdata["union_variants"])
        ref_.extend(cdata["ref"])
        alt_.extend(cdata["alt"])
        if compute_missing_unknown:
            missing_.extend(cdata["missing"])
            unknown_.extend(cdata["unknown"])
        del cdata
    log_step_end("Merge chunk results", step_t)

    # 4) save classification
    step_t = log_step_start("Save classification results")
    save_classification_matrices(
        save_dir=args.save_dir,
        union_variants=final_variants,
        barcode_path=args.barcode_path,
        ref_classifications=ref_,
        alt_classifications=alt_,
        missing_classifications=missing_ if compute_missing_unknown else None,
        unknown_classifications=unknown_ if compute_missing_unknown else None
    )
    log_step_end("Save classification results", step_t)

    logging.info("=== All steps completed successfully ===")

if __name__=="__main__":
    main()
