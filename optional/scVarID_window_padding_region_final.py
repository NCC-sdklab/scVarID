#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Single-file scVarID-like code (region-based, no nested parallel):
 - Window + padding -> chunk
 - For each chunk, we region-fetch from BAM (no until_eof) => partial read
 - track read_name -> nth occurrence => read_unique_name
 - classification.py style logic for SNV/indel
 - Save tmp chunk .pkl
 - Final merge -> HDF5 output
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

############################################################################
# 0. Step logs
############################################################################
def log_step_start(step_name):
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    end_t = time.time()
    elapsed = end_t - start_time
    hh = int(elapsed // 3600)
    mm = int((elapsed % 3600)//60)
    ss = int(elapsed %60)
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {hh:02d}h {mm:02d}m {ss:02d}s ===")

############################################################################
# 1. parse args & logging
############################################################################
def parse_arguments():
    parser = argparse.ArgumentParser(description='scVarID single-file: region-limited, no nested parallel.')
    parser.add_argument('--variant-files', nargs='+', required=True,
                        help='SAMPLE, VCF pairs.')
    parser.add_argument('--bam-path', required=True, help='BAM file path.')
    parser.add_argument('--save-dir', required=True, help='Output directory.')
    parser.add_argument('--barcode-path', required=False, help='(Optional) Barcode file.')
    parser.add_argument('--num-cores', type=int, default=4, help='Classification parallel cores.')
    parser.add_argument('--ref-alt-only', action='store_true', help='Skip missing/unknown.')
    parser.add_argument('--chromosomes', nargs='+', required=False, help='Chroms to process.')
    parser.add_argument('--window-size', type=int, default=10_000_000, help='Window size (bp).')
    parser.add_argument('--padding', type=int, default=5000, help='Padding (bp).')

    args = parser.parse_args()
    if len(args.variant_files)%2 !=0:
        parser.error('Must provide pairs: sample_name file_path.')
    vcf_files_with_origins=[]
    for i in range(0,len(args.variant_files),2):
        sample = args.variant_files[i]
        fpath  = args.variant_files[i+1]
        vcf_files_with_origins.append((fpath, sample))

    return args, vcf_files_with_origins

def setup_logging(save_dir):
    os.makedirs(save_dir, exist_ok=True)
    log_file = os.path.join(save_dir,'processing.log')
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s'
    )
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    cformatter= logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    console.setFormatter(cformatter)
    logging.getLogger().addHandler(console)

############################################################################
# 2. fetch chromosome windows
############################################################################
def get_chromosome_lengths(bam_path, target_chroms=None):
    with pysam.AlignmentFile(bam_path,"rb") as bam:
        all_chroms=list(bam.references)
        all_lens= list(bam.lengths)
    cdict=dict(zip(all_chroms, all_lens))
    if target_chroms:
        filtered={}
        for c in target_chroms:
            if c in cdict:
                filtered[c] = cdict[c]
        return filtered
    else:
        return cdict

def get_windows_for_chromosome(chrom, length_, wsize):
    start=1
    while start<=length_:
        end = min(start+wsize-1, length_)
        yield (start,end)
        start=end+1

############################################################################
# 3. Variant fetch & filter
############################################################################
def is_vcf_file(fp):
    return fp.endswith('.vcf') or fp.endswith('.vcf.gz') or fp.endswith('.bcf')

def extract_variants_vcf_region(f, origin, chrom, start, end):
    vs=[]
    try:
        with pysam.VariantFile(f) as vf:
            for rec in vf.fetch(chrom, start, end):
                if "PASS" in rec.filter:
                    vs.append({
                       "Chromosome": rec.chrom,
                       "POS": rec.pos,
                       "REF": rec.ref,
                       "ALT": ",".join(rec.alts) if rec.alts else "",
                       "ORIGIN": origin
                    })
    except Exception as e:
        logging.warning(f"Error reading {f}:{chrom}:{start}-{end} => {e}")
    return vs

def extract_variants_txt_region(f, origin, chrom, start, end):
    vs=[]
    import pandas as pd
    try:
        df = pd.read_csv(f, sep='\t', header=None, names=["Chromosome","POS","REF","ALT"])
        sub= df[(df["Chromosome"]==chrom) & (df["POS"]>=start)&(df["POS"]<=end)]
        for _, row in sub.iterrows():
            vs.append({
               "Chromosome": row["Chromosome"],
               "POS": row["POS"],
               "REF": row["REF"],
               "ALT": row["ALT"],
               "ORIGIN": origin
            })
    except Exception as e:
        logging.warning(f"Error reading {f}: {e}")
    return vs

def filter_variants(vs):
    union_map={}
    for v in vs:
        key=(v["Chromosome"], v["POS"], v["REF"], v["ALT"])
        if key not in union_map:
            union_map[key] = v
        else:
            oldo= union_map[key]["ORIGIN"].split(',')
            newo= v["ORIGIN"].split(',')
            union_map[key]["ORIGIN"] = ",".join(set(oldo+newo))
    merged= list(union_map.values())
    for i,var in enumerate(merged,1):
        var["VINDEX"] = f"vi_{i}"
    return merged

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    allv=[]
    for (fp, origin) in vcf_files_with_origins:
        if is_vcf_file(fp):
            chunkv=extract_variants_vcf_region(fp, origin, chrom, start, end)
        elif fp.endswith('.txt'):
            chunkv=extract_variants_txt_region(fp, origin, chrom, start, end)
        else:
            chunkv=[]
            logging.warning(f"Unsupported variant file: {fp}")
        allv.extend(chunkv)
    unionv=filter_variants(allv)
    return unionv

############################################################################
# 4. classification logic
############################################################################

def prime_read(read):
    r_name= read.query_name
    try:
        bc= read.get_tag("CB")
    except KeyError:
        bc=None
    seqs=[]
    poss=[]
    ops=[]
    start_ = read.reference_start+1
    qpos=0
    cigar_ops= re.findall(r'(\d+)([MIDNSHP=X])', read.cigarstring)
    Q= read.query_sequence
    for length_, op_ in cigar_ops:
        length_=int(length_)
        if op_ in "M=X":
            for _ in range(length_):
                seqs.append(Q[qpos])
                poss.append(start_)
                ops.append(op_)
                start_+=1
                qpos+=1
        elif op_=='I':
            for _ in range(length_):
                seqs.append(Q[qpos])
                poss.append(-1)
                ops.append(op_)
                qpos+=1
        elif op_ in "DNP":
            for _ in range(length_):
                seqs.append('-')
                poss.append(start_)
                ops.append(op_)
                start_+=1
        elif op_=='S':
            qpos+=length_
        elif op_=='H':
            pass
    return r_name, bc, seqs, poss, ops

def handle_snv(sequences, positions, pos_, ref_, alt_, read_name, bc_):
    for p_, s_ in zip(positions, sequences):
        if p_==pos_:
            if s_==ref_:
                return 'ref', bc_, pos_, read_name
            elif s_==alt_:
                return 'alt', bc_, pos_, read_name
            elif s_=='-':
                return 'missing', bc_, pos_, read_name
            else:
                return 'unknown', bc_, pos_, read_name
    return 'unknown', bc_, pos_, read_name

def handle_insertion(seqs, poss, ops, pos_, ref_, alt_, r_name, bc_):
    try:
        idx = poss.index(pos_)
        if seqs[idx]=='-':
            return 'missing', bc_, pos_, r_name
        return 'unknown', bc_, pos_, r_name
    except:
        return 'unknown', bc_, pos_, r_name

def handle_deletion(seqs, poss, ops, pos_, ref_, alt_, r_name, bc_):
    try:
        idx= poss.index(pos_)
        if seqs[idx]=='-':
            return 'missing', bc_, pos_, r_name
        return 'unknown', bc_, pos_, r_name
    except:
        return 'unknown', bc_, pos_, r_name

def classify_variant(sequences, positions, operations, pos_, ref_, alt_, variant_type, operation, r_name, bc_):
    if variant_type=='deletion':
        return handle_deletion(sequences, positions, operations, pos_, ref_, alt_, r_name, bc_)
    elif variant_type=='insertion':
        return handle_insertion(sequences, positions, operations, pos_, ref_, alt_, r_name, bc_)
    else:
        return handle_snv(sequences, positions, pos_, ref_, alt_, r_name, bc_)

def process_read(read, variants_dict, compute_missing_unknown=True):
    r_name, bc_, seqs, poss, ops = prime_read(read)
    chrom= read.reference_name
    read_start= read.reference_start+1
    read_end= read.reference_end
    refC=set()
    altC=set()
    missC=set()
    unkC=set()

    visited=set()
    for p_ in poss:
        key=(chrom,p_)
        if key in variants_dict and key not in visited:
            visited.add(key)
            for v_ in variants_dict[key]:
                # v_ => {'variant': (chrom,pos,ref,alt), 'type': ???}
                var= v_['variant']
                vtype= v_['type']
                classification, bcX, posX, readN = classify_variant(
                    seqs, poss, ops, p_,
                    var[2], var[3], vtype,
                    ops[poss.index(p_)],
                    r_name, bc_
                )
                froz= frozenset({
                  'read_name': readN,
                  'cell_barcode': bcX,
                  'variant': f"{var[0]}:{var[1]}_{var[2]}->{var[3]}"
                }.items())
                if classification=='ref':
                    refC.add(froz)
                elif classification=='alt':
                    altC.add(froz)
                elif classification=='missing' and compute_missing_unknown:
                    missC.add(froz)
                elif classification=='unknown' and compute_missing_unknown:
                    unkC.add(froz)

    return refC, altC, missC, unkC

def create_variants_dict(sel_vars):
    vd={}
    for v_ in sel_vars:
        chrom= v_['Chromosome']
        pos= int(v_['POS'])
        ref_= v_['REF']
        alt_list = v_['ALT'].split(',')
        key=(chrom,pos)
        if key not in vd:
            vd[key]=[]
        for alt_ in alt_list:
            if len(ref_)>len(alt_):
                vtype='deletion'
            elif len(ref_)<len(alt_):
                vtype='insertion'
            else:
                vtype='snv'
            vd[key].append({'variant':(chrom,pos,ref_,alt_),'type':vtype})
    return vd

############################################################################
# 5. chunk-limited fetch => read occurrence => classification
############################################################################

def process_read_chunk_in_region(
    bam_path,
    chunk_chrom,
    chunk_start,
    chunk_end,
    read_unique_names_set,
    variants_dict,
    read_mapping,
    compute_missing_unknown=True
):
    """
    region fetch => chunk_chrom, chunk_start, chunk_end
    track (read_name -> nth occurrence) => check read_uniq
    if in read_unique_names_set => classification
    """
    refC=[]
    altC=[]
    missC=[]
    unkC=[]

    read_name_counts= defaultdict(int)  # track nth occurrence per read_name
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chunk_chrom, chunk_start-1, chunk_end):
            r_name= read.query_name
            cnt= read_name_counts[r_name]
            # read_mapping[r_name] => list of read_unique
            if r_name in read_mapping:
                # if we still have an entry left:
                if cnt < len(read_mapping[r_name]):
                    read_uniq = read_mapping[r_name][cnt]
                else:
                    read_uniq = f"{r_name}_{cnt}"
            else:
                read_uniq = f"{r_name}_{cnt}"

            read_name_counts[r_name]+=1

            # check if read_uniq in read_unique_names_set
            if read_uniq in read_unique_names_set:
                # classification
                r_ref, r_alt, r_miss, r_unk = process_read(read, variants_dict, compute_missing_unknown)
                refC.extend(r_ref)
                altC.extend(r_alt)
                missC.extend(r_miss)
                unkC.extend(r_unk)

                read_unique_names_set.remove(read_uniq)
                if not read_unique_names_set:
                    break

    return refC, altC, missC, unkC

def process_bam_data_region(
    bam_path,
    chunk_chrom,
    chunk_start,
    chunk_end,
    selected_read_unique_names,
    variants_dict,
    read_mapping,
    compute_missing_unknown=True,
    num_cores=1
):
    """
    실제 분류 (단순히 직렬).
    (원하면 joblib Parallel로 분할 가능, but let's keep it simple.)
    """
    # 그냥 하나의 세션 => region fetch 한 번
    read_unique_names_set = set(selected_read_unique_names)
    refC, altC, missC, unkC = process_read_chunk_in_region(
        bam_path,
        chunk_chrom,
        chunk_start,
        chunk_end,
        read_unique_names_set,
        variants_dict,
        read_mapping,
        compute_missing_unknown
    )
    return refC, altC, missC, unkC

############################################################################
# 6. process_one_chunk_and_dump (메인에서 sequential call)
############################################################################
def process_one_chunk_and_dump(
    chrom, win_start, win_end,
    vcf_files_with_origins,
    bam_path,
    padding=5000,
    compute_missing_unknown=True,
    num_cores=1,
    tmp_dir=None
):
    chunk_str= f"{chrom}:{win_start}-{win_end}"
    logging.info(f"[CHUNK] Start {chunk_str}")

    # 1) variants
    union_variants = process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(union_variants)} variants")
    if not union_variants:
        return None

    # 2) read fetch + mapping
    df_reads = extract_reads_with_padding(bam_path, chrom, win_start, win_end, padding=padding)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(df_reads)} reads")
    if df_reads.empty:
        return None

    read_mapping = create_read_mapping(df_reads)
    read_uniq_list, final_variants = overlap_reads_variants(df_reads, union_variants, win_start, win_end)

    # 3) build variant dict
    variants_dict = create_variants_dict(final_variants)

    # 4) classification
    # region-based => fetch + read occurrence
    refC, altC, missC, unkC = process_bam_data_region(
        bam_path,
        chrom,
        win_start,
        win_end + padding,  # fetch end + padding to be safe
        read_uniq_list,
        variants_dict,
        read_mapping,
        compute_missing_unknown=compute_missing_unknown,
        num_cores=num_cores
    )

    # stats
    m_reads=len(read_uniq_list)
    t_reads=len(df_reads)
    p_reads=100.0*m_reads/t_reads if t_reads>0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable reads: {m_reads} ({p_reads:.2f}%)")

    m_vars=len(final_variants)
    t_vars=len(union_variants)
    p_vars=100.0*m_vars/t_vars if t_vars>0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable variants: {m_vars} ({p_vars:.2f}%)")

    chunk_data={
        "union_variants": final_variants,
        "ref": refC,
        "alt": altC,
        "missing": missC,
        "unknown": unkC
    }
    if not tmp_dir:
        tmp_dir="./tmp_chunks"
    os.makedirs(tmp_dir, exist_ok=True)
    out_pkl = os.path.join(tmp_dir, f"chunk_{chrom}_{win_start}_{win_end}.pkl")
    joblib.dump(chunk_data, out_pkl)
    logging.info(f"[CHUNK] {chunk_str} => {out_pkl}")
    return out_pkl

def extract_reads_with_padding(bam_path, chrom, start, end, padding=5000):
    """
    Fetch reads in region [start-padding, end+padding].
    Return a DataFrame with columns: 
      [Read_Name, Read_Unique_Name, Cell_Barcode, Chromosome, Start_Pos, End_Pos]
    (Here 'Read_Unique_Name' is formed by read_name + an internal counter)
    """
    import pysam
    import pandas as pd
    data = []
    read_name_counts = {}

    fetch_start = max(1, start - padding)
    fetch_end   = end + padding

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam.fetch(chrom, fetch_start-1, fetch_end):
            try:
                cb = r.get_tag("CB")
            except KeyError:
                cb = None
            if r.query_name not in read_name_counts:
                read_name_counts[r.query_name] = 0
            else:
                read_name_counts[r.query_name] += 1
            uid = read_name_counts[r.query_name]
            read_unique = f"{r.query_name}_{uid}"

            data.append([
                r.query_name,
                read_unique,
                cb,
                chrom,
                r.reference_start,
                r.reference_end
            ])
    df = pd.DataFrame(data, columns=[
        "Read_Name","Read_Unique_Name","Cell_Barcode","Chromosome","Start_Pos","End_Pos"
    ])
    return df


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
    variants=[]
    for v_ in union_variants:
        # {Chromosome,POS,REF,ALT}
        key=f"{v_['Chromosome']}:{v_['POS']}_{v_['REF']}->{v_['ALT']}"
        variants.append(key)
    variants=list(set(variants))
    variants.sort()
    var_to_idx={v:i for i,v in enumerate(variants)}

    # gather barcodes
    all_barcodes=set()
    def gather_barcodes(cls_):
        if cls_:
            for c_ in cls_:
                d_= dict(c_)
                bc_= d_["cell_barcode"]
                if bc_ is not None:
                    all_barcodes.add(bc_)

    gather_barcodes(ref_classifications)
    gather_barcodes(alt_classifications)
    gather_barcodes(missing_classifications)
    gather_barcodes(unknown_classifications)
    barcodes= sorted(all_barcodes)
    bc_to_idx={b:i for i,b in enumerate(barcodes)}

    def build_csr_and_save(cls_data, name):
        if not cls_data:
            return
        from collections import defaultdict
        data=[]
        row_inds=[]
        col_inds=[]
        dmap=defaultdict(int)
        for c_ in cls_data:
            dd= dict(c_)
            var_= dd["variant"]
            bc_ = dd["cell_barcode"]
            if var_ not in var_to_idx or bc_ not in bc_to_idx:
                continue
            r= var_to_idx[var_]
            co= bc_to_idx[bc_]
            dmap[(r,co)] +=1

        for (r,co), val in dmap.items():
            row_inds.append(r)
            col_inds.append(co)
            data.append(val)

        mat = csr_matrix((data,(row_inds,col_inds)), shape=(len(variants),len(barcodes)))
        out_h5= os.path.join(save_dir,f"{name}_matrix.h5")
        with h5py.File(out_h5,"w") as hf:
            hf.create_dataset("data", data=mat.data)
            hf.create_dataset("indices", data=mat.indices)
            hf.create_dataset("indptr", data=mat.indptr)
            hf.attrs['shape']=mat.shape

    build_csr_and_save(ref_classifications, "ref")
    build_csr_and_save(alt_classifications, "alt")
    if missing_classifications:
        build_csr_and_save(missing_classifications,"missing")
    if unknown_classifications:
        build_csr_and_save(unknown_classifications,"unknown")

    joblib.dump((variants, barcodes), os.path.join(save_dir,"variant_barcode_mappings.pkl"))

############################################################################
# 8. main
############################################################################
def main():
    args, vcf_files_with_origins = parse_arguments()

    setup_logging(args.save_dir)
    logging.info("Command Line: " + " ".join(sys.argv))
    logging.info("=== [scVarID] Program Started ===")

    num_cores= args.num_cores
    compute_missing_unknown = not args.ref_alt_only

    t0= log_step_start("Chromosome length retrieval")
    import pysam
    chrom_len = get_chromosome_lengths(args.bam_path, args.chromosomes)
    if not chrom_len:
        logging.error("No valid chrom found => abort")
        sys.exit(1)
    log_step_end("Chromosome length retrieval", t0)

    wsize= args.window_size
    pad= args.padding
    logging.info(f"[INFO] window_size={wsize}, padding={pad}, classification cores={num_cores}")

    tmp_dir= os.path.join(args.save_dir,"tmp_chunks")
    os.makedirs(tmp_dir, exist_ok=True)

    # generate chunk
    t0= log_step_start("Generate chunk list")
    all_chunks=[]
    for c_, l_ in chrom_len.items():
        start_=1
        while start_<=l_:
            e_= min(start_+ wsize -1, l_)
            all_chunks.append((c_,start_, e_))
            start_= e_+1
    logging.info(f"Total chunk count: {len(all_chunks)}")
    log_step_end("Generate chunk list", t0)

    # process chunk in sequential => no nested parallel
    t0=log_step_start("Process chunk (sequential region-based)")

    chunk_files=[]
    for (chrom_, st_, en_) in all_chunks:
        pkl_path= process_one_chunk_and_dump(
            chrom=chrom_,
            win_start=st_,
            win_end=en_,
            vcf_files_with_origins=vcf_files_with_origins,
            bam_path=args.bam_path,
            padding=pad,
            compute_missing_unknown=compute_missing_unknown,
            num_cores=num_cores,  # if you want internal parallel, see process_bam_data_region => you can adapt
            tmp_dir=tmp_dir
        )
        chunk_files.append(pkl_path)

    log_step_end("Process chunk (sequential region-based)", t0)

    # merge
    t0=log_step_start("Merge chunk results")
    final_variants=[]
    ref_=[]
    alt_=[]
    miss_=[]
    unk_=[]
    for f_ in chunk_files:
        if f_ is None:
            continue
        cdata= joblib.load(f_)
        final_variants.extend(cdata["union_variants"])
        ref_.extend(cdata["ref"])
        alt_.extend(cdata["alt"])
        if compute_missing_unknown:
            miss_.extend(cdata["missing"])
            unk_.extend(cdata["unknown"])
        del cdata
    log_step_end("Merge chunk results", t0)

    # save
    t0=log_step_start("Save classification results")
    save_classification_matrices(
        save_dir=args.save_dir,
        union_variants=final_variants,
        barcode_path=args.barcode_path,
        ref_classifications=ref_,
        alt_classifications=alt_,
        missing_classifications=miss_ if compute_missing_unknown else None,
        unknown_classifications=unk_ if compute_missing_unknown else None
    )
    log_step_end("Save classification results", t0)

    logging.info("=== All steps completed successfully ===")


if __name__=="__main__":
    main()
