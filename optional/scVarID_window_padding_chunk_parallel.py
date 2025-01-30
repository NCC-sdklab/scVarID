#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Single-file scVarID pipeline:
 - Window + padding => chunk
 - **Chunk-level parallel** => no nested parallel warning
 - region-based BAM fetch => no until_eof
 - classification.py style: handle_snv/insertion/deletion
 - tmp/*.pkl => final HDF5

Usage example:
 python scVarID_window_padding_chunk_parallel.py \
   --variant-files sample1 sample1.vcf sample2 sample2.vcf \
   --bam-path sample.bam \
   --save-dir out_dir \
   --chromosomes chr1 \
   --window-size 10000000 \
   --padding 500 \
   --num-cores 30
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
# 0. Step log
############################################################################
def log_step_start(step_name):
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    end_t= time.time()
    elapsed= end_t - start_time
    hh = int(elapsed//3600)
    mm = int((elapsed%3600)//60)
    ss = int(elapsed%60)
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {hh:02d}h {mm:02d}m {ss:02d}s ===")

############################################################################
# 1. Parse args, logging
############################################################################
def parse_arguments():
    parser= argparse.ArgumentParser(description='scVarID pipeline: chunk-level parallel + region-based + single-file.')
    parser.add_argument('--variant-files', nargs='+', required=True,
                        help='Pairs of sample_name, file_path')
    parser.add_argument('--bam-path', required=True, help='BAM file path')
    parser.add_argument('--save-dir', required=True, help='Output directory')
    parser.add_argument('--barcode-path', required=False, help='(Optional) Barcode file')
    parser.add_argument('--num-cores', type=int, default=4, help='Chunk-level parallel cores')
    parser.add_argument('--ref-alt-only', action='store_true', help='Skip missing/unknown if set')
    parser.add_argument('--chromosomes', nargs='+', required=False, help='Chroms to process')
    parser.add_argument('--window-size', type=int, default=10_000_000, help='Chunk size')
    parser.add_argument('--padding', type=int, default=500, help='Padding around chunk')

    args= parser.parse_args()
    if len(args.variant_files)%2 !=0:
        parser.error('Must provide sample_name + vcf_file pairs in --variant-files.')
    vcf_files_with_origins=[]
    for i in range(0, len(args.variant_files), 2):
        sample= args.variant_files[i]
        f_ = args.variant_files[i+1]
        vcf_files_with_origins.append((f_, sample))
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
    import sys
    console= logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    cformatter= logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    console.setFormatter(cformatter)
    logging.getLogger().addHandler(console)

############################################################################
# 2. Chrom windows
############################################################################
def get_chromosome_lengths(bam_path, target_chroms=None):
    with pysam.AlignmentFile(bam_path,"rb") as bam:
        all_chroms= list(bam.references)
        all_lens= list(bam.lengths)
    cdict= dict(zip(all_chroms, all_lens))
    if target_chroms:
        filtered={}
        for c_ in target_chroms:
            if c_ in cdict:
                filtered[c_]= cdict[c_]
        return filtered
    else:
        return cdict

def get_windows_for_chromosome(chrom, length_, wsize):
    start=1
    while start<= length_:
        end= min(start+wsize-1, length_)
        yield (start,end)
        start= end+1

############################################################################
# 3. VCF region fetch
############################################################################
def is_vcf_file(fp):
    return fp.endswith('.vcf') or fp.endswith('.vcf.gz') or fp.endswith('.bcf')

def extract_variants_vcf_region(vcf_file, origin, chrom, start, end):
    vs=[]
    try:
        with pysam.VariantFile(vcf_file) as vf:
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
        logging.warning(f"Error reading {vcf_file}:{chrom}:{start}-{end} => {e}")
    return vs

def extract_variants_txt_region(txt_file, origin, chrom, start, end):
    vs=[]
    import pandas as pd
    try:
        df = pd.read_csv(txt_file, sep='\t', header=None, names=["Chromosome","POS","REF","ALT"])
        sub= df[(df["Chromosome"]==chrom)&(df["POS"]>=start)&(df["POS"]<=end)]
        for _, row in sub.iterrows():
            vs.append({
              "Chromosome": row["Chromosome"],
              "POS": row["POS"],
              "REF": row["REF"],
              "ALT": row["ALT"],
              "ORIGIN": origin
            })
    except Exception as e:
        logging.warning(f"Error reading {txt_file}: {e}")
    return vs

def filter_variants(vs):
    union_map={}
    for v_ in vs:
        key=(v_["Chromosome"], v_["POS"], v_["REF"], v_["ALT"])
        if key not in union_map:
            union_map[key]=v_
        else:
            oldo= union_map[key]["ORIGIN"].split(',')
            newo= v_["ORIGIN"].split(',')
            union_map[key]["ORIGIN"]= ",".join(set(oldo+newo))
    merged= list(union_map.values())
    for i,var in enumerate(merged,1):
        var["VINDEX"]=f"vi_{i}"
    return merged

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    allv=[]
    for (fp, origin) in vcf_files_with_origins:
        if is_vcf_file(fp):
            chunkv= extract_variants_vcf_region(fp, origin, chrom, start, end)
        elif fp.endswith('.txt'):
            chunkv= extract_variants_txt_region(fp, origin, chrom, start, end)
        else:
            chunkv=[]
            logging.warning(f"Unsupported variant file: {fp}")
        allv.extend(chunkv)
    unionv= filter_variants(allv)
    return unionv

############################################################################
# 4. region-based read fetch (with padding)
############################################################################
def extract_reads_with_padding(bam_path, chrom, start, end, padding=500):
    import pysam
    from collections import defaultdict
    data=[]
    read_name_counts= defaultdict(int)
    fs= max(1, start-padding)
    fe= end+padding
    with pysam.AlignmentFile(bam_path,"rb") as bam:
        for read in bam.fetch(chrom, fs-1, fe):
            rname= read.query_name
            try:
                bc= read.get_tag("CB")
            except KeyError:
                bc=None
            cnt= read_name_counts[rname]
            read_unique= f"{rname}_{cnt}"
            read_name_counts[rname]+=1

            data.append([
                rname, read_unique, bc,
                chrom, read.reference_start, read.reference_end
            ])
    df= pd.DataFrame(data, columns=[
      "Read_Name","Read_Unique_Name","Cell_Barcode","Chromosome",
      "Start_Pos","End_Pos"
    ])
    return df

############################################################################
# 5. classification.py style
############################################################################

def prime_read(read):
    rname= read.query_name
    try:
        bc_= read.get_tag("CB")
    except KeyError:
        bc_=None
    seqs=[]
    poss=[]
    ops=[]
    start_ = read.reference_start+1
    qpos=0
    import re
    cigar_ops= re.findall(r'(\d+)([MIDNSHP=X])', read.cigarstring)
    Q= read.query_sequence
    for (length_, op_) in cigar_ops:
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
    return rname, bc_, seqs, poss, ops

def handle_snv(sequences, positions, pos_, ref_, alt_, rname, bc_):
    for p_, s_ in zip(positions, sequences):
        if p_==pos_:
            if s_==ref_:
                return 'ref', bc_, pos_, rname
            elif s_==alt_:
                return 'alt', bc_, pos_, rname
            elif s_=='-':
                return 'missing', bc_, pos_, rname
            else:
                return 'unknown', bc_, pos_, rname
    return 'unknown', bc_, pos_, rname

def handle_insertion(seqs, poss, ops, p_, ref_, alt_, rname, bc_):
    try:
        idx= poss.index(p_)
        if seqs[idx]=='-':
            return 'missing', bc_, p_, rname
        return 'unknown', bc_, p_, rname
    except:
        return 'unknown', bc_, p_, rname

def handle_deletion(seqs, poss, ops, p_, ref_, alt_, rname, bc_):
    try:
        idx= poss.index(p_)
        if seqs[idx]=='-':
            return 'missing', bc_, p_, rname
        return 'unknown', bc_, p_, rname
    except:
        return 'unknown', bc_, p_, rname

def classify_variant(sequences, positions, operations, pos_, ref_, alt_, vtype, operation, rname, bc_):
    if vtype=='deletion':
        return handle_deletion(seqs=sequences, poss=positions, ops=operations, p_=pos_, ref_=ref_, alt_=alt_, rname=rname, bc_=bc_)
    elif vtype=='insertion':
        return handle_insertion(seqs=sequences, poss=positions, ops=operations, p_=pos_, ref_=ref_, alt_=alt_, rname=rname, bc_=bc_)
    else:
        return handle_snv(sequences, positions, pos_, ref_, alt_, rname, bc_)

def process_read(read, variants_dict, compute_missing_unknown=True):
    rname, bc_, seqs, poss, ops= prime_read(read)
    chrom= read.reference_name
    visited= set()

    refC=set()
    altC=set()
    missC=set()
    unkC=set()

    for p_ in poss:
        key=(chrom,p_)
        if key in variants_dict and key not in visited:
            visited.add(key)
            for v_ in variants_dict[key]:
                # v_ => { 'variant':(chrom,pos,ref,alt), 'type':?}
                var= v_['variant']
                vtype= v_['type']
                classification, bcX, posX, readN = classify_variant(
                    seqs, poss, ops, p_,
                    var[2], var[3], vtype,
                    ops[poss.index(p_)],
                    rname, bc_
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
    """
    sel_vars => list of {Chromosome, POS, REF, ALT...}
    """
    vd={}
    for v_ in sel_vars:
        chrom= v_['Chromosome']
        pos= int(v_['POS'])
        ref_= v_['REF']
        alt_list= v_['ALT'].split(',')
        key= (chrom, pos)
        if key not in vd:
            vd[key]=[]
        for alt_ in alt_list:
            if len(ref_)>len(alt_):
                t_='deletion'
            elif len(ref_)<len(alt_):
                t_='insertion'
            else:
                t_='snv'
            vd[key].append({'variant':(chrom,pos,ref_,alt_),'type':t_})
    return vd

############################################################################
# 6. Chunk classification (single-thread) => region fetch
############################################################################

def classify_chunk_reads_region(
    bam_path,
    chrom,
    start,
    end,
    read_uniq_list,
    variants_dict,
    read_mapping,
    compute_missing_unknown=True,
    padding=0
):
    """
    region-limited fetch => track read_name->count => read_uniq => classification
    """
    refC=[]
    altC=[]
    missC=[]
    unkC=[]

    read_uniq_set= set(read_uniq_list)
    read_name_counts= defaultdict(int)

    fs= max(1, start-padding)
    fe= end+padding

    with pysam.AlignmentFile(bam_path,"rb") as bam:
        for r_ in bam.fetch(chrom, fs-1, fe):
            rname= r_.query_name
            cnt= read_name_counts[rname]
            if rname in read_mapping:
                if cnt< len(read_mapping[rname]):
                    read_uniq= read_mapping[rname][cnt]
                else:
                    read_uniq= f"{rname}_{cnt}"
            else:
                read_uniq= f"{rname}_{cnt}"
            read_name_counts[rname]+=1

            if read_uniq in read_uniq_set:
                # do classification
                rr,ra,rm,ru= process_read(r_, variants_dict, compute_missing_unknown)
                refC.extend(rr)
                altC.extend(ra)
                missC.extend(rm)
                unkC.extend(ru)

                read_uniq_set.remove(read_uniq)
                if not read_uniq_set:
                    break
    return refC, altC, missC, unkC

############################################################################
# 7. process_one_chunk_and_dump => parallel
############################################################################

def process_one_chunk_and_dump(
    chrom, win_start, win_end,
    vcf_files_with_origins,
    bam_path,
    padding=500,
    compute_missing_unknown=True,
    tmp_dir=None
):
    """
    1) VCF region => union_variants
    2) BAM region => read_uniq
    3) build variant_dict
    4) classification => region fetch again?
    5) store pkl
    """
    chunk_str= f"{chrom}:{win_start}-{win_end}"
    logging.info(f"[CHUNK] Start {chunk_str}")

    # 1) VCF
    union_variants= process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(union_variants)} variants")
    if not union_variants:
        return None

    # 2) read fetch + mapping
    df_reads= extract_reads_with_padding(bam_path, chrom, win_start, win_end, padding=padding)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(df_reads)} reads")
    if df_reads.empty:
        return None

    read_mapping= {}
    from collections import defaultdict
    group_map = df_reads.groupby("Read_Name")["Read_Unique_Name"].apply(list).to_dict()
    read_mapping= dict(group_map)

    read_uniq_list= df_reads["Read_Unique_Name"].tolist()

    # 3) build variant dict
    # filter variants by chunk range
    final_variants= [v for v in union_variants if (v["POS"]>=win_start and v["POS"]<=win_end)]
    variants_dict= create_variants_dict(final_variants)

    # 4) classification => region-limited fetch (start..end + padding)
    refC, altC, missC, unkC= classify_chunk_reads_region(
        bam_path= bam_path,
        chrom= chrom,
        start= win_start,
        end= win_end,
        read_uniq_list= read_uniq_list,
        variants_dict= variants_dict,
        read_mapping= read_mapping,
        compute_missing_unknown= compute_missing_unknown,
        padding= padding
    )

    m_reads= len(read_uniq_list)
    t_reads= len(df_reads)
    p_reads= 100.0*m_reads/t_reads if t_reads>0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable reads: {m_reads} ({p_reads:.2f}%)")

    m_vars= len(final_variants)
    t_vars= len(union_variants)
    p_vars= 100.0*m_vars/t_vars if t_vars>0 else 0
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
    out_pkl= os.path.join(tmp_dir, f"chunk_{chrom}_{win_start}_{win_end}.pkl")
    joblib.dump(chunk_data, out_pkl)
    logging.info(f"[CHUNK] {chunk_str} => {out_pkl}")
    return out_pkl

############################################################################
# 8. Save final matrix
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
        key= f"{v_['Chromosome']}:{v_['POS']}_{v_['REF']}->{v_['ALT']}"
        variants.append(key)
    variants= list(set(variants))
    variants.sort()
    var_to_idx={v:i for i,v in enumerate(variants)}

    all_bc= set()
    def gather_bc(cls_):
        if cls_:
            for c_ in cls_:
                d_= dict(c_)
                bc_= d_["cell_barcode"]
                if bc_ is not None:
                    all_bc.add(bc_)
    gather_bc(ref_classifications)
    gather_bc(alt_classifications)
    gather_bc(missing_classifications)
    gather_bc(unknown_classifications)
    barcodes= sorted(all_bc)
    bc_to_idx= {b:i for i,b in enumerate(barcodes)}

    from collections import defaultdict

    def build_csr_and_save(cls_data, name):
        if not cls_data:
            return
        data=[]
        row_inds=[]
        col_inds=[]
        dmap= defaultdict(int)
        for c_ in cls_data:
            dd= dict(c_)
            var_= dd["variant"]
            bc_= dd["cell_barcode"]
            if var_ not in var_to_idx or bc_ not in bc_to_idx:
                continue
            r= var_to_idx[var_]
            co= bc_to_idx[bc_]
            dmap[(r,co)]+=1
        for (r,co),val in dmap.items():
            row_inds.append(r)
            col_inds.append(co)
            data.append(val)

        mat= csr_matrix((data,(row_inds,col_inds)), shape=(len(variants), len(barcodes)))
        out_h5= os.path.join(save_dir, f"{name}_matrix.h5")
        with h5py.File(out_h5,"w") as hf:
            hf.create_dataset("data", data=mat.data)
            hf.create_dataset("indices", data=mat.indices)
            hf.create_dataset("indptr", data=mat.indptr)
            hf.attrs['shape']= mat.shape

    build_csr_and_save(ref_classifications,"ref")
    build_csr_and_save(alt_classifications,"alt")
    if missing_classifications:
        build_csr_and_save(missing_classifications,"missing")
    if unknown_classifications:
        build_csr_and_save(unknown_classifications,"unknown")

    joblib.dump((variants, barcodes), os.path.join(save_dir, "variant_barcode_mappings.pkl"))

############################################################################
# 9. main => chunk-level parallel
############################################################################
def main():
    args, vcf_files_with_origins= parse_arguments()
    setup_logging(args.save_dir)
    logging.info("Command Line: " + " ".join(sys.argv))
    logging.info("=== [scVarID] Program Started ===")

    from joblib import Parallel, delayed

    compute_missing_unknown= not args.ref_alt_only

    step_t= log_step_start("Chrom length retrieval")
    chrom_len= get_chromosome_lengths(args.bam_path, args.chromosomes)
    if not chrom_len:
        logging.error("No valid chromosome => abort.")
        sys.exit(1)
    log_step_end("Chrom length retrieval", step_t)

    wsize= args.window_size
    pad  = args.padding
    num_cores= args.num_cores

    logging.info(f"[INFO] chunk-level parallel => n_jobs={num_cores}, window_size={wsize}, padding={pad}")

    # create chunk list
    step_t= log_step_start("Generate chunk list")
    all_chunks=[]
    for c_, length_ in chrom_len.items():
        start_=1
        while start_<=length_:
            e_= min(start_+ wsize -1, length_)
            all_chunks.append((c_, start_, e_))
            start_= e_+1
    logging.info(f"Total chunk count: {len(all_chunks)}")
    log_step_end("Generate chunk list", step_t)

    tmp_dir= os.path.join(args.save_dir,"tmp_chunks")
    os.makedirs(tmp_dir, exist_ok=True)

    # chunk-level parallel
    step_t= log_step_start("Parallel chunk processing")
    chunk_files= Parallel(n_jobs=num_cores)(
        delayed(process_one_chunk_and_dump)(
            chrom= c[0],
            win_start= c[1],
            win_end= c[2],
            vcf_files_with_origins= vcf_files_with_origins,
            bam_path= args.bam_path,
            padding= pad,
            compute_missing_unknown= compute_missing_unknown,
            tmp_dir= tmp_dir
        )
        for c in all_chunks
    )
    log_step_end("Parallel chunk processing", step_t)

    # merge
    step_t= log_step_start("Merge chunk results")
    final_variants=[]
    ref_=[]
    alt_=[]
    missing_=[]
    unknown_=[]
    for cf_ in chunk_files:
        if cf_ is None:
            continue
        cdata= joblib.load(cf_)
        final_variants.extend(cdata["union_variants"])
        ref_.extend(cdata["ref"])
        alt_.extend(cdata["alt"])
        if compute_missing_unknown:
            missing_.extend(cdata["missing"])
            unknown_.extend(cdata["unknown"])
        del cdata
    log_step_end("Merge chunk results", step_t)

    # save
    step_t= log_step_start("Save classification results")
    save_classification_matrices(
        save_dir= args.save_dir,
        union_variants= final_variants,
        barcode_path= args.barcode_path,
        ref_classifications= ref_,
        alt_classifications= alt_,
        missing_classifications= missing_ if compute_missing_unknown else None,
        unknown_classifications= unknown_ if compute_missing_unknown else None
    )
    log_step_end("Save classification results", step_t)

    logging.info("=== All steps completed successfully ===")


if __name__=="__main__":
    main()
