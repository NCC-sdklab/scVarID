#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
scVarID pipeline:
 - Window+padding => chunk-level parallel
 - PyRanges-based overlap => 'mappable reads', 'mappable variants'
 - handle_insertion/deletion for classification
 - final HDF5 output
Fixes KeyError: 'CHROM' by renaming 'Chromosome' -> 'CHROM' after PyRanges steps.
"""

import os
import sys
import time
import logging
import argparse
import re
from collections import defaultdict

import pysam
import pyranges as pr
import h5py
import joblib
import pandas as pd
from scipy.sparse import csr_matrix
from joblib import Parallel, delayed

############################################################################
# 0. Step logging
############################################################################

def log_step_start(step_name):
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    end_ = time.time()
    elapsed= end_ - start_time
    hh= int(elapsed//3600)
    mm= int((elapsed%3600)//60)
    ss= int(elapsed%60)
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {hh:02d}h {mm:02d}m {ss:02d}s ===")

############################################################################
# 1. Parse args + logging
############################################################################

def parse_arguments():
    parser= argparse.ArgumentParser(description='scVarID chunk parallel + mappable reads/variants + handle_insertion/deletion.')
    parser.add_argument('--variant-files', nargs='+', required=True,
                        help='Pairs: sampleName1 file1.vcf.gz sampleName2 file2.vcf.gz ...')
    parser.add_argument('--bam-path', required=True, help='Path to BAM file')
    parser.add_argument('--save-dir', required=True, help='Output directory')
    parser.add_argument('--barcode-path', required=False, help='(Optional) barcodes file')
    parser.add_argument('--num-cores', type=int, default=4, help='Parallel chunk processes')
    parser.add_argument('--ref-alt-only', action='store_true', help='Skip missing/unknown if set')
    parser.add_argument('--chromosomes', nargs='+', required=False, help='Which chromosomes to process')
    parser.add_argument('--window-size', type=int, default=10_000_000, help='Chunk size (bp)')
    parser.add_argument('--padding', type=int, default=500, help='Padding around chunk boundary (bp)')

    args= parser.parse_args()
    if len(args.variant_files)%2!=0:
        parser.error('Must provide sample_name + vcf_file pairs in --variant-files.')

    vcf_files_with_origins=[]
    for i in range(0, len(args.variant_files), 2):
        sample= args.variant_files[i]
        f_= args.variant_files[i+1]
        vcf_files_with_origins.append((f_, sample))

    return args, vcf_files_with_origins

def setup_logging(save_dir):
    os.makedirs(save_dir, exist_ok=True)
    log_file= os.path.join(save_dir,'processing.log')
    logging.basicConfig(
        filename= log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(message)s'
    )
    console= logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    cformatter= logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    console.setFormatter(cformatter)
    logging.getLogger().addHandler(console)

############################################################################
# 2. VCF region fetch + filter
############################################################################

def is_vcf_file(fpath):
    return fpath.endswith('.vcf') or fpath.endswith('.vcf.gz') or fpath.endswith('.bcf')

def extract_variants_vcf_region(vcf_file, origin, chrom, start, end):
    import pysam
    vs=[]
    try:
        with pysam.VariantFile(vcf_file) as vf:
            for rec in vf.fetch(chrom, start, end):
                if "PASS" in rec.filter:
                    c_= rec.chrom
                    p_= rec.pos
                    ref_= rec.ref
                    alts_= ",".join(rec.alts) if rec.alts else ""
                    var= {
                      "CHROM": c_,
                      "POS": p_,
                      "REF": ref_,
                      "ALT": alts_,
                      "QUAL": rec.qual if rec.qual else ".",
                      "FILTER": ";".join(rec.filter) if rec.filter else ".",
                      "ORIGIN": origin
                    }
                    vs.append(var)
    except Exception as e:
        logging.warning(f"extract_variants_vcf_region error in {vcf_file}: {e}")
    logging.info(f"Extracted {len(vs)} variants from region {chrom}:{start}-{end} in {vcf_file}")
    return vs

def extract_variants_txt_region(txt_file, origin, chrom, start, end):
    import pandas as pd
    vs=[]
    try:
        df= pd.read_csv(txt_file, sep='\t', header=None,
                        names=["CHROM","POS","REF","ALT"],
                        dtype={"CHROM":str,"POS":int,"REF":str,"ALT":str})
        if df.shape[1]!=4:
            raise ValueError(f"TXT file {txt_file} must have 4 columns.")
        sub= df[(df["CHROM"]==chrom)&(df["POS"]>=start)&(df["POS"]<=end)]
        for _, row in sub.iterrows():
            var= {
              "CHROM": row["CHROM"],
              "POS": row["POS"],
              "REF": row["REF"],
              "ALT": row["ALT"],
              "QUAL": ".",
              "FILTER": ".",
              "ORIGIN": origin
            }
            vs.append(var)
    except Exception as e:
        logging.warning(f"extract_variants_txt_region error in {txt_file}: {e}")
    logging.info(f"Extracted {len(vs)} variants from region {chrom}:{start}-{end} in {txt_file}")
    return vs

def filter_variants(variants):
    union_map= {}
    for v_ in variants:
        ref= v_["REF"]
        alt= v_["ALT"]
        key= (v_["CHROM"], v_["POS"], ref, alt)

        # exclude multi-commas
        if "," in ref or "," in alt:
            continue
        if (len(ref)>1 or len(alt)>1) and ref[0]!= alt[0]:
            continue

        if key not in union_map:
            union_map[key]= v_
        else:
            # merge ORIGIN
            exist_o= union_map[key]["ORIGIN"].split(", ")
            new_o= v_["ORIGIN"]
            combo= exist_o + [new_o]
            # keep order normal, tumor, SCREADCOUNT
            origin_ordered=[]
            for tag_ in ["normal","tumor","SCREADCOUNT"]:
                if tag_ in combo and tag_ not in origin_ordered:
                    origin_ordered.append(tag_)
            union_map[key]["ORIGIN"]= ", ".join(origin_ordered)

    unioned= list(union_map.values())
    for i,v_ in enumerate(unioned,1):
        v_["VINDEX"]= f"vi_{i}"
    logging.info(f"Filtered variants count: {len(unioned)}")
    return unioned

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    all_vs=[]
    for (fpath,origin) in vcf_files_with_origins:
        if is_vcf_file(fpath):
            vs= extract_variants_vcf_region(fpath, origin, chrom, start, end)
        elif fpath.endswith('.txt'):
            vs= extract_variants_txt_region(fpath, origin, chrom, start, end)
        else:
            vs=[]
            logging.warning(f"Skipping unsupported file: {fpath}")
        all_vs.extend(vs)
    union_= filter_variants(all_vs)
    return union_

############################################################################
# 3. Read fetch + PyRanges overlap => mappable reads/variants
############################################################################

def parse_cigar(cigar_string):
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)]

def calculate_end_position_cigar(read):
    ref_pos= read.reference_start
    ops= parse_cigar(read.cigarstring)
    for length,op_ in ops:
        if op_ in "M=XDN":
            ref_pos += length
    return ref_pos

def extract_reads_info_region(bam_path, chrom, start, end):
    data=[]
    from collections import defaultdict
    read_name_counts= defaultdict(int)
    with pysam.AlignmentFile(bam_path,"rb") as bam:
        for read in bam.fetch(chrom, start-1, end):
            try:
                bc= read.get_tag("CB")
            except KeyError:
                continue
            rname= read.query_name
            c_= read_name_counts[rname]
            read_unique= f"{rname}_{c_}"
            read_name_counts[rname]+=1

            read_end= calculate_end_position_cigar(read)
            data.append([
              rname, read_unique, bc, chrom,
              read.reference_start, read_end
            ])
    df= pd.DataFrame(data, columns=[
      "Read_Name","Read_Unique_Name","Cell_Barcode","Chromosome",
      "Start_Position","End_Position"
    ])
    logging.info(f"Extracted {len(df)} reads from region {chrom}:{start}-{end} in {bam_path}")
    return df

def annotate_reads_with_variants_pyranges(df_reads, union_variants):
    if not len(union_variants):
        df_reads["VINDEX"]=""
        return df_reads

    df_var= pd.DataFrame(union_variants).copy()
    def calc_end(row):
        ref_len= len(row["REF"])
        alt_len= len(row["ALT"])
        st= row["POS"]-1
        if ref_len> alt_len:
            return st+ ref_len
        elif ref_len< alt_len:
            return st+ alt_len
        else:
            return st+ ref_len

    df_var["Start"]= df_var["POS"]-1
    df_var["End"]= df_var.apply(calc_end, axis=1)
    # rename CHROM->Chromosome
    df_var.rename(columns={"CHROM":"Chromosome"}, inplace=True)

    var_pr= pr.PyRanges(df_var[["Chromosome","Start","End","VINDEX"]])

    tmp_rd= df_reads[["Read_Unique_Name","Chromosome","Start_Position","End_Position"]].copy()
    tmp_rd.rename(columns={"Start_Position":"Start","End_Position":"End"}, inplace=True)
    rd_pr= pr.PyRanges(tmp_rd)

    overlap= rd_pr.join(var_pr)
    odf= overlap.df
    # group by read => all VINDEX
    read_vindex= odf.groupby("Read_Unique_Name")["VINDEX"].apply(lambda x: ",".join(x)).reset_index()

    df_reads= df_reads.merge(read_vindex, on="Read_Unique_Name", how="left")
    df_reads["VINDEX"]= df_reads["VINDEX"].fillna("")

    total_= len(df_reads)
    mappable_= (df_reads["VINDEX"]!="").sum()
    perc= (100*mappable_/total_) if total_>0 else 0
    logging.info(f"# mappable reads: {mappable_} ({perc:.2f}%)")

    return df_reads

def annotate_variants_with_reads_pyranges(union_variants, df_reads):
    if not union_variants:
        return union_variants

    df_var= pd.DataFrame(union_variants).copy()
    def calc_end(row):
        ref_len= len(row["REF"])
        alt_len= len(row["ALT"])
        st= row["POS"]-1
        if ref_len> alt_len:
            return st+ ref_len
        elif ref_len< alt_len:
            return st+ alt_len
        else:
            return st+ ref_len

    df_var["Start"]= df_var["POS"]-1
    df_var["End"]= df_var.apply(calc_end, axis=1)
    df_var.rename(columns={"CHROM":"Chromosome"}, inplace=True)

    var_pr= pr.PyRanges(df_var[["Chromosome","Start","End","VINDEX"]])

    tmp_rd= df_reads[["Chromosome","Start_Position","End_Position","Read_Unique_Name"]].copy()
    tmp_rd.rename(columns={"Start_Position":"Start","End_Position":"End"}, inplace=True)
    rd_pr= pr.PyRanges(tmp_rd)

    overlap= var_pr.join(rd_pr)
    odf= overlap.df
    # group by variant => read list
    var_reads= odf.groupby("VINDEX")["Read_Unique_Name"].apply(lambda x: ",".join(x)).reset_index()

    df_var= df_var.merge(var_reads, on="VINDEX", how="left")
    df_var["Read_Names"]= df_var["Read_Unique_Name"].fillna("")
    df_var.drop(columns=["Start","End","Read_Unique_Name"], inplace=True)

    # **중요**: 여기서 'Chromosome' -> 'CHROM' 되돌려야
    df_var.rename(columns={"Chromosome":"CHROM"}, inplace=True)

    updated_vars= df_var.to_dict(orient="records")

    total_var= len(updated_vars)
    mappable_var= sum(1 for vv in updated_vars if vv["Read_Names"]!="")
    p_= (100*mappable_var/ total_var) if total_var>0 else 0
    logging.info(f"# mappable variants: {mappable_var} ({p_:.2f}%)")

    return updated_vars

def process_reads_and_variants_pyranges(df_reads, union_variants):
    df_reads2= annotate_reads_with_variants_pyranges(df_reads, union_variants)
    updated_vars= annotate_variants_with_reads_pyranges(union_variants, df_reads2)
    return df_reads2, updated_vars

############################################################################
# 4. Classification handle_insertion / handle_deletion
############################################################################

def prime_read_for_classify(read):
    rname= read.query_name
    try:
        bc= read.get_tag("CB")
    except KeyError:
        bc= None

    seqs=[]
    poss=[]
    ops=[]
    start_= read.reference_start+1
    qpos=0
    cigar_ops= re.findall(r'(\d+)([MIDNSHP=X])', read.cigarstring)
    Q= read.query_sequence

    for (length_, op_) in cigar_ops:
        length_= int(length_)
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
            qpos+= length_
        elif op_=='H':
            pass
    return rname, bc, seqs, poss, ops

def handle_snv(sequences, positions, pos_, ref_, alt_, rname, bc):
    for p_, s_ in zip(positions, sequences):
        if p_== pos_:
            if s_== ref_:
                return 'ref', bc, pos_, rname
            elif s_== alt_:
                return 'alt', bc, pos_, rname
            elif s_=='-':
                return 'missing', bc, pos_, rname
            else:
                return 'unknown', bc, pos_, rname
    return 'unknown', bc, pos_, rname

def handle_insertion(sequences, positions, ops, pos_, ref_, alt_, rname, bc):
    # minimal => fallback => unknown
    return 'unknown', bc, pos_, rname

def handle_deletion(sequences, positions, ops, pos_, ref_, alt_, rname, bc):
    # minimal => fallback => unknown
    return 'unknown', bc, pos_, rname

def classify_variant(seqs, poss, ops, pos_, ref_, alt_, variant_type, op_, rname, bc):
    if variant_type=='deletion':
        return handle_deletion(seqs, poss, ops, pos_, ref_, alt_, rname, bc)
    elif variant_type=='insertion':
        return handle_insertion(seqs, poss, ops, pos_, ref_, alt_, rname, bc)
    else:
        return handle_snv(seqs, poss, pos_, ref_, alt_, rname, bc)

def create_variants_dict(union_vars):
    vd={}
    for v_ in union_vars:
        # now we have v_["CHROM"], v_["POS"], v_["REF"], v_["ALT"]
        chrom= v_["CHROM"]
        pos= int(v_["POS"])
        ref_= v_["REF"]
        alt_list= v_["ALT"].split(',')
        key= (chrom, pos)
        if key not in vd:
            vd[key]= []
        for alt_ in alt_list:
            if len(ref_)> len(alt_):
                t_='deletion'
            elif len(ref_)< len(alt_):
                t_='insertion'
            else:
                t_='snv'
            vd[key].append({'variant':(chrom,pos,ref_,alt_), 'type':t_})
    return vd

def classify_chunk_reads_bam(
    bam_path,
    df_reads,
    updated_union_variants,
    compute_missing_unknown=True
):
    """
    single-thread classification
    read -> prime_read_for_classify -> classify each variant => ref/alt/missing/unknown
    """
    variants_dict= create_variants_dict(updated_union_variants)
    read_uniq_list= df_reads["Read_Unique_Name"].tolist()
    read_uniq_set= set(read_uniq_list)
    read_mapping= df_reads.groupby("Read_Name")["Read_Unique_Name"].apply(list).to_dict()

    refC=[]
    altC=[]
    missC=[]
    unkC=[]

    from collections import defaultdict
    read_name_counts= defaultdict(int)

    with pysam.AlignmentFile(bam_path,"rb") as bam:
        for read in bam.fetch(until_eof=True):
            rnm= read.query_name
            c_= read_name_counts[rnm]
            if rnm in read_mapping:
                if c_< len(read_mapping[rnm]):
                    read_uniq= read_mapping[rnm][c_]
                else:
                    read_uniq= f"{rnm}_{c_}"
            else:
                read_uniq= f"{rnm}_{c_}"

            read_name_counts[rnm]+=1

            if read_uniq in read_uniq_set:
                rname, bc, seqs, poss, ops= prime_read_for_classify(read)
                visited= set()
                chrom= read.reference_name
                for p_ in poss:
                    key= (chrom, p_)
                    if key in variants_dict and key not in visited:
                        visited.add(key)
                        for vv in variants_dict[key]:
                            var= vv["variant"] # (chr,pos,ref,alt)
                            vtype= vv["type"]
                            classification, bcX, pX, rdN= classify_variant(
                                seqs, poss, ops, p_,
                                var[2], var[3], vtype,
                                ops[poss.index(p_)],
                                rname, bc
                            )
                            froz= frozenset({
                              'read_name': rdN,
                              'cell_barcode': bcX,
                              'variant': f"{var[0]}:{var[1]}_{var[2]}->{var[3]}"
                            }.items())
                            if classification=='ref':
                                refC.append(froz)
                            elif classification=='alt':
                                altC.append(froz)
                            elif classification=='missing' and compute_missing_unknown:
                                missC.append(froz)
                            elif classification=='unknown' and compute_missing_unknown:
                                unkC.append(froz)

                read_uniq_set.remove(read_uniq)
                if not read_uniq_set:
                    break

    return refC, altC, missC, unkC

############################################################################
# 5. process_one_chunk_and_dump => chunk-level
############################################################################

def process_one_chunk_and_dump(
    chrom,
    win_start,
    win_end,
    vcf_files_with_origins,
    bam_path,
    padding=500,
    compute_missing_unknown=True,
    tmp_dir=None
):
    chunk_str= f"{chrom}:{win_start}-{win_end}"
    logging.info(f"[CHUNK] Start {chunk_str}")

    # 1) VCF => union_variants
    union_variants= process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(union_variants)} variants")
    if not union_variants:
        return None

    # 2) Reads => region-based
    df_reads= extract_reads_info_region(bam_path, chrom, win_start, win_end)
    if df_reads.empty:
        logging.info(f"[CHUNK] {chunk_str} no reads => skip.")
        return None
    logging.info(f"[CHUNK] {chunk_str} extracted {len(df_reads)} reads")

    # 3) PyRanges overlap => mappable reads => mappable variants
    df_reads2, updated_union_variants= process_reads_and_variants_pyranges(df_reads, union_variants)

    # 4) classification => handle_insertion/deletion => store
    refC, altC, missC, unkC= classify_chunk_reads_bam(
        bam_path,
        df_reads2,
        updated_union_variants,
        compute_missing_unknown= compute_missing_unknown
    )

    chunk_data= {
      "union_variants": updated_union_variants,
      "ref": refC,
      "alt": altC,
      "missing": missC,
      "unknown": unkC
    }
    if not tmp_dir:
        tmp_dir= "./tmp_chunks"
    os.makedirs(tmp_dir, exist_ok=True)

    out_pkl= os.path.join(tmp_dir, f"chunk_{chrom}_{win_start}_{win_end}.pkl")
    joblib.dump(chunk_data, out_pkl)
    logging.info(f"[CHUNK] {chunk_str} => {out_pkl}")

    return out_pkl

############################################################################
# 6. Save final matrix
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
    variants=[]
    for v_ in union_variants:
        # we have v_["CHROM"], v_["POS"], v_["REF"], v_["ALT"]
        key= f"{v_['CHROM']}:{v_['POS']}_{v_['REF']}->{v_['ALT']}"
        variants.append(key)
    variants= list(set(variants))
    variants.sort()
    var_to_idx= {v:i for i,v in enumerate(variants)}

    all_barcodes= set()
    def gather_barcodes(cls_):
        if cls_:
            for c_ in cls_:
                dd= dict(c_)
                bc_= dd["cell_barcode"]
                if bc_:
                    all_barcodes.add(bc_)

    gather_barcodes(ref_classifications)
    gather_barcodes(alt_classifications)
    gather_barcodes(missing_classifications)
    gather_barcodes(unknown_classifications)

    barcodes= sorted(all_barcodes)
    bc_to_idx= {b:i for i,b in enumerate(barcodes)}

    from collections import defaultdict

    def build_csr_and_save(cls_data, name):
        if not cls_data:
            return
        dmap= defaultdict(int)
        for c_ in cls_data:
            dd= dict(c_)
            var_= dd["variant"]
            bc_ = dd["cell_barcode"]
            if var_ not in var_to_idx or bc_ not in bc_to_idx:
                continue
            r= var_to_idx[var_]
            co= bc_to_idx[bc_]
            dmap[(r,co)] +=1

        data=[]
        row_inds=[]
        col_inds=[]
        for (r,co), val in dmap.items():
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

    joblib.dump((variants, barcodes), os.path.join(save_dir,"variant_barcode_mappings.pkl"))

############################################################################
# 7. main => chunk-level parallel
############################################################################

def main():
    args, vcf_files_with_origins= parse_arguments()
    setup_logging(args.save_dir)
    logging.info("Command Line: " + " ".join(sys.argv))
    logging.info("=== [scVarID] Program Started ===")

    compute_missing_unknown= not args.ref_alt_only
    num_cores= args.num_cores
    wsize= args.window_size
    pad  = args.padding

    # 1) get chrom lengths
    step_t= log_step_start("Chromosome length retrieval")
    with pysam.AlignmentFile(args.bam_path,"rb") as bam_:
        all_chroms= list(bam_.references)
        all_lens= list(bam_.lengths)
    chrom_len= dict(zip(all_chroms, all_lens))

    # filter if needed
    if args.chromosomes:
        filtered={}
        for c_ in args.chromosomes:
            if c_ in chrom_len:
                filtered[c_]= chrom_len[c_]
        chrom_len= filtered
    if not chrom_len:
        logging.error("No valid chromosome => abort.")
        sys.exit(1)
    log_step_end("Chromosome length retrieval", step_t)

    logging.info(f"[INFO] chunk-level parallel => n_jobs={num_cores}, window_size={wsize}, padding={pad}")

    # 2) build chunk list
    step_t= log_step_start("Generate chunk list")
    all_chunks=[]
    for c_, length_ in chrom_len.items():
        st=1
        while st<= length_:
            e_= min(st+ wsize-1, length_)
            all_chunks.append((c_, st, e_))
            st= e_+1
    logging.info(f"Total chunk count: {len(all_chunks)}")
    log_step_end("Generate chunk list", step_t)

    tmp_dir= os.path.join(args.save_dir,"tmp_chunks")
    os.makedirs(tmp_dir, exist_ok=True)

    # 3) parallel chunk
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

    # 4) merge
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

    # 5) save
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
