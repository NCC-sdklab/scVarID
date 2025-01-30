#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
scVarID pipeline:
 - Chunk-level parallel
 - region-based VCF/BAM fetch
 - parse_cigar + prime_read (from read_processing.py)
 - handle_insertion, handle_deletion, handle_snv (from classification.py)
 - tmp/*.pkl intermediate => final HDF5
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
# 0. Step logging
############################################################################
def log_step_start(step_name):
    logging.info(f"=== Step Start: {step_name} ===")
    return time.time()

def log_step_end(step_name, start_time):
    end_t= time.time()
    elapsed= end_t- start_time
    hh= int(elapsed//3600)
    mm= int((elapsed%3600)//60)
    ss= int(elapsed%60)
    logging.info(f"=== Step End: {step_name} | Elapsed Time: {hh:02d}h {mm:02d}m {ss:02d}s ===")

############################################################################
# 1. Arguments & logging
############################################################################
def parse_arguments():
    parser= argparse.ArgumentParser(description='scVarID chunk-level parallel + region-based + integrated classification logic.')
    parser.add_argument('--variant-files', nargs='+', required=True,
                        help='Pairs: SAMPLE, FILE (e.g. sample1 file1.vcf sample2 file2.vcf')
    parser.add_argument('--bam-path', required=True, help='BAM file path')
    parser.add_argument('--save-dir', required=True, help='Output directory')
    parser.add_argument('--barcode-path', required=False, help='(Optional) Barcode file')
    parser.add_argument('--num-cores', type=int, default=4, help='Number of chunk-level parallel processes')
    parser.add_argument('--ref-alt-only', action='store_true', help='If set, skip missing/unknown classification')
    parser.add_argument('--chromosomes', nargs='+', required=False, help='Chromosomes to process')
    parser.add_argument('--window-size', type=int, default=10_000_000, help='Chunk size (bp)')
    parser.add_argument('--padding', type=int, default=5000, help='Padding around chunk boundary (bp)')
    args= parser.parse_args()

    if len(args.variant_files)%2 !=0:
        parser.error("Must provide sample_name + vcf_file pairs in --variant-files.")

    vcf_files_with_origins=[]
    for i in range(0, len(args.variant_files), 2):
        sample= args.variant_files[i]
        f_    = args.variant_files[i+1]
        vcf_files_with_origins.append((f_, sample))

    return args, vcf_files_with_origins

def setup_logging(save_dir):
    os.makedirs(save_dir, exist_ok=True)
    log_file= os.path.join(save_dir, 'processing.log')
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
# 2. Variant processing
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
                    qual_= rec.qual if rec.qual else "."
                    filt_= ";".join(rec.filter) if rec.filter else "."
                    var= {
                      "CHROM": c_,
                      "POS": p_,
                      "REF": ref_,
                      "ALT": alts_,
                      "QUAL": qual_,
                      "FILTER": filt_,
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
            raise ValueError(f"The TXT file {txt_file} must have 4 columns.")
        sub= df[(df["CHROM"]== chrom)&(df["POS"]>= start)&(df["POS"]<= end)]
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

        if "," in ref or "," in alt:
            continue
        if (len(ref)>1 or len(alt)>1) and ref[0]!= alt[0]:
            continue

        if key not in union_map:
            union_map[key]= v_
        else:
            # merge ORIGIN
            existing_origin= union_map[key]["ORIGIN"].split(", ")
            new_o= v_["ORIGIN"]
            combo= existing_origin+ [new_o]
            origins_in_order=[]
            for or_ in ["normal","tumor","SCREADCOUNT"]:
                if or_ in combo and or_ not in origins_in_order:
                    origins_in_order.append(or_)
            union_map[key]["ORIGIN"]= ", ".join(origins_in_order)

    unioned= list(union_map.values())
    for i, v_ in enumerate(unioned, 1):
        v_["VINDEX"]= f"vi_{i}"
    logging.info(f"Filtered variants count: {len(unioned)}")
    return unioned

def process_vcf_files_region(vcf_files_with_origins, chrom, start, end):
    all_vs=[]
    for (fpath, origin) in vcf_files_with_origins:
        if is_vcf_file(fpath):
            v_ = extract_variants_vcf_region(fpath, origin, chrom, start, end)
        elif fpath.endswith('.txt'):
            v_ = extract_variants_txt_region(fpath, origin, chrom, start, end)
        else:
            v_=[]
            logging.warning(f"Skipping unsupported variant file: {fpath}")
        all_vs.extend(v_)
    union_vs= filter_variants(all_vs)
    return union_vs

############################################################################
# 3. Read fetch + classification logic
############################################################################

def extract_reads_with_padding(bam_path, chrom, start, end, padding=5000):
    data=[]
    read_name_counts= defaultdict(int)
    fs= max(1, start-padding)
    fe= end+ padding

    with pysam.AlignmentFile(bam_path,"rb") as bam:
        for read in bam.fetch(chrom, fs-1, fe):
            rname= read.query_name
            try:
                bc= read.get_tag("CB")
            except KeyError:
                bc=None
            c_= read_name_counts[rname]
            read_unique= f"{rname}_{c_}"
            read_name_counts[rname]+=1

            data.append([
                rname, read_unique, bc, chrom,
                read.reference_start, read.reference_end
            ])
    df= pd.DataFrame(data, columns=[
      "Read_Name","Read_Unique_Name","Cell_Barcode","Chromosome",
      "Start_Pos","End_Pos"
    ])
    logging.info(f"Extracted {len(df)} reads from region {chrom}:{start}-{end} in {bam_path}")
    return df

########## from read_processing.py (parse_cigar, prime_read) ##########
def parse_cigar(cigar):
    # parse strings like "20M2I5D..." => list of (length, op)
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

def prime_read(read):
    """
    parse read => (read_name, cell_barcode, sequences, positions, ops)
    using parse_cigar
    """
    read_name= read.query_name
    try:
        bc= read.get_tag("CB")
    except KeyError:
        bc= None

    sequences=[]
    positions=[]
    ops=[]
    current_pos= read.reference_start+1
    qpos=0
    cigar_ops= parse_cigar(read.cigarstring)
    Q= read.query_sequence

    for (length_, op_) in cigar_ops:
        if op_ in "M=X":
            for _ in range(length_):
                sequences.append(Q[qpos])
                positions.append(current_pos)
                ops.append(op_)
                current_pos+=1
                qpos+=1
        elif op_=='I':
            for _ in range(length_):
                sequences.append(Q[qpos])
                positions.append(-1)  # insertion pos
                ops.append(op_)
                qpos+=1
        elif op_ in "DNP":
            for _ in range(length_):
                sequences.append('-')
                positions.append(current_pos)
                ops.append(op_)
                current_pos+=1
        elif op_=='S':
            # softclip => skip these bases
            qpos+= length_
        elif op_=='H':
            # hardclip => do nothing
            pass

    return read_name, bc, sequences, positions, ops


########## from classification.py (handle_insertion, etc.) ##########

def handle_snv(sequences, positions, pos_, ref_, alt_, read_name, cell_barcode):
    for p_, s_ in zip(positions, sequences):
        if p_== pos_:
            if s_== ref_:
                return 'ref', cell_barcode, pos_, read_name
            elif s_== alt_:
                return 'alt', cell_barcode, pos_, read_name
            elif s_=='-':
                return 'missing', cell_barcode, pos_, read_name
            else:
                return 'unknown', cell_barcode, pos_, read_name
    return 'unknown', cell_barcode, pos_, read_name

def handle_insertion(sequences, positions, operations, pos_, ref_, alt_, read_name, cell_barcode):
    """
    Simple insertion logic:
    if length(alt) > length(ref) => possibly insertion
    fallback => unknown
    """
    try:
        pos_index= positions.index(pos_)

        # if read base='-' => missing
        if sequences[pos_index]=='-':
            return 'missing', cell_barcode, pos_, read_name

        # if we want to do minimal checking => same as old
        # or if length is 1 => check quickly...
        # (More advanced logic can be placed here if needed)
    except ValueError:
        return 'unknown', cell_barcode, pos_, read_name
    except IndexError:
        return 'unknown', cell_barcode, pos_, read_name

    # fallback => unknown
    return 'unknown', cell_barcode, pos_, read_name

def handle_deletion(sequences, positions, operations, pos_, ref_, alt_, read_name, cell_barcode):
    """
    Simple deletion logic
    """
    try:
        pos_index= positions.index(pos_)
        if sequences[pos_index]=='-':
            return 'missing', cell_barcode, pos_, read_name
    except ValueError:
        return 'unknown', cell_barcode, pos_, read_name
    except IndexError:
        return 'unknown', cell_barcode, pos_, read_name

    # fallback => unknown
    return 'unknown', cell_barcode, pos_, read_name

def classify_variant(seqs, poss, ops, pos_, ref_, alt_, variant_type, operation, rname, bc_):
    if variant_type=='deletion':
        return handle_deletion(seqs, poss, ops, pos_, ref_, alt_, rname, bc_)
    elif variant_type=='insertion':
        return handle_insertion(seqs, poss, ops, pos_, ref_, alt_, rname, bc_)
    else:
        return handle_snv(seqs, poss, pos_, ref_, alt_, rname, bc_)

def process_read(read, variants_dict, compute_missing_unknown=True):
    """
    parse read => for each variant position => classify => (ref,alt,missing,unknown)
    """
    rname, bc_, seqs, poss, ops= prime_read(read)
    chrom= read.reference_name

    visited= set()
    refC=set()
    altC=set()
    missC=set()
    unkC=set()

    for p_ in poss:
        key= (chrom,p_)
        if key in variants_dict and key not in visited:
            visited.add(key)
            for v_ in variants_dict[key]:
                var= v_['variant']  # (chrom,pos,ref,alt)
                vtype= v_['type']
                classification, bcX, posX, readN= classify_variant(
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

def create_variants_dict(selected_variants):
    vd={}
    for var in selected_variants:
        chrom= var['CHROM']
        pos= int(var['POS'])
        ref_= var['REF']
        alt_str= var['ALT']
        alt_list= alt_str.split(',')
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
            vd[key].append({'variant': (chrom,pos,ref_,alt_), 'type':t_})
    return vd


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
    region-limited fetch => parse read => if read_uniq in read_uniq_list => classification
    """
    refC=[]
    altC=[]
    missC=[]
    unkC=[]
    read_uniq_set= set(read_uniq_list)
    read_name_counts= defaultdict(int)

    fs= max(1, start-padding)
    fe= end+ padding

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
# 4. process_one_chunk_and_dump => chunk-level parallel function
############################################################################

def process_one_chunk_and_dump(
    chrom, win_start, win_end,
    vcf_files_with_origins,
    bam_path,
    padding=5000,
    compute_missing_unknown=True,
    tmp_dir=None
):
    chunk_str= f"{chrom}:{win_start}-{win_end}"
    logging.info(f"[CHUNK] Start {chunk_str}")

    # 1) VCF region => union_variants
    union_variants= process_vcf_files_region(vcf_files_with_origins, chrom, win_start, win_end)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(union_variants)} variants")
    if not union_variants:
        return None

    # 2) read fetch => read_mapping
    df_reads= extract_reads_with_padding(bam_path, chrom, win_start, win_end, padding=padding)
    logging.info(f"[CHUNK] {chunk_str} extracted {len(df_reads)} reads")
    if df_reads.empty:
        return None

    # build read_mapping
    group_map= df_reads.groupby("Read_Name")["Read_Unique_Name"].apply(list).to_dict()
    read_mapping= dict(group_map)
    read_uniq_list= df_reads["Read_Unique_Name"].tolist()

    # 3) filter variants by chunk range
    final_variants= [v for v in union_variants if (v["POS"]>= win_start and v["POS"]<= win_end)]

    # create variants_dict
    variants_dict= create_variants_dict(final_variants)

    # 4) classification => region-limited
    refC, altC, missC, unkC= classify_chunk_reads_region(
        bam_path,
        chrom,
        win_start,
        win_end,
        read_uniq_list,
        variants_dict,
        read_mapping,
        compute_missing_unknown,
        padding= padding
    )

    # stats
    m_reads= len(read_uniq_list)
    t_reads= len(df_reads)
    p_reads= (100.0*m_reads/t_reads) if t_reads>0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable reads: {m_reads} ({p_reads:.2f}%)")

    m_vars= len(final_variants)
    t_vars= len(union_variants)
    p_vars= (100.0*m_vars/t_vars) if t_vars>0 else 0
    logging.info(f"[CHUNK] {chunk_str} # mappable variants: {m_vars} ({p_vars:.2f}%)")

    # save chunk
    chunk_data= {
      "union_variants": final_variants,
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
# 5. Save final matrix
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
            dmap[(r,co)] +=1

        for (r,co), val in dmap.items():
            row_inds.append(r)
            col_inds.append(co)
            data.append(val)

        mat= csr_matrix((data, (row_inds,col_inds)), shape=(len(variants), len(barcodes)))
        out_h5= os.path.join(save_dir, f"{name}_matrix.h5")
        with h5py.File(out_h5,"w") as hf:
            hf.create_dataset("data", data=mat.data)
            hf.create_dataset("indices", data=mat.indices)
            hf.create_dataset("indptr", data=mat.indptr)
            hf.attrs['shape']= mat.shape

    build_csr_and_save(ref_classifications, "ref")
    build_csr_and_save(alt_classifications, "alt")
    if missing_classifications:
        build_csr_and_save(missing_classifications, "missing")
    if unknown_classifications:
        build_csr_and_save(unknown_classifications, "unknown")

    joblib.dump((variants, barcodes), os.path.join(save_dir,"variant_barcode_mappings.pkl"))

############################################################################
# 6. main => chunk-level parallel
############################################################################

def main():
    args, vcf_files_with_origins= parse_arguments()
    setup_logging(args.save_dir)
    logging.info("Command Line: " + " ".join(sys.argv))
    logging.info("=== [scVarID] Program Started ===")

    compute_missing_unknown= not args.ref_alt_only

    # step: read chromosome lengths
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

    wsize= args.window_size
    pad= args.padding
    num_cores= args.num_cores
    logging.info(f"[INFO] chunk-level parallel => n_jobs={num_cores}, window_size={wsize}, padding={pad}")

    # create chunk list
    step_t= log_step_start("Generate chunk list")
    all_chunks=[]
    for c_, length_ in chrom_len.items():
        st=1
        while st<= length_:
            e_= min(st+ wsize -1, length_)
            all_chunks.append((c_, st, e_))
            st= e_+1
    logging.info(f"Total chunk count: {len(all_chunks)}")
    log_step_end("Generate chunk list", step_t)

    # parallel chunk
    tmp_dir= os.path.join(args.save_dir,"tmp_chunks")
    os.makedirs(tmp_dir, exist_ok=True)

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

    # save matrix
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
