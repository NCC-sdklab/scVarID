# classification.py

import logging
import re
from collections import defaultdict
from joblib import Parallel, delayed
import pysam
import tempfile
import numpy as np

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

def process_read_chunk(
    bam_file_path, 
    read_unique_names_set, 
    variants_dict, 
    read_mapping, 
    compute_missing_unknown=True
):
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
                ref, alt, missing, unknown = process_read(
                    read, variants_dict, compute_missing_unknown
                )
                ref_classifications.extend(ref)
                alt_classifications.extend(alt)
                missing_classifications.extend(missing)
                unknown_classifications.extend(unknown)

                read_unique_names_set.remove(read_unique_name)
                if not read_unique_names_set:
                    break

    return ref_classifications, alt_classifications, missing_classifications, unknown_classifications

# def process_bam_data_parallel(
#     bam_path, 
#     selected_read_unique_names, 
#     selected_variants, 
#     read_mapping, 
#     num_cores=4, 
#     compute_missing_unknown=True
# ):
#     """
#     Process BAM data in parallel and return classification results.
#     """
#     variants_dict = create_variants_dict(selected_variants)
#     read_unique_names_chunks = split_list(selected_read_unique_names, num_cores)

#     # Initialize result lists
#     ref_classifications = []
#     alt_classifications = []
#     missing_classifications = []
#     unknown_classifications = []

#     # Parallel processing
#     results = Parallel(n_jobs=num_cores, backend='multiprocessing')(
#         delayed(process_read_chunk)(
#             bam_path, 
#             set(chunk), 
#             variants_dict, 
#             read_mapping, 
#             compute_missing_unknown
#         ) for chunk in read_unique_names_chunks
#     )

#     for result in results:
#         ref_classifications.extend(result[0])
#         alt_classifications.extend(result[1])
#         missing_classifications.extend(result[2])
#         unknown_classifications.extend(result[3])

#     return ref_classifications, alt_classifications, missing_classifications, unknown_classifications
def process_bam_data_parallel(
    bam_path, 
    selected_read_unique_names, 
    selected_variants, 
    read_mapping, 
    num_cores=4, 
    compute_missing_unknown=True
):
    variants_dict = create_variants_dict(selected_variants)

    # 메모리 공유를 위한 청크 분할
    chunk_factor = 4  # 코어당 4개 청크
    #chunk_factor = 8  # 코어당 8개 청크
    chunks = np.array_split(
        selected_read_unique_names, 
        num_cores * chunk_factor
    )

    # 병렬 처리 설정
    results = Parallel(
        n_jobs=num_cores,
        backend='loky',  # 메모리 공유 백엔드
        #temp_folder='/dev/shm',  # 공유 메모리 사용
        temp_folder='/tmp',  # 공유 메모리 사용
        max_nbytes='256M',  # 직렬화 크기 제한
        prefer="processes"  # 프로세스 기반 병렬화
    )(
        delayed(process_read_chunk)(
            bam_path, 
            set(chunk), 
            variants_dict, 
            read_mapping, 
            compute_missing_unknown
        ) for chunk in chunks
    )

    # 결과 통합
    ref_classifications = []
    alt_classifications = []
    missing_classifications = []
    unknown_classifications = []
    
    for res in results:
        ref_classifications.extend(res[0])
        alt_classifications.extend(res[1])
        missing_classifications.extend(res[2])
        unknown_classifications.extend(res[3])
        
    return ref_classifications, alt_classifications, missing_classifications, unknown_classifications

def process_window(bam_path, chrom, start, end, variants_dict, read_mapping, compute_missing_unknown):
    """개별 window 처리 함수"""
    window_variants = {k: v for k, v in variants_dict.items() if k[0] == chrom and start <= k[1] <= end}
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        reads = list(bam.fetch(chrom, start, end))
    
    # 임시 파일 사용
    with tempfile.NamedTemporaryFile(mode='w+t', delete=True) as tmp:
        for read in reads:
            tmp.write(f"{read.to_string()}\n")
        tmp.seek(0)
        
        return process_read_chunk(
            tmp.name, 
            set(selected_reads), 
            window_variants, 
            read_mapping, 
            compute_missing_unknown
        )

def process_read(read, variants_dict, compute_missing_unknown=True):
    ref_classifications_local = set()
    alt_classifications_local = set()
    missing_classifications_local = set()
    unknown_classifications_local = set()

    read_name, cell_barcode, sequences, positions, operations = prime_read(read)
    chrom = read.reference_name

    read_start = read.reference_start + 1
    read_end = read.reference_end

    processed_positions = set()

    for pos in positions:
        key = (chrom, pos)

        if key in variants_dict and key not in processed_positions:
            processed_positions.add(key)

            for variant_data in variants_dict[key]:
                variant = variant_data['variant']
                variant_type = variant_data['type']

                ref_len = len(variant[2])
                alt_len = len(variant[3])

                if pos < read_start or pos > read_end - (alt_len - 1) or pos > read_end - (ref_len - 1):
                    continue

                classification, barcode, variant_id, read_id = classify_variant(
                    sequences, positions, operations, pos, variant[2], variant[3], variant_type,
                    operations[positions.index(pos)], read_name, cell_barcode
                )

                classification_data = {
                    'read_name': read_id,
                    'cell_barcode': barcode,
                    'variant': f"{variant[0]}:{variant[1]}_{variant[2]}->{variant[3]}",
                }

                if classification == 'ref':
                    ref_classifications_local.add(frozenset(classification_data.items()))
                elif classification == 'alt':
                    alt_classifications_local.add(frozenset(classification_data.items()))
                elif classification == 'missing' and compute_missing_unknown:
                    missing_classifications_local.add(frozenset(classification_data.items()))
                elif classification == 'unknown' and compute_missing_unknown:
                    unknown_classifications_local.add(frozenset(classification_data.items()))

    return (
        ref_classifications_local,
        alt_classifications_local,
        missing_classifications_local,
        unknown_classifications_local,
    )

def create_variants_dict(selected_variants):
    variants_dict = {}
    for variant in selected_variants:
        chrom = variant['Chromosome']
        pos = int(variant['POS'])
        ref = variant['REF']
        alt_list = variant['ALT'].split(',')  # Handle cases where ALT values are multiple

        for alt in alt_list:
            if len(ref) > len(alt):
                variant_type = 'deletion'
            elif len(ref) < len(alt):
                variant_type = 'insertion'
            else:
                variant_type = 'snv'

            key = (chrom, pos)
            if key not in variants_dict:
                variants_dict[key] = []
            variants_dict[key].append({'variant': (chrom, pos, ref, alt), 'type': variant_type})
    return variants_dict

def classify_variant(sequences, positions, operations, pos, ref, alt, variant_type, operation, read_name, cell_barcode):
    if variant_type == 'deletion':
        return handle_deletion(sequences, positions, operations, pos, ref, alt, read_name, cell_barcode)
    elif variant_type == 'insertion':
        return handle_insertion(sequences, positions, operations, pos, ref, alt, read_name, cell_barcode)
    else:  # SNV or other types
        return handle_snv(sequences, positions, pos, ref, alt, read_name, cell_barcode)

def prime_read(read):
    read_name = read.query_name
    try:
        cell_barcode = read.get_tag("CB")
    except KeyError:
        cell_barcode = None

    sequences, positions, operations = [], [], []
    current_position = read.reference_start + 1
    current_read_pos = 0

    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', read.cigarstring)

    for length, op in cigar_ops:
        length = int(length)
        if op in "M=X":
            for _ in range(length):
                sequences.append(read.query_sequence[current_read_pos])
                positions.append(current_position)
                operations.append(op)
                current_position += 1
                current_read_pos += 1
        elif op == "I":
            for _ in range(length):
                sequences.append(read.query_sequence[current_read_pos])
                positions.append(-1)  # Use -1 as the standard representation for insertion positions
                operations.append(op)
                current_read_pos += 1
        elif op in "DNP":
            for _ in range(length):
                sequences.append('-')
                positions.append(current_position)
                operations.append(op)
                current_position += 1
        elif op == "S":
            current_read_pos += length
        elif op == "H":
            pass

    return read_name, cell_barcode, sequences, positions, operations

def handle_snv(sequences, positions, pos, ref, alt, read_name, cell_barcode):
    for p, seq in zip(positions, sequences):
        if p == pos:
            if seq == ref:
                return 'ref', cell_barcode, pos, read_name
            elif seq == alt:
                return 'alt', cell_barcode, pos, read_name
            elif seq == '-':
                return 'missing', cell_barcode, pos, read_name
            else:
                return 'unknown', cell_barcode, pos, read_name
    return 'unknown', cell_barcode, pos, read_name

def handle_insertion(sequences, positions, operations, pos, ref, alt, read_name, cell_barcode):
    try:
        pos_index = positions.index(pos)

        # Check for missing status first
        if sequences[pos_index] == '-':
            return 'missing', cell_barcode, pos, read_name

        # Check for ref status
        elif (
            all(op == '=' for op in operations[pos_index : pos_index + len(ref)]) and
            operations[pos_index + len(ref)] == '=' 
        ):
            return 'ref', cell_barcode, pos, read_name

        # Check for alt status
        elif (
            all(op == '=' for op in operations[pos_index : pos_index + len(ref)]) and
            sequences[pos_index + len(ref) : pos_index + len(alt) - len(ref) + 1] == list(alt[len(ref):]) and
            all(op == 'I' for op in operations[pos_index + len(ref) : pos_index + len(alt) - len(ref) + 1]) and
            operations[pos_index + len(alt) - len(ref) + 1] == '='
        ):
            return 'alt', cell_barcode, pos, read_name

    except IndexError:
        return 'unknown', cell_barcode, pos, read_name

    return 'unknown', cell_barcode, pos, read_name

def handle_deletion(sequences, positions, operations, pos, ref, alt, read_name, cell_barcode):
    try:
        pos_index = positions.index(pos)

        # Check for missing status first
        if sequences[pos_index] == '-':
            return 'missing', cell_barcode, pos, read_name

        # Check for ref status
        elif(
            all(op == '=' for op in operations[pos_index : pos_index + len(ref)]) and
            operations[pos_index + len(ref)] == '='
        ):
            return 'ref', cell_barcode, pos, read_name

        # Check for alt status
        elif (
            all(op == '=' for op in operations[pos_index : pos_index + len(alt)]) and
            all(seq == '-' for seq in sequences[pos_index + len(alt) : pos_index + len(ref) - len(alt) + 1]) and
            all(op == 'D' for op in operations[pos_index + len(alt) : pos_index + len(ref) - len(alt) + 1]) and
            operations[pos_index + len(ref) - len(alt) + 1] == '='
        ):
            return 'alt', cell_barcode, pos, read_name

    except IndexError:
        return 'unknown', cell_barcode, pos, read_name

    return 'unknown', cell_barcode, pos, read_name

# classification.py에 추가할 보조 함수

def create_variants_dict(selected_variants):
    variants_dict = defaultdict(list)
    for variant in selected_variants:
        chrom = variant['Chromosome']
        pos = int(variant['POS'])
        ref = variant['REF']
        alt_list = variant['ALT'].split(',')
        
        for alt in alt_list:
            key = (chrom, pos)
            variant_type = (
                'deletion' if len(ref) > len(alt) else
                'insertion' if len(ref) < len(alt) else 
                'snv'
            )
            variants_dict[key].append({
                'variant': (chrom, pos, ref, alt),
                'type': variant_type
            })
    return variants_dict
