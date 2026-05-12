import os, sys, time
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
from functools import wraps
import dask.bag as db


# Costanti di configurazione
KMER_LEN_MIN = 8
KMER_LEN_MAX = 12
VALID_BASES = "ATCGN"
COMPLEMENTARY = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
IUPAC_CODES = {
    'A':['A'], 'T':['T'], 'C':['C'], 'G':['G'], 'N':['A', 'T', 'C', 'G'],
    'R':['A', 'G'], 'Y':['C', 'T'], 'S':['G', 'C'], 'W':['A', 'T'],
    'K':['G', 'T'], 'M':['A', 'C'], 'B':['C', 'G', 'T'],
    'D':['A', 'G', 'T'], 'H':['A', 'C', 'T'], 'V':['A', 'C', 'G']
}
Q_SCORE_DICT_33 = {chr(i):i for i in range(33, 85)}  # Phred-33
Q_SCORE_DICT_64 = {chr(i):i for i in range(64, 116)}  # Phred-64
CHECK_OFF_SCORE_33_63 = set(chr(i) for i in range(33, 63))
CHECK_OFF_SCORE_83_103 = set(chr(i) for i in range(85, 103))


#######################################################################################
#################       AGGIORNALO PER IL TUO PERCORSO        #########################
FOLD_PATH = "C:/Users/andre/Desktop/my_project/combo_pf/prjct_gamma/" 
#FOLD_PATH = "C:/Users/admin/Desktop/my_project/combo_pf/prjct_gamma/"
#######################################################################################


# Decoratore per misurare il tempo
def speedtest(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed_time = time.perf_counter() - start_time
        hours = int(elapsed_time / 3600)
        minutes = int((elapsed_time - hours * 3600) / 60)
        seconds = round(elapsed_time - hours * 3600 - minutes * 60, 3)
        print(f"Exec time of '{func.__name__}':{hours:02}:{minutes:02}:{seconds:03.3f} [hh:mm:ss]")
        return result
    return wrapper

# Validazione delle sequenze e quality scores
def is_valid_seq(sequence:str) -> bool:
    return all(char in VALID_BASES for char in sequence)

def is_valid_qs(qscore:str, phred:int = 33) -> bool:
    score_dict = Q_SCORE_DICT_33 if phred == 33 else Q_SCORE_DICT_64
    return all(char in score_dict for char in qscore)

# Conversione dei quality scores
def conv_qs(qscore_list:list[str]) -> list[list[int]]:
    def detect_phred(qscores:str) -> int:
        combined = ''.join(qscores)
        if any(c in CHECK_OFF_SCORE_33_63 for c in combined) and not any(c in CHECK_OFF_SCORE_83_103 for c in combined):
            print("Detected sequencing method: Sanger (Phred-33)")
            return 33
        print("Detected sequencing method:Illumina (Phred-64)")
        return 64

    phred = detect_phred(qscore_list)
    score_dict = Q_SCORE_DICT_33 if phred == 33 else Q_SCORE_DICT_64
    offset = phred
    return [[score_dict[char] - offset for char in qs] for qs in qscore_list]

# Reverse e complemento delle sequenze
def rev_seqlist(seq_list:list[str]) -> list[str]:
    return [''.join(reversed(seq)) for seq in seq_list]

def compl_seqlist(seq_list:list[str]) -> list[str]:
    return [''.join(COMPLEMENTARY.get(base, 'N') for base in seq) for seq in seq_list]

# Parsing FASTA
def fa_pars_1cpu(file_path:str) -> tuple[list[str], list[str]]:
    idseq_list, seq_list = [], []
    try:
        fasta = list(db.read_text(file_path))
        if not fasta or not fasta[0].startswith('>'):
            raise ValueError(f"Corrupted FASTA file:No ID at row 0 in '{file_path}'")
        
        i = 0
        while i < len(fasta):
            line = fasta[i].strip().upper()
            if line.startswith(">"):
                idseq_list.append(line.lstrip('>').strip())
                seq = []
                i += 1
                while i < len(fasta) and not fasta[i].startswith('>'):
                    curr_line = fasta[i].strip().upper()
                    if not is_valid_seq(curr_line):
                        raise ValueError(f"Corrupted FASTA file:Invalid sequence at row {i} in '{file_path}'")
                    seq.append(curr_line)
                    i += 1
                seq_list.append(''.join(seq))
            else:
                i += 1
        if len(idseq_list) != len(seq_list):
            raise ValueError(f"Corrupted FASTA file:Mismatch between IDs and sequences in '{file_path}'")
    except Exception as e:
        print(f"Error parsing FASTA:{e}")
        sys.exit(1)
    return idseq_list, seq_list

# Funzione di parsing per i chunk FASTQ (spostata al livello superiore)
def fq_pars_chunk(file_path:str, chunk:list[str], start_idx:int) -> tuple[list[str], list[str], list[str], list[str]]:
    acc_num_list, seq_list, comment_list, qscore_list = [], [], [], []
    for i in range(0, len(chunk), 4):
        try:
            if not chunk[i].startswith('@'):
                print(f"Corrupted FASTQ:No ID at row {start_idx + i} in '{file_path}'")
                return [], [], [], []
            acc_num_list.append(chunk[i].lstrip('@').strip())
            seq = chunk[i+1].strip()
            if not is_valid_seq(seq):
                print(f"Corrupted FASTQ:Invalid sequence at row {start_idx + i + 1} in '{file_path}'")
                return [], [], [], []
            seq_list.append(seq)
            if not chunk[i+2].startswith('+'):
                print(f"Corrupted FASTQ:No comment at row {start_idx + i + 2} in '{file_path}'")
                return [], [], [], []
            comment_list.append(chunk[i+2].strip())
            qs = chunk[i + 3].strip()
            if not (is_valid_qs(qs, 33) or is_valid_qs(qs, 64)) or len(qs) != len(seq):
                print(f"Corrupted FASTQ: Invalid quality score at row {start_idx + i + 3} in '{file_path}'")
                return [], [], [], []
            qscore_list.append(qs)
        except IndexError:
            print(f"Corrupted FASTQ:Incomplete record at row {start_idx + i} in '{file_path}'")
            return [], [], [], []
    return acc_num_list, seq_list, comment_list, qscore_list

# Parsing FASTQ in parallelo
@speedtest
def fq_pars_hpc(file_path:str) -> tuple[list[str], list[str], list[str], list[str]]:
    fastq_file = list(db.read_text(file_path))
    if not fastq_file or len(fastq_file) % 4 != 0:
        raise ValueError(f"Corrupted FASTQ:File '{file_path}' has incomplete records")
    
    n_cpu = max(1, cpu_count()-1)
    chunk_size = (len(fastq_file)//(n_cpu*4))*4 or 4  # Assicura chunk validi
    chunks = [fastq_file[i:i+chunk_size] for i in range(0, len(fastq_file), chunk_size)]
    
    with Pool(processes=n_cpu) as pool:
        results = pool.starmap(fq_pars_chunk, [(file_path, chunk, i*chunk_size) for i, chunk in enumerate(chunks)])
    
    acc_num_list, seq_list, comment_list, qscore_list = [], [], [], []
    for res in results:
        if not res[0]: # Se un chunk fallisce, skippa il risultato
            continue
        acc_num_list.extend(res[0])
        seq_list.extend(res[1])
        comment_list.extend(res[2])
        qscore_list.extend(res[3])
    
    print(f"Parsed '{file_path[len(FOLD_PATH):]}': {len(seq_list)} sequences")
    return acc_num_list, seq_list, comment_list, qscore_list

# Controllo dell'expected error
def check_ee(qscore_num_list:list[list[int]], max_err:float=2.0) -> list[int]:
    ee_list = [1 if sum([10**(q/-10) for q in qs]) <= max_err else 0 for qs in qscore_num_list]
    good_seq = sum(ee_list)
    total_seq = len(ee_list)
    curr_max_err = max_err
    
    while good_seq < total_seq*0.65:
        print(f"Sequences with EE <= {curr_max_err}: {good_seq}/{total_seq}")
        try:
            choice = input("Keep current sequences (0) or increase EE threshold (1)? ").strip()
            if choice == "0":
                break
            elif choice == "1":
                curr_max_err += 0.5
                ee_list = [1 if sum(10**(q/-10) for q in qs)<=curr_max_err else 0 for qs in qscore_num_list]
                good_seq = sum(ee_list)
            else:
                print("Invalid input, please enter 0 or 1")
        except KeyboardInterrupt:
            print("Process interrupted by user")
            break
    return ee_list

# Esclusione delle sequenze duplicate o con errore alto
def exclude_seq(acc_num_list:list[str], seq_list:list[str], comment_list:list[str], 
                qscore_num_list:list[list[int]], err_list:list[int]) -> list[list]:
    mean_qs = [sum(qs) / len(qs) for qs in qscore_num_list]
    valid_seq = [[acc_num_list[i], seq_list[i], comment_list[i], qscore_num_list[i]] 
                 for i in range(len(seq_list)) if err_list[i] == 1]
    seq_to_idx = {}
    for i, seq in enumerate(valid_seq):
        seq_str = seq[1]
        if seq_str in seq_to_idx:
            prev_idx = seq_to_idx[seq_str]
            if mean_qs[i] > mean_qs[prev_idx]:
                seq_to_idx[seq_str] = i
        else:
            seq_to_idx[seq_str] = i
    
    return [valid_seq[idx] for idx in sorted(seq_to_idx.values())]

# Merging delle paired-end reads
def merge_pe(f1_seq_list:list[str], f2_rc_seq_list:list[str], min_overlap:int=10) -> list[list]:
    merged_sequences_info = []
    for i, seq1 in enumerate(f1_seq_list):
        for j, seq2 in enumerate(f2_rc_seq_list):
            len_seq1 = len(seq1)
            for k in range(len_seq1-min_overlap):
                overlap = seq1[k:]
                if seq2.startswith(overlap):
                    merged_sequences_info.append([seq1[:k] + seq2, k, i, j])
                    break
    return merged_sequences_info

# Indice di Jaccard pesato
def weight_jcrd_idx(seq1:str, seq2:str, len_kmer:int) -> float:
    def weight_kmers(seq:str, len_kmer:int, alpha:float=0.25, max_n:int=len_kmer//2) -> dict[str, float]:
        kmers = {}
        for i in range(len(seq) - len_kmer + 1):
            kmer = seq[i:i + len_kmer]
            n_count = kmer.count('N')
            if n_count < max_n:
                kmers[kmer] = max(kmers.get(kmer, 0), alpha ** n_count)
        return kmers
    
    k1, k2 = weight_kmers(seq1, len_kmer), weight_kmers(seq2, len_kmer)
    union_keys = set(k1.keys()) | set(k2.keys())
    num = sum(min(k1.get(k, 0), k2.get(k, 0)) for k in union_keys)
    den = sum(max(k1.get(k, 0), k2.get(k, 0)) for k in union_keys)
    return num/den if den > 0 else 0.0

# Indice di Jaccard non pesato
def jcrd_idx(seq1:str, seq2:str, len_kmer:int) -> float:
    k1 = set(seq1[i:i+len_kmer] for i in range(len(seq1)-len_kmer+1))
    k2 = set(seq2[i:i+len_kmer] for i in range(len(seq2)-len_kmer+1))
    intersection = k1 & k2
    union = k1 | k2
    return len(intersection) / len(union) 

# Inizializzazione dei file di output
def file_init(fold_path:str, n_pool_cpu:int) -> None:
    folder_path = f"{fold_path}cluster_file"
    os.makedirs(folder_path, exist_ok=True)
    for file_idx in range(n_pool_cpu):
        open(f"{folder_path}/OTU_clustered_seq_{file_idx}.txt", "wb").close()

# Clustering OTU
def otu_cluster(args:tuple[str, tuple[list[str], list[str]], list[list], list[list[str]], int, int, int]) -> None:
    fold_path, ref_seq_fa, merged_sequences_info, f1f2_acc_num_list, file_idx, idx_range = args
    start_idx = file_idx * idx_range
    end_idx = min(start_idx+idx_range, len(ref_seq_fa[0]))
    
    with open(f"{fold_path}cluster_file/OTU_clustered_seq_{file_idx}.txt", "w", encoding='utf-8') as f:
        for rf in range(start_idx, end_idx):
            f.write(f"{ref_seq_fa[0][rf]}\n")
            for ms in merged_sequences_info:
                for kmer_len in range(KMER_LEN_MIN, KMER_LEN_MAX + 1):
                    #method = "W" if "N" in ref_seq_fa[1][rf] else "N"
                    method = "N"
                    j_idx = weight_jcrd_idx(ms[0], ref_seq_fa[1][rf], kmer_len) if method == "W" else \
                            jcrd_idx(ms[0], ref_seq_fa[1][rf], kmer_len)
                    if j_idx > 0.0:
                        id_seq1 = f1f2_acc_num_list[0][ms[2]]
                        id_seq2 = f1f2_acc_num_list[1][ms[3]]
                        f.write(f"[k={kmer_len}] [{method}ji={j_idx:.3f}] [{id_seq1}||{id_seq2}]\t\t{ms[0]}\n")
                        break

# Main
if __name__ == "__main__":
    mp.freeze_support()  # Supporto per Windows e Linux
    
    # File di input
    F1 = f"{FOLD_PATH}merging_test1.fq.gz"
    F2 = f"{FOLD_PATH}merging_test2.fq.gz"
    REF_SEQ_FILE = f"{FOLD_PATH}ref.fa"

    # Parsing FASTQ
    print("Parsing F1...")
    f1_acc_num_list, f1_seq_list, f1_comment_list, f1_qscore_list = fq_pars_hpc(F1)
    f1_qscore_num_list = conv_qs(f1_qscore_list)
    f1_err_list = check_ee(f1_qscore_num_list)

    print("\nParsing F2...")
    f2_acc_num_list, f2_seq_list, f2_comment_list, f2_qscore_list = fq_pars_hpc(F2)
    f2_qscore_num_list = conv_qs(f2_qscore_list)
    f2_err_list = check_ee(f2_qscore_num_list)

    # Filtraggio delle sequenze
    print("\nFiltering sequences...")
    f1_valid_seq = exclude_seq(f1_acc_num_list, f1_seq_list, f1_comment_list, f1_qscore_num_list, f1_err_list)
    f2_valid_seq = exclude_seq(f2_acc_num_list, f2_seq_list, f2_comment_list, f2_qscore_num_list, f2_err_list)
    print(f"Valid F1 sequences: {len(f1_valid_seq)}\nValid F2 sequences: {len(f2_valid_seq)}")

    # Merging delle paired-end reads
    print("\nMerging paired-end reads...")
    f1_seq_list = [seq[1] for seq in f1_valid_seq]
    f2_rc_seq_list = rev_seqlist(compl_seqlist([seq[1] for seq in f2_valid_seq]))
    merged_sequences_info = merge_pe(f1_seq_list, f2_rc_seq_list)
    print(f"Number of merged reads: {len(merged_sequences_info)}")

    print (merged_sequences_info)


    # Parsing FASTA di riferimento
    print("\nParsing reference FASTA...")
    ref_seq_fa = fa_pars_1cpu(REF_SEQ_FILE)
    print(f"Number of reference sequences: {len(ref_seq_fa[0])}")

    # Clustering OTU
    print("\nStarting OTU clustering...")
    n_cpu = max(1, cpu_count()-1)
    idx_range = len(ref_seq_fa[0])//n_cpu + (len(ref_seq_fa[0])%n_cpu > 0)
    f1f2_acc_num_list = [[seq[0] for seq in f1_valid_seq], [seq[0] for seq in f2_valid_seq]]
    
    file_init(FOLD_PATH, n_cpu)
    with Pool(processes=n_cpu) as pool:
        pool.map(otu_cluster, [(FOLD_PATH, ref_seq_fa, merged_sequences_info, f1f2_acc_num_list, i, idx_range) 
                              for i in range(n_cpu)])
    print("OTU clustering completed.")