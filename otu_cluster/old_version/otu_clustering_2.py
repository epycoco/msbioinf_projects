from multiprocessing import Pool, cpu_count, Process, Queue
import numpy as np

from config import *
from itertools import product
import dask.bag as db
import sys, os
import json, subprocess, tempfile

from itertools import combinations
from functools import reduce
import operator

#FOLD_PATH = "C:/Users/admin/Desktop/my_project/combo_pf/prjct_gamma/"
FOLD_PATH = "C:/Users/andre/Desktop/my_project/combo_pf/prjct_gamma/"
sys.path.append(FOLD_PATH)

# consente la verifica della sequenza controllando che non ci siano caratteri diversi da quelli attesi
def is_valid_seq(sequence:str) -> bool:
    return all(char in VALID_BASES for char in sequence)

# consente la verifica del quality score controllando che non ci siano caratteri diversi da quelli attesi 
def is_valid_qs33(qscore:str) -> bool:
    return all(char in Q_SCORE_DICT_33 for char in qscore)

# consente la verifica del quality score controllando che non ci siano caratteri diversi da quelli attesi
def is_valid_qs64(qscore:str) -> bool:
    return all(char in Q_SCORE_DICT_64 for char in qscore)

# converte la lista di stringhe quality score in una lista di liste di interi, analizzando quale metodo di sequenziamento è stato usato, per consentire la corretta conversione
def conv_qs (qscore_list):
    def off_score(qscore_list:str):
        str_qscore = ''.join(qscore_list)
        if any(char in CHECK_OFF_SCORE_33_63 for char in str_qscore) and not any(char in CHECK_OFF_SCORE_83_103 for char in str_qscore):
            print("Sequencing Method: Sanger")
            return True #off score 33
        print("Sequencing Method: Illumina")
        return False #off score 64
    def conv_qs33 (qscore_list:list[str]) -> list[list[int]]:
        len_qscore_list = len(qscore_list)
        qscore_num_list = []  
        for i in range(len_qscore_list):
            qscore_num_list_i = []
            for j in range(len(qscore_list[i])):
                qscore_num_list_i.append(Q_SCORE_DICT_33[qscore_list[i][j]] - 33)  
            qscore_num_list.append(qscore_num_list_i) 
        return qscore_num_list
    def conv_qs64 (qscore_list:list[str]) -> list[list[int]]:
        len_qscore_list = len(qscore_list)
        qscore_num_list = []  
        for i in range(len_qscore_list):
            qscore_num_list_i = [] 
            for j in range(len(qscore_list[i])):
                qscore_num_list_i.append(Q_SCORE_DICT_64[qscore_list[i][j]] - 64)  
            qscore_num_list.append(qscore_num_list_i) 
        return qscore_num_list
    
    if off_score(qscore_list):
        return conv_qs33(qscore_list)
    else:
        return conv_qs64(qscore_list)

# prende una lista di sequenze per ottenere una lista di sequenze revertite
def rev_seqlist (seq_list:list[str]):
    len_seqlist = len(seq_list)
    rev_seqlist = [""]*len_seqlist
    for i in range (len_seqlist):
        rev_seqlist[i] = ''.join(list(reversed(seq_list[i])))
    return rev_seqlist

# prende una lista di sequenze per ottenere una lista di sequenze complementari
def compl_seqlist(seq_list:list[str]):
    len_seqlist = len(seq_list)
    complement_seqlist = [""]*len_seqlist
    for i in range (len_seqlist):
        complement_seqlist[i] = ''.join([COMPLEMENTARY[base] for base in seq_list[i]])
    return complement_seqlist



def fa_pars_1cpu (file_path) -> tuple[str, str]:
    idseq_list, seq_list = [], []
    fasta = list(db.read_text(file_path))
    len_fasta = len(fasta)

    if not (str(fasta[0])).startswith('>'):
        print(f"Corrupted File. No correct ID at row 0 in file '{file_path}'")
        sys.exit("Programm Closed.")
    
    for i in range (len_fasta):
        l = (str(fasta[i]).strip('\n')).upper()
        if l.startswith(">"):
            idseq_list.append(l.strip('>').strip('\n'))
            seq = []
            for j in range (i+1, len_fasta):
                l = str(fasta[j]).strip('\n').upper()
                if is_valid_seq(l):
                    seq.append(l)
                elif l.startswith(">"):
                    i == j
                    seq_list.append (''.join(map(str, seq)))
                    break
                elif not is_valid_seq(l):
                    print(f"Corrupted File. No correct sequence at row {j} in file '{file_path}'")
                    sys.exit("Programm Closed.")

    return idseq_list, seq_list

# sottofunzione di 'fq_pars_hpc' utile per la suddivisione in chunk del parsing parallelizzato di un fastq
def fq_pars_chunk(file_path:str, fastq_file:list[str], idx_cpu:int, range_per_pool:int):  
    acc_num_list, seq_list, comment_list, qscore_list = [], [], [], []
    c = 0

    for i in range(range_per_pool):
        idx = i + idx_cpu * range_per_pool
        line = fastq_file[idx].strip()
        if c == 0:  # Accession number
            if line.startswith("@"):
                acc_num_list.append(line.strip("@").strip('\n'))
            else:
                print(f"Corrupted File. No correct ID at row {idx} in file '{file_path}'\n")
                return [], [], [], []
        elif c == 1:  # Sequence
            if is_valid_seq(line):
                seq_list.append(line)
            else:
                print(f"Corrupted File. No correct sequence at row {idx} in file '{file_path}'\n")
                return [], [], [], []
        elif c == 2:  # Comment
            if line.startswith("+"):
                comment_list.append(line)
            else:
                print(f"Corrupted File. No correct comment at row {idx} in file '{file_path}'\n")
                return [], [], [], []
        elif c == 3:  # Quality score
            seql = fastq_file[idx-2].strip()
            if is_valid_qs33(line) or is_valid_qs64(line):
                if len(line) != len(seql):
                    print(f"The length of sequence (line {idx-2}) does not correspond to the length of quality score (line {idx})\n")
                    return [], [], [], []
                else:
                    qscore_list.append(line)
            else:
                print(f"Corrupted File. Invalid char in quality score at row {idx} in file '{file_path}'")
                return [], [], [], []
        c += 1
        if c == 4:
            c = 0
    return acc_num_list, seq_list, comment_list, qscore_list

# parsing parallelizzato di un fastq 
@speedtest
def fq_pars_hpc(file_path:str, CPU_PARS=cpu_count()-1) -> tuple[list[str], list[str], list[str], list[str]]:
    fastq_file = list(db.read_text(file_path))
    len_fastq = len(fastq_file)
    print(f"{file_path}\nNum file rows: {len_fastq}")
    
    ## Definisco il numero di sequenze che ogni core deve analizzare
    if len_fastq % CPU_PARS == 0:
        range_per_pool = len_fastq // CPU_PARS
    else:
        # Nel caso in cui il numero di CPU non dovesse essere compatibile con il numero di sequenze presenti nel file,
        #   RICALCOLO il valore di CPU_PARS (numero di CPU da usare per il parsing parallelizzato)
        def cpus_to_use(len_fastq_file):
            def prime_fktr(n):
                fktrs = []
                while n % 2 == 0:
                    fktrs.append(2)
                    n //= 2
                i = 3
                while i * i <= n:
                    while n % i == 0:
                        fktrs.append(i)
                        n //= i
                    i += 2
                if n > 2:
                    fktrs.append(n)
                return fktrs
            def product_calc(fktrs):
                products = set()
                for r in range(2, len(fktrs) + 1):  # Consideriamo combinazioni di almeno 2 fattori
                    for combo in combinations(fktrs, r):
                        prdct = reduce(operator.mul, combo)
                        products.add(prdct)
                return sorted(products)  
            
            fktr_len_fastq = prime_fktr(len_fastq_file)
            combos = product_calc(fktr_len_fastq)
            # Seleziono il valore più alto possibile senza superare CPU_PARS
            valid_products = [p for p in combos if p <= CPU_PARS]
            return max(valid_products) if valid_products else 4  # Ritorno il massimo valore valido o 4 se nessuno è valido

        CPU_PARS = cpus_to_use(len_fastq)
        range_per_pool = len_fastq // CPU_PARS 

    print (f"Num of CPU: {CPU_PARS}")

    len_lists = len_fastq // 4
    acc_num_list = [""] * len_lists
    seq_list = [""] * len_lists
    comment_list = [""] * len_lists
    qscore_list = [""] * len_lists

    pool = Pool(processes=CPU_PARS)
    results = [pool.apply_async(fq_pars_chunk, args=(file_path, fastq_file, idx_cpu, range_per_pool)) for idx_cpu in range(CPU_PARS)]
    pool.close()
    pool.join()

    # Raccolgo i risultati dai processi figli
    acc_num_list, seq_list, comment_list, qscore_list = [], [], [], []
    for res in results:
        acc_num, seq, comment, qscore = res.get()
        acc_num_list.extend(acc_num)
        seq_list.extend(seq)
        comment_list.extend(comment)
        qscore_list.extend(qscore)

    print(f"Parsing COMPLETED")
    print(f"Number of sequences: {len(seq_list)}")
    return acc_num_list, seq_list, comment_list, qscore_list

# Funzione per controllare l'expected error.
# Permette all'utente di calidare anche sequenze che hanno un expected error maggiore di 2.
#   Il programma però non tiene in considerazione 
def check_ee(qscore_num_list:list[list[int]]):
    def exp_errq(qscore_num_list: list[list[int]], max_err=2) -> list[list[float]]:
        n_seq = len(qscore_num_list)
        expected_error_list = [[] for _ in range(n_seq)]
        for i in range(n_seq):
            expected_error_list[i] = sum(map(lambda q: 10**(q/-10), qscore_num_list[i]))
            if expected_error_list[i]<=max_err:
                expected_error_list[i] = 1
            else:
                expected_error_list[i] = 0
        return expected_error_list

    err_list = exp_errq(qscore_num_list)
    n_seq = len(err_list)
    good_seq = err_list.count(1)
    ee = 2.0
    while True:
        if good_seq < n_seq*0.65:
            ee += 0.5
            response = input (f"The number of seqs with an Expected Error <= {ee-0.5}: {good_seq}/{n_seq}\nDo you want to:\n- [Press 0] keep them \n- [Press 1] try with an EE<={ee}\nDigit your response: ").strip()
            if response == "0":
                break
            elif response == "1":
                err_list = exp_errq(qscore_num_list, max_err=ee)
                good_seq = err_list.count(1)
            else:
                while response not in ["0", "1"]:
                    response = input ("Press only '0' or '1': ").strip()
                if response == "0":
                    break
                elif response == "1":
                    err_list = exp_errq(qscore_num_list, max_err=ee)
                    good_seq = err_list.count(1)
        else:
            break

    return err_list


# Questa funzione esclude tutte le sequenze che presentano un errore maggiore di 2
#   e tutte le sequeze duplicate con un quality score più basso

def exclude_seq (acc_num_list, seq_list, comment_list, qscore_num_list, err_list):
    def mean_qscore (qscore_num_list):
        n_seq = len(qscore_num_list)
        mean_qscore_list = [0 for _ in range (n_seq)]
        for s in range (n_seq):
            mean_qscore_list[s] = sum(qscore_num_list[s])/len(qscore_num_list[s])
        return mean_qscore_list
    
    n_seq = len(seq_list)
    valid_seq = []
    for s in range (n_seq):
        if err_list[s]==1:
            valid_seq.append([acc_num_list[s], seq_list[s], comment_list[s], qscore_num_list[s]])
            
    n_seq = len(valid_seq)
    dupl_idx_list = [None for _ in range(n_seq)]
    for s1 in range (n_seq-1):
        if dupl_idx_list [s1] is None: # se la seq1 in esame NON è stata assegnata si procede
            dupl_idx_list[s1] = s1
            for s2 in range (s1+1, n_seq):
                if dupl_idx_list [s2] is None: # Se la seq2 in esame NON è stata assegnata si procede
                    if valid_seq[s1][1]==valid_seq[s2][1]:
                        dupl_idx_list [s2] = s1
    print ("\tduplicated sequences found.")
    
    qscore_num_list = [qsnum[3] for qsnum in valid_seq]
    mean_qs_list = mean_qscore(qscore_num_list)
    print ("\tmean quality score calculated.")

    dupl_idx_set = {} # creo un dizionario per raggruppare sotto un unica chiave (idx della seq)
                      # tutte le sequenze uguali riferendomi ad esse tramite il loro idx di lista
    for i in range(n_seq):
        kidx = dupl_idx_list[i]
        dupl_idx_set.setdefault(kidx, []).append(i)
    
    valid_idx = []
    for kidx, indices in dupl_idx_set.items():
        best_idx = max(indices, key=lambda x: mean_qs_list[x]) # trova la seq con qscore medio maggiore
        valid_idx.append(best_idx)
    
    valid_idx.sort()
    final_valid_seq = [valid_seq[idx] for idx in valid_idx]
    return final_valid_seq
  

# Funzione che permette il merging delle PE reads.
# restituisce una lista di liste: 
#   [0] merged sequence
#   [1] idx che indica il punto in cui inizia la sovrapposizione nella read 1
#   [2] index seq 1 (indice della lista di partenza 'seq1_list')
#   [3] index seq 2
def merge_pe_process(F1_seq_part, F2rc_seq_list, queue:Queue):
    def merging(seq1: str, seq2: str, min_overlap=20):
        len_seq1 = len(seq1)
        for i in range(len_seq1 - min_overlap):
            overlap = seq1[i:]
            if seq2.startswith(overlap):
                return f"{seq1[:i]}{seq2}", i
        return None

    merged_sequences_info = []
    for i, seq1 in enumerate(F1_seq_part):  
        for j, seq2 in enumerate(F2rc_seq_list):
            merg_seq = merging(seq1, seq2)
            if merg_seq:
                merged_sequences_info.append([merg_seq[0], merg_seq[1], i, j])

    queue.put(merged_sequences_info)

@speedtest
def merge_pe_hpc(F1_seq_list, F2rc_seq_list, CPU=cpu_count()-1):
    queue = Queue()
    processes = []
    chunk_size = len(F1_seq_list)//CPU  

    for i in range(CPU):
        start_idx = i * chunk_size
        end_idx = (i + 1) * chunk_size if i < CPU - 1 else len(F1_seq_list)
        process = Process(target=merge_pe_process, args=(F1_seq_list[start_idx:end_idx], F2rc_seq_list, queue))
        processes.append(process)
        process.start()

    # Raccolta dei risultati
    merged_sequences_info = []
    for _ in range(CPU):
        merged_sequences_info.extend(queue.get())

    for process in processes:
        process.join()

    return merged_sequences_info


def combo_seq(seq:str):
    base_options = [IUPAC_CODES[base] for base in seq]
    combo_seq_list = [''.join(c) for c in product(*base_options)] # itertools.product genera tutte le combinazioni
    return list(set(combo_seq_list))


def weight_jcrd_idx (seq1:list, seq2:list, len_kmer:int): 
    def weight_kmers_set (seq:str, len_kmer:int, alpha=0.25, max_n=len_kmer//2):
        weighted_kmers = {}
        seq_len = len(seq)

        for i in range(seq_len - len_kmer + 1):
            new_kmer = seq[i:i+len_kmer]
            num_N = new_kmer.count("N")
            if num_N >= max_n:
                continue
            weight = alpha**num_N
            weighted_kmers[new_kmer] = max(weighted_kmers.get(new_kmer, 0), weight)
        return weighted_kmers
    
    kmer_s1 = weight_kmers_set(seq1, len_kmer) # merge seq
    kmer_s2 = weight_kmers_set(seq2, len_kmer) # refer seq

    union_keys = set(kmer_s1.keys()) | set(kmer_s2.keys())

    num = 0.0  # numeratore
    den = 0.0  # denominatore

    for x in union_keys:
        wA = kmer_s1.get(x, 0)
        wB = kmer_s2.get(x, 0)
        num += min(wA, wB)
        den += max(wA, wB)

    if den == 0:
        return 0.0
    else:
        return num/den

def jcrd_idx (seq1:list, seq2:list, len_kmer:int): 
    def set_of_kmers (seq:str, len_kmer:int):
        len_seq = len(seq)
        set_kmer = {}
        for i in range (len_seq-len_kmer+1):
            new_kmer = seq[i:i+len_kmer]
            if new_kmer in set_kmer:
                set_kmer[new_kmer] += 1
            else: 
                set_kmer[new_kmer] = 1
        return set_kmer
    
    set_kmer_s1 = set(set_of_kmers(seq1, len_kmer).keys())
    set_kmer_s2 = set(set_of_kmers(seq2, len_kmer).keys())

    intersection = set_kmer_s1.intersection(set_kmer_s2)
    union = set_kmer_s1.union(set_kmer_s2)

    if len(union) > 0:
        jaccard_index = len(intersection) / len(union)
        return jaccard_index
    else:
        return 0
    


def file_init(FOLD_PATH: str, n_pool_cpu: int):
    folder_path = f"{FOLD_PATH}cluster_file"
    if not os.path.isdir(folder_path): # controllo se esiste il path, sennò lo creo
        os.mkdir(folder_path)
        for file_idx in range(n_pool_cpu): # creo i file 
            with open(f"{folder_path}/OTU_clustered_seq_{file_idx}.txt", "wb") as file:
                pass

def OTU_cluster(FOLD_PATH, ref_seq_fa, merged_sequences_info, F1F2_acc_num_list, file_idx, N_CPU, idx_range):
    n_merge_seq = len(merged_sequences_info)
    n_ref_seq = len(ref_seq_fa)
    start_idx = file_idx * idx_range
    end_idx = start_idx + idx_range if file_idx < N_CPU else n_ref_seq
    method = ""
    with open(f"{FOLD_PATH}cluster_file/OTU_clustered_seq_{file_idx}.txt", "wb") as file:
        for rf in range(start_idx, end_idx):  
            ref_seq = f"{ref_seq_fa[0][rf]}\n".encode('utf-8')
            file.write(ref_seq)
            for ms in range(n_merge_seq): 
                for kmer_len in range (KMER_LEN_MIN, KMER_LEN_MAX+1):
                    if "N" in ref_seq_fa[1][rf]:  
                        jaccard_idx = weight_jcrd_idx(merged_sequences_info[ms][0], ref_seq_fa[0][rf], kmer_len)
                        method = "W"
                    else:
                        method = "N"
                        jaccard_idx = jcrd_idx(merged_sequences_info[ms][0], ref_seq_fa[0][rf], kmer_len)
                    if jaccard_idx > 0.0:
                        id_seq1 = F1F2_acc_num_list[0][merged_sequences_info[ms][2]]
                        id_seq2 = F1F2_acc_num_list[0][merged_sequences_info[ms][3]]
                        result = f"[k={kmer_len}] [{method}ji={jaccard_idx}] [{id_seq1}||{id_seq2}]\t\t{merged_sequences_info[ms][0]}\n".encode('utf-8')
                        file.write(result)
                        break


@speedtest
def filter_file1(queue, F1_acc_num_list, F1_seq_list, F1_comment_list, F1_qscore_num_list, F1_err_list):
    print("START FILTERING file 1 sequences...")
    F1_valid_seq = exclude_seq(F1_acc_num_list, F1_seq_list, F1_comment_list, F1_qscore_num_list, F1_err_list)
    queue.put(F1_valid_seq if F1_valid_seq else [])  # Assicura sempre un valore
    print("\tNumber of valid F1 seq:", len(F1_valid_seq))

@speedtest
def filter_file2(queue, F2_acc_num_list, F2_seq_list, F2_comment_list, F2_qscore_num_list, F2_err_list):
    print("\nSTART FILTERING file 2 sequences...")
    F2_valid_seq = exclude_seq(F2_acc_num_list, F2_seq_list, F2_comment_list, F2_qscore_num_list, F2_err_list)
    queue.put(F2_valid_seq if F2_valid_seq else [])  # Assicura sempre un valore
    print("\tNumber of valid F2 seq:", len(F2_valid_seq))

if __name__ == "__main__":
    #F1 = f"{FOLD_PATH}BAM_1012_136_R1.fastq.gz"
    #F2 = f"{FOLD_PATH}BAM_1012_136_R2.fastq.gz"
    F1 = f"{FOLD_PATH}merging_test1.fq.gz"
    F2 = f"{FOLD_PATH}merging_test2.fq.gz"
    REF_SEQ_FILE = f"{FOLD_PATH}ref.fa"

    F1_acc_num_list, F1_seq_list, F1_comment_list, F1_qscore_list = fq_pars_hpc(F1)
    F1_qscore_num_list = conv_qs(F1_qscore_list)
    F1_err_list = check_ee(F1_qscore_num_list)
    print('\n')

    F2_acc_num_list, F2_seq_list, F2_comment_list, F2_qscore_list = fq_pars_hpc(F2)
    F2_qscore_num_list = conv_qs (F2_qscore_list)
    F2_err_list = check_ee(F2_qscore_num_list)
    print('\n')
            
    queue1 = Queue()
    queue2 = Queue()

    p1 = Process(target=filter_file1, args=(queue1, F1_acc_num_list, F1_seq_list, F1_comment_list, F1_qscore_num_list, F1_err_list))
    p2 = Process(target=filter_file2, args=(queue2, F2_acc_num_list, F2_seq_list, F2_comment_list, F2_qscore_num_list, F2_err_list))

    p1.start()
    p2.start()

    p1.join()
    print("1")
    
    p2.join()
    print("a")

    # Assicura che ci siano dati nelle code prima di prelevarli
    if not queue1.empty():
        F1_valid_seq = queue1.get()
    else:
        F1_valid_seq = []
    
    print("b")

    if not queue2.empty():
        F2_valid_seq = queue2.get()
    else:
        F2_valid_seq = []

    print("c")

    print("\nFiltering completato su entrambi i file.")
    print("Numero di sequenze filtrate F1:", len(F1_valid_seq))
    print("Numero di sequenze filtrate F2:", len(F2_valid_seq))



    """
    with open("C:/Users/andre/Desktop/my_project/combo_pf/prjct_gamma/F1_tab_F2.txt", "w") as file:
        for seq1, seq2 in zip(F1_seq_list, F2_seq_list):
            file.write(f"{seq1}\t{seq2}\n")
    """

    # MERGING DELLE SEQUENZE
    print("\nSTART MERGING OF P.E. READS")
    F1_seq_list = [seq[1] for seq in F1_valid_seq]
    F2_seq_list = [seq[1] for seq in F2_valid_seq]
    F2_rc_seq_list = rev_seqlist(compl_seqlist(F2_seq_list))
    print("\treverse complementary F2 sequences calculated.")
    #   start
    merged_sequences_info = merge_pe_hpc(F1_seq_list, F2_rc_seq_list)
    num_merged_seq = len(merged_sequences_info)
    print(f"Numero di reads mergiate: {num_merged_seq}")


    # PARSING FASTA SEQ DI RIFERIMENTO
    ref_seq_fa = fa_pars_1cpu(REF_SEQ_FILE) # [0]id || [1]seq
    n_ref_seq = len(ref_seq_fa[0])


    # CLUSTERING
    # Stampa su file:
    #   -- id sequenza di riferimento --
    #   -- elenco seq mergiate -- {id_merg_seq (idx_seq1.idx_seq2): seq_mergiata)}
    N_CPU = cpu_count()-1
    idx_range = n_ref_seq//N_CPU
    file_init(FOLD_PATH, N_CPU)

    F1_acc_num_list = [an[0] for an in F1_valid_seq]
    F2_acc_num_list = [an[0] for an in F2_valid_seq]
    F1F2_acc_num_list = [F1_acc_num_list, F2_acc_num_list]
    
    pool = Pool(processes=N_CPU)
    for file_idx in range(N_CPU):
        print(f"Start CLUSTER {file_idx}")
        pool.apply_async(OTU_cluster, args=(FOLD_PATH, ref_seq_fa, merged_sequences_info, F1F2_acc_num_list, file_idx, N_CPU, idx_range))
    pool.close()
    pool.join()
    

    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as ref_seq_file:
        json.dump(ref_seq_fa, ref_seq_file)
        ref_seq_file_path = ref_seq_file.name.replace("\\", "/")

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as merged_seq_file:
        json.dump(merged_sequences_info, merged_seq_file)
        merged_seq_file_path = merged_seq_file.name.replace("\\", "/")

    # Avvio dei processi in nuovi terminali visibili
    for file_idx in range(N_CPU):
        print(f"Start CLUSTER {file_idx}")

        if os.name == "nt":
            cmd = (
                f"start cmd.exe /k \"cd /d {FOLD_PATH} & "
                f"{sys.executable} -c \"import sys; sys.path.append(r'{FOLD_PATH}'); import otu_clustering, json; "
                f"ref_seq_fa = json.load(open(r'{ref_seq_file_path}')); "
                f"merged_sequences_info = json.load(open(r'{merged_seq_file_path}')); "
                f"otu_clustering.OTU_cluster(r'{FOLD_PATH}', ref_seq_fa, merged_sequences_info, {file_idx}, {N_CPU}, {idx_range})\" & pause\""
            )
            subprocess.Popen(cmd, shell=True)

    input("Premi INVIO per terminare...")
    """