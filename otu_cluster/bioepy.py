import os, sys, time
from functools import wraps
import pandas as pd
import dask.bag as db
import regex as re
from multiprocessing import cpu_count, Pool
CPU = cpu_count()-1

VALID_BASES = 'ATCGN'
COMPLEMENTARY = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
RNA_TRANSCRIPT = {'A':'A', 'T':'U', 'C':'C', 'G':'G'}
IUPAC_CODE = {
    'A':['A'], 'T':['T'], 'C':['C'], 'G':['G'], 'N':['A', 'T', 'C', 'G'],
    'R':['A', 'G'], 'Y':['C', 'T'], 'S':['G', 'C'], 'W':['A', 'T'],
    'K':['G', 'T'], 'M':['A', 'C'], 'B':['C', 'G', 'T'],
    'D':['A', 'G', 'T'], 'H':['A', 'C', 'T'], 'V':['A', 'C', 'G']}

CODON_DIZ = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',

    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',

    'TAT':'Y', 'TAC':'Y', 'TAA':'s', 'TAG':'s',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',

    'TGT':'C', 'TGC':'C', 'TGA':'s', 'TGG':'W',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
}

Q_SCORE_DICT_33 = {chr(i):i for i in range(33, 85)}  # Phred-33
Q_SCORE_DICT_64 = {chr(i):i for i in range(64, 116)}  # Phred-64
CHECK_OFF_SCORE_33_63 = set(chr(i) for i in range(33, 63))
CHECK_OFF_SCORE_83_103 = set(chr(i) for i in range(85, 103))

# @decoratore
def speedtest(func):
    """
    sped test decorator for any function
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        hours = int((elapsed_time) / 3600)
        minutes = int((elapsed_time - hours*3600) / 60)
        seconds = round((elapsed_time - hours*3600 - minutes*60), 3)
        if minutes==0 and hours==0:
            print(f"Exec time of '{func.__name__}': {seconds:02} s")
        else:
            print(f"Total Time: {hours:02}:{minutes:02}:{seconds:02}[hh:mm:ss]\n\n")
        #print(f"Execution time of '{func.__name__}': {elapsed_time:.6f} seconds")
        return result
    return wrapper


def check_fastq_path (fold_path:str, file_name:str) -> tuple[str, str]:
    """
    Check the file existance and if it is a fastq file.
    Allow to change the path of file in case of syntax error.
    -------------------------------------------------------------------
    Args:
        fold_path (str): Fold path which contains all files needed
        file_name (str): The name of the fastq file
    Return:
        file_path (str): The (new) whole file path
        file_name (str): The (new) file name
    """

    file_path = f'{fold_path}{file_name}'

    if os.path.isfile(file_path):
        if file_path.endswith('.fastq','.fq','.fastq.gz','fq.gz','.fastq.bz2','.fq.bz2',
                                '.fastq.xz', '.fq.xz', '.fastq.zip', 'fq.zip'):
            return file_path, file_name
        else:
            while True:
                file_name = input ('The file selected IS NOT a fastq file. '
                                   'Please past here the correct name of your fastq file: ').strip()
                file_path = f'{fold_path}{file_name}'
                if file_path.endswith('.fastq','.fq','.fastq.gz','fq.gz','.fastq.bz2','.fq.bz2',
                                '.fastq.xz', '.fq.xz', '.fastq.zip', 'fq.zip'):
                    return file_path, file_name
                print()
    else:
        while True:
            file_name = input ('The file selected NOT exists. '
                               'Please past here the name of your fastq file: ').strip()
            file_path = f'{fold_path}{file_name}'
            if file_path.endswith('.fastq','.fq','.fastq.gz','fq.gz','.fastq.bz2','.fq.bz2',
                                '.fastq.xz', '.fq.xz', '.fastq.zip', 'fq.zip'):
                return file_path, file_name
            print()

def check_fasta_path (fold_path:str, file_name:str) -> tuple[str, bool]:
    """
    Check the file existance and if it is a fasta/q file.
    Allow to change the path of file in case of syntax error.
    ----------------------------------------------------------
    Args:
        fold_path (str): Fold path which contains all files needed
        file_name (str): The name of the fastq file
    Return:
        file_path (str): The (new) whole file path
        file_name (str): The (new) file name
    """

    file_path = f'{fold_path}{file_name}'

    if os.path.isfile(file_path):
        if file_path.endswith('.fasta','.fa','.fasta.gz','fa.gz','.fasta.bz2','.fa.bz2',
                                '.fasta.xz', '.fa.xz', '.fasta.zip', 'fa.zip'):
            return file_path, file_name
        else:
            while True:
                file_name = input ('The file selected IS NOT a fastq file. '
                                   'Please past here the correct name of your fasta file: ').strip()
                if file_path.endswith('.fasta','.fa','.fasta.gz','fa.gz','.fasta.bz2','.fa.bz2',
                                        '.fasta.xz', '.fa.xz', '.fasta.zip', 'fa.zip'):
                    return file_path, file_name
                print()
    else:
        while True:
            file_name = input ('The file selected IS NOT a fastq file. '
                                   'Please past here the correct name of your fasta file: ').strip()
            if file_path.endswith('.fasta','.fa','.fasta.gz','fa.gz','.fasta.bz2','.fa.bz2',
                                    '.fasta.xz', '.fa.xz', '.fasta.zip', 'fa.zip'):
                return file_path, file_name
            print()


def is_valid_seq(sequence: str) -> bool:
    """
    Checks if a DNA sequence is valid.
    A sequence is considered valid if all characters in it belong 
    to the set of valid nucleotide bases (VALID_BASES).
    --------------------------------------------------------------
    Args:
        sequence (str): The DNA sequence to verify.
    Returns:
        bool: True if the sequence is valid, False otherwise.
    """
    return all(char in VALID_BASES for char in sequence)

def is_valid_qs(qscore:str, phred:int = 33) -> bool:
    score_dict = Q_SCORE_DICT_33 if phred == 33 else Q_SCORE_DICT_64
    return all(char in score_dict for char in qscore)


def rev_seqlist(seq_list: list[str]) -> list[str]:
    """
    Reverses the order of nucleotide bases in each sequence of a list.
    -------------------------------------------------------------------
    Args:
        seq_list (list[str]): List of DNA sequences.
    Returns:
        list[str]: List of reversed sequences.
    """
    return [''.join(reversed(seq)) for seq in seq_list]

def compl_seqlist(seq_list: list[str]) -> list[str]:
    """
    Computes the complementary sequence for each sequence in the list.
    The function uses a complementary base dictionary (COMPLEMENTARY) 
    to replace each base with its complement.
    -------------------------------------------------------------------
    Args:
        seq_list (list[str]): List of DNA sequences.
    Returns:
        list[str]: List of complementary sequences.
    """
    return [''.join(COMPLEMENTARY.get(base, 'N') for base in seq) for seq in seq_list]


"""
def weight_jcrd_idx(seq1:str, seq2:str, len_kmer:int) -> float:
    def weight_kmers(seq:str, len_kmer:int, alpha:float=0.25, max_n:int=2) -> dict[str, float]:
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
"""

def weight_jcrd_idx(seq1: str, seq2: str, len_kmer: int) -> float:
    def weight_kmers(seq: str, len_kmer: int, alpha: float = 0.25, max_n: int = 2) -> dict[str, float]:
        """
        Generate weighted k-mers for a sequence, assigning a weight based on the presence of 'N' characters.
        """
        kmers = {}
        for i in range(len(seq) - len_kmer + 1):
            kmer = seq[i:i + len_kmer]
            n_count = kmer.count('N')
            if n_count < max_n:
                kmers[kmer] = max(kmers.get(kmer, 0), alpha ** n_count)  # Assign weight based on 'N' presence
        return kmers

    k1, k2 = weight_kmers(seq1, len_kmer), weight_kmers(seq2, len_kmer)
    intersection_keys = set(k1.keys()) & set(k2.keys())  # Intersezione tra i k-mer delle due sequenze
    num = sum(min(k1.get(k, 0), k2.get(k, 0)) for k in intersection_keys)  # Somma dei minimi (peso dei match)
    den = sum(k1.get(k, 0) for k in k1.keys())  # Normalizzazione rispetto ai k-mer di seq1

    return num/den if den > 0 else 0.0 

def jcrd_idx(seq1:str, seq2:str, len_kmer:int) -> float:
    k1 = set(seq1[i:i+len_kmer] for i in range(len(seq1)-len_kmer+1))
    k2 = set(seq2[i:i+len_kmer] for i in range(len(seq2)-len_kmer+1))
    intersection = k1 & k2
    #union = k1 | k2
    return len(intersection)/len(k1) 

class ParsingFasta:
    """
    FastA file parsing class.

    Class attributes:
    - file_type (str): 'FastA'.

    Istance attributes:
    - ids (list[str]): List of sequence identifiers.
    - seqs (list[str]): List of sequences.
    """

    file_type = 'FastA'

    def __init__(self, file_path:str, file_name:str):
        self.file_path = file_path
        self.file_name = file_name
        self.ids = []
        self.seqs = []

        self.parse_fasta()
        if len(self.ids)==1:
            self.ids = self.ids[0]
            self.seqs = self.seqs[0]


    def parse_fasta(self):
        """
        Parsing of a FastA file.
        """
        try:
            fasta_data = list(db.read_text(self.file_path))
            if not fasta_data or not fasta_data[0].startswith('>'):
                raise ValueError(f"Error: The FastA file '{self.file_path}' has not a valid header.")
                
            i = 0
            while i < len(fasta_data):
                line = fasta_data[i].strip()
                if line.startswith('>'):
                    self.ids.append(line.lstrip('>'))
                    seq = []
                    i += 1
                    while i < len(fasta_data) and not fasta_data[i].startswith('>'):
                        seq.append(fasta_data[i].strip().upper())
                        i += 1
                    self.seqs.append(''.join(seq))
                else:
                    i += 1
        except Exception as e:
            print(f"Error while parsing file '{self.file_path}': {e}")
            sys.exit(1)


class ParsingFastq:
    """
    FastQ file parsing class.

    Class attributes:
        file_type (str): 'FastQ'
    """

    file_type = 'FastQ'

    def __init__(self, file_path:str, file_name:str):
        """
        Parsing of a FastA file.
        ---------------------------------------------------------------
        Args:
            file_path (str): Percorso completo del file FastQ.
            file_name (str): Nome del file FastQ.

        Istance attributes:
            ids (list[str]): List of sequence identifiers.
            seqs (list[str]): List of sequences.
            comments (list[str]): List of comments associated with sequences.
            qss (list[str]): List of quality scores.
            qss_num (list[list[int]]): Quality scores converted to numeric.
            exp_errs (list[bool]): List of 
        """
        self.file_path = file_path
        self.file_name = file_name

        self.ids = []
        self.seqs = []
        self.comments = []
        self.qss = []
        self.qss_num = []
        self.exp_errs = []
        self.valid_seqs = []

        self.pars_fastq()
        

    def pars_fastq(self):
        """
        Parsing of a FastA file.
        Parallelized process: it uses all the cpu available minus one (left one for OS subprocess)
        """
        fastq_data = list(db.read_text(self.file_path))
        if not fastq_data or len(fastq_data) % 4 != 0:
            raise ValueError(f"Errore: Il file FASTQ '{self.file_path}' ha record incompleti.")

        num_cpu = max(1, cpu_count()-1)  # Usa tutti i core meno uno
        chunk_size = (len(fastq_data)//(num_cpu*4))*4 or 4
        chunks = [fastq_data[i:i+chunk_size] for i in range(0, len(fastq_data), chunk_size)]

        with Pool(processes=num_cpu) as pool:
            results = pool.starmap(self.pars_fastq_chunk, [(chunk, i*chunk_size, self.file_name) for i, chunk in enumerate(chunks)])

        for res in results:
            self.ids.extend(res[0])
            self.seqs.extend(res[1])
            self.comments.extend(res[2])
            self.qss.extend(res[3])

        print (f'File name: {self.file_name}')
        print (f'Number of sequence: {len(self.seqs)}')
    
        self.qss_num = self.convert_qss(self.qss)
        self.exp_errs = self.check_ee(self.qss_num)
        self.valid_seqs = self.exclude_seq(self.ids, self.seqs, self.qss_num, self.exp_errs)
        print()

    def pars_fastq_chunk(self, chunk, start_idx, file_name):
        """
        Chunk worker of the 'pars_fastq' function
        """
        ids, seqs, comments, qss = [], [], [], []
        try:
            for i in range(0, len(chunk), 4):
                try:
                    if not chunk[i].startswith('@'):
                        raise ValueError(f'Id Error: NOT valid ID (file: {file_name} || line: {start_idx+i})')
                    ids.append(chunk[i].strip()[1:-2])

                    seq = chunk[i+1].strip()
                    if not is_valid_seq(seq):
                        print(f'Quality Score Error: NOT valid char in sequence (file: {file_name} || line: {start_idx+i+1})')
                    seqs.append(seq)

                    if not chunk[i+2].startswith('+'):
                        raise ValueError(f'Comment Error: NOT valid COMMENT (file: {file_name} || line: {start_idx+i+2})')
                    comments.append(chunk[i+2].strip())

                    qs = chunk[i+3].strip()
                    if len(qs) != len(seq):
                        raise ValueError(f'Quality Score Error: the length of quality score does not correspond at the correspective sequence (file: {file_name} || line: {start_idx+i+3})')
                    if not (is_valid_qs(qs, 33) or is_valid_qs(qs, 64)):
                        print(f'Quality Score Error: NOT valid char in quality score (file: {file_name} || line: {start_idx+i+3})')
                    qss.append(qs)

                except IndexError:
                    raise ValueError(f'Parsing Error: NOT valid record (line: {start_idx+i})')
          
        except Exception as e:
            print (f"Error while parsing file '{self.file_name}': {e}")
            sys.exit(1)

        return ids, seqs, comments, qss
    

    def convert_qss(self, qscore_list:list[str]) -> list[list[int]]:
        def detect_phred(qss:str) -> int:
            combined = ''.join(qss)
            if any(c in CHECK_OFF_SCORE_33_63 for c in combined) and not any(c in CHECK_OFF_SCORE_83_103 for c in combined):
                print('Detected sequencing method: Sanger (Phred-33)')
                return 33
            print('Detected sequencing method: Illumina (Phred-64)')
            return 64

        phred = detect_phred(qscore_list)
        score_dict = Q_SCORE_DICT_33 if phred==33 else Q_SCORE_DICT_64
        return [[score_dict[char]-phred for char in qs] for qs in qscore_list]
    
    def check_ee(self, qss_num, max_err:float=2.0) -> list[int]:
        """
        Evaluate the quality of sequences based on the Expected Error (EE)
        and allows filtering out the sequence with high error.

        If the EE is less than or equal to 'max_err', the sequence is considered valid.
        The user can increase the 'max_err' threshold if less than 65% of the sequences
        satisfy the initial criterion.
        -------------------------------------------------------------------------------------
        Arguments:
            qss_num (list[list[int]]): List of lists containing the numerical quality scores for each sequence.
            max_err (float, optional): Maximum Expected Error threshold. Default is 2.0.

        Returns:
            list[int]: Binary list where:
                - '1' indicates that the sequence is accepted (EE ≤ max_err).
                - '0' indicates that the sequence is discarded."
        """
        ee_list = [1 if sum([10**(q/-10) for q in qs]) <= max_err else 0 for qs in qss_num]
        good_seq = sum(ee_list)
        total_seq = len(ee_list)
        curr_max_err = max_err
        
        while good_seq <= total_seq*0.85:
            print(f'Sequences with EE <= {curr_max_err}: {good_seq}/{total_seq}')
            try:
                choice = input('Keep current sequences (0) or increase EE threshold (1)? ').strip()
                if choice == '0':
                    break
                elif choice == '1':
                    curr_max_err += 0.5
                    ee_list = [1 if sum(10**(q/-10) for q in qs) <= curr_max_err else 0 for qs in qss_num]
                    good_seq = sum(ee_list)
                else:
                    print('Invalid input, please enter 0 or 1')
            except KeyboardInterrupt:
                print('Process interrupted by user')
                break
        return ee_list
    
    def exclude_seq(self, ids, seqs, qss_num, exp_errs) -> list[list]:
        """
        Filters the sequence based on the Expected Error (EE) and removes duplicates by choosing
        the sequence with the best average quality score.

        Filtering occurs in two steps:
            Only sequences with 'exp_errs[i]==1' (valid sequence) are selected.
            If there are duplicates (same sequence), the one with the highest average quality score is retained.
        --------------------------------------------------------------------------------------------------------
        Arguments:
            ids (list[str]): List of sequence identifiers.
            seqs (list[str]): List of DNA/RNA sequences.
            qss_num (list[list[int]]): List of lists containing numerical quality scores for each sequence.
            exp_errs (list[int]): Binary list indicating whether a sequence is valid (1) or not (0).
        Returns:
            list[list[str]]: List of filtered and unique sequences [[id_1, seq_1], [id_2, seq_2], ...]
        """
        
        num_seqs = len(seqs)
        mean_qs = [sum(qs) / len(qs) for qs in qss_num]
        valid_seq = [[ids[i], seqs[i]] for i in range(num_seqs) if exp_errs[i] == 1]
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
    

class MergingPE:
    def __init__ (self, valid_seqs_f1:list[list[str]], valid_seqs_f2:list[list[str]]):
        """
        """
        self.valid_seqs_f1 = {seq[0]: seq[1] for seq in valid_seqs_f1}
        self.valid_rc_seqs_f2 = {seq[0]: compl_seqlist(rev_seqlist([seq[1]]))[0] for seq in valid_seqs_f2}

        self.merged_sequences = self.merge_pe (self.valid_seqs_f1, self.valid_rc_seqs_f2)
        print (f'Number of Merged Sequence: {len(self.merged_sequences)}')

    def merge_pe(self, valid_seq_f1:dict, rc_valid_seq_f2:dict, min_overlap:int=8, max_err:int=3) -> list[list]:
        merged_sequences = []
        counter = 0

        # itero sugli id (chiavi del dizionario)
        for seq_id, seq1 in valid_seq_f1.items():
            if seq_id in rc_valid_seq_f2:  # Only if the accession number are equal, try to merge
                counter += 1
                seq2 = rc_valid_seq_f2[seq_id]

                for overlap in range (len(seq1), min_overlap-1, -1):
                    pattern = re.compile(f'({seq1[-overlap:]}){{e<={max_err}}}')
                    match = re.search(pattern, seq2)
                    if match:
                        overlap_start = seq2.find(match.group(1))  # Get position of overlap
                        merged_seq = seq1[:-overlap] + seq2[overlap_start:]  # Merge sequences
                        merged_sequences.append([seq_id, merged_seq, overlap_start])
                        break
        print ('counter id:', counter)

        return merged_sequences
    

class ClusterOTU:
    """
    Class to perform clustering of OTU (Operational Taxonomic Units) sequences.

    This class performs clustering of reference sequences with merged sequences based on
    the Jaccard index, either standard or weighted, over different k-mer length ranges.

    Attributes:
        kmer_min (int): Minimum length of the k-mer to consider for clustering.
        kmer_max (int): Maximum length of the k-mer to consider for clustering.

    Methods:
    - __init__: Initializes the class and starts clustering for each k-mer value in the specified range.
    - otu_clustering: Divides the work among the processors and starts parallel clustering.
    - otu_clustering_chunk: Performs clustering on a subset of the sequences.
    - otu_file_init: Initializes the folders and files to store the clustering results.
    """

    def __init__(self, fold_path:str, ref_seqs:list[list[str]], merg_seqs:list[list[str]], kmer_min:int=8, kmer_max:int=12):
        """
        Initializes OTU clustering on a range of k-mer lengths.

        Iterates over all k-mer values between `kmer_min` and `kmer_max`, creating a directory for each value
        and starting the clustering process.

        Args:
            fold_path (str): Path to the folder where to save the clustering results.
            ref_seqs (list[list[str]]): List of reference sequences with their identifiers.
            merg_seqs (list[list[str]]): List of merged sequences with their identifiers.
            kmer_min (int, optional): Minimum k-mer length for clustering (default: 8).
            kmer_max (int, optional): Maximum k-mer length for clustering (default: 12).
        """
        self.kmer_min = kmer_min
        self.kmer_max = kmer_max

        for kl in range (self.kmer_max, self.kmer_min-1, -1):
            print (f'Starting OTU clustering (kmer length: {kl})...', end='')
            self.otu_file_init(fold_path, kl)
            self.otu_clustering(fold_path, ref_seqs, merg_seqs, kl)
            

    def otu_clustering (self, fold_path, ref_seqs, merg_seqs, kmer_len):
        num_ref_seqs = len(ref_seqs)
        idx_range = num_ref_seqs//CPU + (num_ref_seqs)%CPU > 0
        
        with Pool(processes=CPU) as pool:
            pool.map(self.otu_clustering_chunk, [(fold_path, ref_seqs, merg_seqs, i, idx_range, kmer_len) for i in range(CPU)])
        print(' OTU clustering completed.')
        file_fold_path = f'{fold_path}cluster_file_k{kmer_len}'
        if os.path.isdir(file_fold_path) and not os.listdir(file_fold_path):
            os.rmdir(file_fold_path) 

    def otu_clustering_chunk(self, args:tuple[list[list[str]], list[list[str]], int, int, int]) -> None:
        fold_path, ref_seqs, merg_seqs, file_idx, idx_range, kmer_len = args
        start_idx = file_idx*idx_range
        end_idx = min(start_idx+idx_range, len(ref_seqs[0]))
        file_path = f'{fold_path}cluster_file_k{kmer_len}/OTU_clustered_seq_{file_idx}.txt'
        
        with open(file_path, 'w', encoding='utf-8') as f:
            for rf in range(start_idx, end_idx):
                flag = True # usefull to understand if the ID of ref seq it is already written or not 
                            # (NOT needed if not merge seq is clustered in that ref seq)
                for ms in merg_seqs:
                    method = 'W' if 'N' in ref_seqs[1][rf] else 'S' # (W)eighted Jaccard Idx || (S)tandard Jaccard Idx
                    j_idx = weight_jcrd_idx(ms[1], ref_seqs[rf][1], kmer_len) if method == 'W' else \
                            jcrd_idx(ms[1], ref_seqs[rf][1], kmer_len)
                    if j_idx >= 0.97:
                        if flag:
                            f.write(f'{ref_seqs[0][rf]}\n')
                            flag = False
                        f.write(f'[{method}ji={j_idx:.3f}] [{ms[0]}]\t\t{ms[1]}\n')
                f.write(f'\n\n')
        
        if os.path.isfile(file_path) and os.path.getsize(file_path) < 5:
            os.remove(file_path)
        
    
    def otu_file_init(self, fold_path:str, kmer_len:int) -> None:
        folder_path = f'{fold_path}cluster_file_k{kmer_len}'
        os.makedirs(folder_path, exist_ok=True)
        for file_idx in range(CPU):
            open(f'{folder_path}/OTU_clustered_seq_{file_idx}.txt', 'wb').close()



class Gene:
    """
    The Gene class represents a reference gene extracted from genomic annotation data.

    Attributes:
        assembly_name (str): Name of the genomic assembly.
        assembly_version (str): Version of the genomic assembly.
        gene_name (str): Name of the gene of interest.
        gene_id (str): Identifier for the gene (e.g. ENSG...).
        gene_info_df (pd.DataFrame): A dataframe containing annotation data (e.g., start, end, feature, etc.) 
                                     filtered for the specified gene.
    """
    assembly_name = "T2T"
    assembly_version = "CHM13.v2.0"

    def __init__(self, gene_name:str, gene_id:str, chr_id:str, chr_seq:str, annot_table_file_path:str):
        """
        Initializes a Gene object by reading annotation data relevant to the specified gene.

        Args:
            gene_name (str): The name of the gene of interest.
            gene_id (str): The identifier of the gene (converted to uppercase).
            chr_id (str): The chromosome identifier (possibly with position metadata).
            chr_seq (str): The full chromosome sequence in which the gene is located.
            annot_table_file_path (str): The path to the tab-delimited annotation file (e.g., derived from GFF/GTF).
        """
        self.gene_name = gene_name
        self.gene_id = gene_id.upper()

        self.gene_info_df = self.find_gene_seqs(annot_table_file_path, self.gene_id, chr_id, chr_seq)
    

    def find_gene_seqs(self, annot_table_file_path:str, gene_id:str, chr_id:str, chr_seq:str):
        """
        Extracts rows from the annotation file that match the provided gene_id, and appends 
        the corresponding nucleotide sequences to those rows.

        This method:
          1. Reads the annotation file into a dataframe.
          2. Filters rows whose gene identifier matches gene_id.
          3. Computes the sequence for each row based on [start, end] relative to chr_seq.
          4. Removes redundant columns (e.g., 'seqname', 'source', 'feature', 'start', 'end') and empty columns.

        Args:
            annot_table_file_path (str): Path to the tab-delimited annotation file.
            gene_id (str): The gene identifier to filter for.
            chr_id (str): Chromosome identifier string (with position info).
            chr_seq (str): The full chromosome sequence used to extract each gene segment.
        Returns:
            pd.DataFrame: A dataframe containing rows relevant to the specified gene_id, 
                          each with a 'sequence' column for the extracted DNA fragment.
        """
        df = pd.read_csv(annot_table_file_path, sep="\t")
        features = df.columns
        name_feat_geneid = [features[i] for i in range (df.shape[1]) \
                   if (features[i].__contains__('gene') and features[i].__contains__('id'))][0]
        
        # In gene_df sono presenti tutte le tuple che hanno 'gene_id' come valore per l'attributo gene_id
        #   e quindi contengono solo quelle tuple che riguardano il gene interessato
        gene_df = df[df[f'{name_feat_geneid}'].str.startswith(f'{gene_id}')].copy()

        # Aggiungo al dataframe <gene_df> l'attributo <sequences>, 
        #   grazie al quale ogni tupla avrà tutta la sequenza nucleotidica corrispondente
        #   (start-end) presa dal chr*.fa contenente la sequenza dell'intero cromosama 
        #   in cui il gene è contenuto.
        #
        # Per prima cosa ricavo gli indici di 'start' ed 'end' da gene_df
        start_list = df['start'].tolist()
        end_list = df['end'].tolist()
        
        # Dal chr_id ricavo lo start point assoluto
        chr_start = int((chr_id.split('chromosome:')[-1]).split(':')[2])

        # Poi genero la lista in cui sono riportate tutte le sequenze facendo riferimento
        seqs_list = [chr_seq[(start_list[i]-chr_start):(end_list[i]-chr_start)] for i in range (gene_df.shape[0])]
        gene_df['sequence'] = seqs_list
        gene_df = gene_df.drop(columns=['seqname', 'source', 'feature', 'start', 'end']) # rimuovo informazioni ridondanti o superflue
        gene_df.dropna(axis=1, how='all', inplace=True) # rimuovo dal dataframe tutte le colonne in cui tutti valori sono nulli

        return gene_df
    

    @classmethod
    def change_assembly_version(cls, new_assembly_version):
        """
        Class method to change the assembly version for all Gene/Transcript instances.

        Args:
            new_assembly_version (str): The new assembly version (e.g. 'CHM13.v3.0').
        """
        cls.assembly_version = new_assembly_version
    
    @classmethod
    def change_assembly_name(cls, new_assembly_name):
        """
        Class method to change the assembly name for all Gene/Transcript instances.

        Args:
            new_assembly_name (str): The new assembly name (e.g. 'T2T').
        """
        cls.assembly_version = new_assembly_name



class Transcript(Gene):
    """
    The Transcript class represents a transcript derived from a Gene. It inherits 
    basic annotation data (e.g., gene_info_df) from Gene, then adds methods for 
    performing transcription (DNA -> RNA) and translation (mRNA -> protein).

    Attributes:
        gene_trans_info_df (pd.DataFrame): A dataframe enriched with mRNA transcripts in different reading frames.
        gene_prot_info_df (pd.DataFrame): A dataframe enriched with protein translations (if 'gene_type' == 'protein_coding').
    """
    def __init__(self, gene_name, gene_id, chr_id, chr_seq, annot_table_file_path):
        """
        Initializes the Transcript class using the specified gene data, then 
        generates the transcripts and protein sequences for each reading frame.

        Args:
            gene_name (str): The name of the gene of interest.
            gene_id (str): The identifier of the gene (converted to uppercase).
            chr_id (str): The chromosome identifier string.
            chr_seq (str): The full chromosome sequence from which the transcripts are derived.
            annot_table_file_path (str): Path to the tab-delimited annotation file.
        """
        super().__init__(gene_name, gene_id, chr_id, chr_seq, annot_table_file_path)

        self.gene_trans_info_df = self.transcription(self.gene_info_df)
        self.gene_prot_info_df = self.translation(self.gene_trans_info_df)

    def transcription (self, gene_info_df:pd.DataFrame):
        """
        Performs the DNA -> RNA transcription for each row in the given dataframe, 
        generating transcripts in all three possible reading frames.

        Specifically, it replaces 'T' with 'U' (via a codon dictionary such as RNA_TRANSCRIPT) 
        and creates three new columns: 'transcript_seq0', 'transcript_seq1', 'transcript_seq2'.

        Args:
            gene_info_df (pd.DataFrame): A dataframe containing at least a 'sequence' column 
                                         with the gene's DNA.

        Returns:
            pd.DataFrame: The original dataframe with three additional columns for mRNA transcripts.
        """
        gene_trans_info_df = gene_info_df.copy()
        seqs_list = gene_info_df['sequence'].tolist()
            
        transcript_list0 = [''.join(RNA_TRANSCRIPT[seq[i]] for i in range(0, len(seq))) for seq in seqs_list]
        transcript_list1 = [''.join(RNA_TRANSCRIPT[seq[i]] for i in range(1, len(seq))) for seq in seqs_list]
        transcript_list2 = [''.join(RNA_TRANSCRIPT[seq[i]] for i in range(2, len(seq))) for seq in seqs_list]
        gene_trans_info_df["transcript_seq0"] = transcript_list0
        gene_trans_info_df["transcript_seq1"] = transcript_list1
        gene_trans_info_df["transcript_seq2"] = transcript_list2

        return gene_trans_info_df


    def translation (self, gene_info_df:pd.DataFrame):
        """
        Translates the DNA (or hypothetical RNA) sequences to proteins in each of the 
        three reading frames, but only for rows with 'gene_type' = 'protein_coding'.

        The resulting DataFrame includes new columns: 'protein_seq0', 'protein_seq1', and 'protein_seq2'.

        Args:
            gene_info_df (pd.DataFrame): A dataframe containing a 'sequence' column 
                                         and, optionally, a 'gene_type' column.

        Returns:
            pd.DataFrame: A copy of the dataframe containing the protein translation columns 
                          (protein_seq0, protein_seq1, protein_seq2) for all protein-coding genes.
        """
        gene_prot_info_df = gene_info_df[gene_info_df['gene_type'].str.startswith('protein_coding')].copy()
        seqs_list = gene_info_df['sequence'].tolist()
            
        protein_list0 = [''.join(CODON_DIZ[seq[i:i+3]] for i in range(0, len(seq)-(len(seq)%3), 3)) for seq in seqs_list]
        protein_list1 = [''.join(CODON_DIZ[seq[i:i+3]] for i in range(1, len(seq)-((len(seq)-1)%3), 3)) for seq in seqs_list]
        protein_list2 = [''.join(CODON_DIZ[seq[i:i+3]] for i in range(2, len(seq)-((len(seq)-2)%3), 3)) for seq in seqs_list]
        gene_prot_info_df["protein_seq0"] = protein_list0
        gene_prot_info_df["protein_seq1"] = protein_list1
        gene_prot_info_df["protein_seq2"] = protein_list2

        return gene_prot_info_df


