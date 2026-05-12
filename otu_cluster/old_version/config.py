from multiprocessing import cpu_count
from functools import wraps
import time

N_CPU = cpu_count()-1


KMER_LEN_MIN = 8
KMER_LEN_MAX = 12
OFF_33 = 33
OFF_64 = 64
CHECK_OFF_SCORE_33_63 = {chr(i): i for i in range(OFF_33, OFF_64-1)}
CHECK_OFF_SCORE_83_103 = {chr(i): i for i in range(OFF_33+52, OFF_64+52)}
Q_SCORE_DICT_33 = {chr(i): i for i in range(OFF_33, OFF_33+52)}
Q_SCORE_DICT_64 = {chr(i): i for i in range(OFF_64, OFF_64+52)}
VALID_BASES = "ATCGN"
COMPLEMENTARY = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

IUPAC_CODES = {'A':['A'], 'T':['T'], 'C':['C'], 'G':['G'],
               'R':['A', 'G'],
               'Y':['C', 'T'],
               'S':['G', 'C'],
               'W':['A', 'T'],
               'K':['G', 'T'],
               'M':['A', 'C'],
               'B':['C', 'G', 'T'],
               'D':['A', 'G', 'T'],
               'H':['A', 'C', 'T'],
               'V':['A', 'C', 'G'], 
               'N':['A', 'T', 'C', 'G']}


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


def eta(batch_size=10_000):
    """
    Decoratore per tracciare il tempo di esecuzione di una funzione con ETA.
    
    Args:
    - batch_size (int): Numero di iterazioni tra ogni calcolo di ETA.

    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            iterator = func(*args, **kwargs)  # La funzione deve restituire un generatore
            
            num_combo_ref = kwargs.get("num_combo_ref", None) or args[1]  
            count = 0
            start = time.time()
            
            for item in iterator:
                yield item  # Restituisce il valore al chiamante
                count += 1
                
                if count % batch_size == 0:
                    elapsed = time.time() - start
                    time_per_iter = elapsed / batch_size
                    eta = num_combo_ref*time_per_iter/batch_size  # Stima del tempo rimanente

                    hours = int(eta / 3600)
                    minutes = int((eta % 3600) / 60)
                    seconds = round((eta % 60), 3)

                    eta_string = f"{hours:02}:{minutes:02}:{seconds:02.2f}"
                    print(f"| {count} || ETpS: {time_per_iter:.2f} || ETA: {eta_string} [hh:mm:ss]")

                    start = time.time()  # Resetta il timer per il prossimo batch
            
            total_time = time.time() - start_time
            hours = int(total_time / 3600)
            minutes = int((total_time % 3600) / 60)
            seconds = round((total_time % 60), 3)
            total_time = f"{hours:02}:{minutes:02}:{seconds:05.2f}"
            print(f"Time of checking all the combo: {total_time} [hh:mm:ss]")

        return wrapper
    return decorator