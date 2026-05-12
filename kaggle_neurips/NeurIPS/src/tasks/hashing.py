import math
import random
from typing import List, Set, Dict, Any, Union


class MinHash:
    """
    Classe MinHash per generare signature MinHash per insiemi.
    
    MinHash è una tecnica per stimare velocemente quanto sono simili due insiemi,
    utilizzando il coefficiente di similarità di Jaccard.
    """
    
    def __init__(self, num_perm: int = 128, seed: int = 1):
        """
        Inizializza una nuova istanza MinHash.
        
        Args:
            num_perm: Numero di permutazioni da usare (default: 128)
            seed: Seed per la generazione di numeri casuali (default: 1)
        """
        # Il numero primo più piccolo maggiore del più grande hash possibile (32 bit int)
        self.prime = 4294967311
        self.max_hash = 2**32 - 1
        self.num_perm = num_perm
        self.seed = seed
        
        # Inizializza i valori hash al valore massimo
        self.hashvalues = [self.max_hash] * self.num_perm
        
        # Inizializza le funzioni di permutazione
        self.perm_a, self.perm_b = self._init_permutations()
        
    def _init_permutations(self) -> tuple:
        """
        Inizializza le funzioni di permutazione per a & b.
        Non riutilizza alcun intero nella creazione delle funzioni.
        """
        random.seed(self.seed)
        used = set()
        perm_a = []
        perm_b = []
        
        # Genera permutazioni uniche per entrambi i parametri a e b
        for _ in range(self.num_perm):
            # Genera un numero casuale unico per a
            a = self._rand_int()
            while a in used:
                a = self._rand_int()
            perm_a.append(a)
            used.add(a)
            
            # Genera un numero casuale unico per b
            b = self._rand_int()
            while b in used:
                b = self._rand_int()
            perm_b.append(b)
            used.add(b)
            
        return perm_a, perm_b
    
    def _rand_int(self) -> int:
        """
        Genera un numero intero casuale >= 0 e <= max_hash.
        Usa lo stesso metodo del codice JavaScript originale.
        """
        # Simula il comportamento di Math.sin(seed++) * maxHash di JavaScript
        x = math.sin(self.seed) * self.max_hash
        self.seed += 1
        return int((x - int(x)) * self.max_hash)
    
    def _hash(self, string: str) -> int:
        """
        Funzione di hash che converte una stringa in un intero a 32 bit senza segno.
        Replica l'algoritmo di hash del codice JavaScript originale.
        
        Args:
            string: La stringa da hashare
            
        Returns:
            Hash a 32 bit della stringa
        """
        hash_val = 0
        if len(string) == 0:
            return hash_val + (2**32 // 2 - 1)
        
        for char in string:
            char_code = ord(char)
            hash_val = ((hash_val << 5) - hash_val) + char_code
            # Converte a intero a 32 bit
            hash_val = hash_val & 0xFFFFFFFF
            
        return hash_val + (2**32 // 2 - 1)
    
    def update(self, item: str) -> None:
        """
        Aggiorna i valori hash interni con un nuovo elemento.
        
        Args:
            item: Elemento da aggiungere al MinHash
        """
        item_hash = self._hash(item)
        
        for i in range(len(self.hashvalues)):
            a = self.perm_a[i]
            b = self.perm_b[i]
            hash_val = (a * item_hash + b) % self.prime
            
            if hash_val < self.hashvalues[i]:
                self.hashvalues[i] = hash_val
    
    def jaccard(self, other: 'MinHash') -> float:
        """
        Stima la similarità di Jaccard con un altro MinHash.
        
        Args:
            other: Altro oggetto MinHash per il confronto
            
        Returns:
            Stima della similarità di Jaccard (0-1)
            
        Raises:
            ValueError: Se i MinHash hanno parametri incompatibili
        """
        if len(self.hashvalues) != len(other.hashvalues):
            raise ValueError("I conteggi dei valori hash differiscono")
        if self.seed != other.seed:
            raise ValueError("I valori seed differiscono")
        
        shared = sum(1 for i in range(len(self.hashvalues)) 
                    if self.hashvalues[i] == other.hashvalues[i])
        
        return shared / len(self.hashvalues)
    
    def digest(self) -> List[int]:
        """
        Restituisce la signature MinHash come lista di valori hash.
        
        Returns:
            Lista dei valori hash della signature
        """
        return self.hashvalues.copy()


class LSHIndex:
    """
    Indice Locality Sensitive Hashing per trovare signature MinHash simili.
    
    Consente di trovare documenti simili in tempo sub-lineare raggruppando
    signature MinHash simili in bucket.
    """
    
    def __init__(self, band_size: int = 4):
        """
        Inizializza un nuovo indice LSH.
        
        Args:
            band_size: Dimensione delle bande per l'hashing (default: 4)
        """
        self.band_size = band_size
        self.index: Dict[str, List[str]] = {}
    
    def insert(self, key: str, minhash: MinHash) -> None:
        """
        Inserisce una signature MinHash nell'indice con una chiave associata.
        
        Args:
            key: Chiave identificativa per il documento
            minhash: Oggetto MinHash da indicizzare
        """
        hashbands = self._get_hashbands(minhash)
        
        for band in hashbands:
            band_key = str(band)
            if band_key not in self.index:
                self.index[band_key] = []
            self.index[band_key].append(key)
    
    def query(self, minhash: MinHash) -> List[str]:
        """
        Trova tutti i documenti simili a una query MinHash.
        
        Args:
            minhash: MinHash query
            
        Returns:
            Lista delle chiavi dei documenti candidati simili
        """
        matches = set()
        hashbands = self._get_hashbands(minhash)
        
        for band in hashbands:
            band_key = str(band)
            if band_key in self.index:
                matches.update(self.index[band_key])
        
        return list(matches)
    
    def _get_hashbands(self, minhash: MinHash) -> List[List[int]]:
        """
        Divide la signature MinHash in bande per l'indicizzazione LSH.
        
        Args:
            minhash: Oggetto MinHash
            
        Returns:
            Lista di bande (sottoliste di valori hash)
        """
        # Cache delle bande se già calcolate
        if hasattr(minhash, '_hashbands'):
            return minhash._hashbands
        
        hashbands = []
        hashvalues = minhash.hashvalues
        
        # Divide i valori hash in bande di dimensione band_size
        for i in range(0, len(hashvalues), self.band_size):
            band = hashvalues[i:i + self.band_size]
            if len(band) == self.band_size:  # Considera solo bande complete
                hashbands.append(band)
        
        # Cache per uso futuro
        minhash._hashbands = hashbands
        return hashbands


# Esempio d'uso e test
if __name__ == "__main__":
    # Dataset di esempio
    s1 = ['minhash', 'is', 'a', 'probabilistic', 'data', 'structure', 'for',
          'estimating', 'the', 'similarity', 'between', 'datasets']
    
    s2 = ['minhash', 'is', 'a', 'probability', 'data', 'structure', 'for',
          'estimating', 'the', 'similarity', 'between', 'documents']
    
    s3 = ['cats', 'are', 'tall', 'and', 'have', 'been',
          'known', 'to', 'sing', 'quite', 'loudly']
    
    print("=== Test MinHash ===")
    
    # Crea hash per ogni lista di parole
    m1 = MinHash()
    m2 = MinHash()
    m3 = MinHash()
    
    # Aggiorna ogni hash
    for word in s1:
        m1.update(word)
    
    for word in s2:
        m2.update(word)
    
    for word in s3:
        m3.update(word)
    
    # Calcola similarità di Jaccard
    print(f"Similarità Jaccard tra s1 e s2: {m1.jaccard(m2):.4f}")
    print(f"Similarità Jaccard tra s1 e s3: {m1.jaccard(m3):.4f}")
    print(f"Similarità Jaccard tra s2 e s3: {m2.jaccard(m3):.4f}")
    
    print("\n=== Test LSH Index ===")
    
    # Crea indice LSH
    index = LSHIndex()
    
    # Inserisce documenti nell'indice
    index.insert('documento1', m1)
    index.insert('documento2', m2)
    index.insert('documento3', m3)
    
    # Query per documenti simili
    matches = index.query(m1)
    print(f"Documenti simili a documento1: {matches}")
    
    matches = index.query(m2)
    print(f"Documenti simili a documento2: {matches}")
    
    matches = index.query(m3)
    print(f"Documenti simili a documento3: {matches}")
    
    print("\n=== Calcolo similarità effettiva di Jaccard ===")
    
    def jaccard_similarity(set1: Set[str], set2: Set[str]) -> float:
        """Calcola la vera similarità di Jaccard tra due insiemi."""
        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))
        return intersection / union if union > 0 else 0.0
    
    # Confronta con la vera similarità di Jaccard
    set1, set2, set3 = set(s1), set(s2), set(s3)
    
    print(f"Vera similarità Jaccard s1-s2: {jaccard_similarity(set1, set2):.4f}")
    print(f"Stima MinHash s1-s2: {m1.jaccard(m2):.4f}")
    print()
    print(f"Vera similarità Jaccard s1-s3: {jaccard_similarity(set1, set3):.4f}")
    print(f"Stima MinHash s1-s3: {m1.jaccard(m3):.4f}")
    print()
    print(f"Vera similarità Jaccard s2-s3: {jaccard_similarity(set2, set3):.4f}")
    print(f"Stima MinHash s2-s3: {m2.jaccard(m3):.4f}")