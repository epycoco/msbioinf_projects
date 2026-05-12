import os
from typing import Union, List
import warnings
from src.config.config import *
from src.logging import lprint, LoggingLevels as ll

import pandas as pd
import numpy as np

import torch
from transformers import AutoTokenizer, AutoModel



hf_token = "" #TODO: use your token here
MODEL_NAME = "DeepChem/ChemBERTa-77M-MLM"
TOKENIZER = AutoTokenizer.from_pretrained(MODEL_NAME, token=hf_token)

def compute_chemberta(df: pd.DataFrame, 
                        model_name: str = MODEL_NAME,
                        max_length: int = 320,
                        batch_size: int = 32,
                        device: str = None,
                        phase: str = '') -> pd.DataFrame:
    """
    Genera embedding ChemBERTa da una Series di SMILES.
    
    Parameters:
    -----------
    smiles : pd.Series
        Serie pandas contenente le stringhe SMILES
    model_name : str, default="DeepChem/ChemBERTa-77M-MTR"
        Nome del modello ChemBERTa da utilizzare
    max_length : int, default=314
        Lunghezza massima delle sequenze tokenizzate
    batch_size : int, default=32
        Dimensione del batch per l'elaborazione
    device : str, optional
        Device da utilizzare ('cpu', 'cuda', 'mps'). Se None, viene scelto automaticamente
    
    Returns:
    --------
    pd.DataFrame
        DataFrame con embedding come colonne (embedding_0, embedding_1, ..., embedding_n)
        L'indice corrisponde all'indice originale della Series di input
    """
    idx, smiles = df['id'].astype(str), df['SMILES'].to_list()
    berta_dir = f'./data/{phase}_chemberta.csv'

    if not os.path.exists(berta_dir):
        if device is None:
            if torch.cuda.is_available():
                device = 'cuda'
            elif torch.backends.mps.is_available():
                device = 'mps'
            else:
                device = 'cpu'
        
        print(f"Utilizzando device: {device}")
        
        # Carica tokenizer e modello
        lprint(ll.INFO, f"Upload of ChemBERTa model...")
        model = AutoModel.from_pretrained(model_name, token=hf_token)
        model = model.to(device)
        model.eval()
        
        embeddings = []
        
        # Elabora in batch
        with torch.no_grad():
            for i in range(0, len(smiles), batch_size):
                batch_smiles = smiles[i:i+batch_size]
                
                try:
                    # Tokenizzazione
                    inputs = TOKENIZER(
                        batch_smiles, 
                        padding=True, 
                        truncation=True, 
                        max_length=max_length,
                        return_tensors="pt"
                    )
                    
                    # Sposta i tensori sul device
                    inputs = {k: v.to(device) for k, v in inputs.items()}
                    
                    # Forward pass
                    outputs = model(**inputs)
                    
                    # Estrai embedding dal token [CLS] (primo token)
                    batch_embeddings = outputs.last_hidden_state[:, 0, :].cpu().numpy()
                    embeddings.extend(batch_embeddings)
                    
                    # Progress update
                    if (i // batch_size + 1) % 10 == 0 or i + batch_size >= len(smiles):
                        print(f"Processati {min(i + batch_size, len(smiles))}/{len(smiles)} SMILES")
                        
                except Exception as e:
                    warnings.warn(f"Errore nel processare batch {i//batch_size + 1}: {str(e)}")
                    # Aggiungi embedding nulli per questo batch
                    batch_size_actual = len(batch_smiles)
                    embedding_dim = 768 if 'embedding_dim' not in locals() else embedding_dim
                    null_embeddings = np.zeros((batch_size_actual, embedding_dim))
                    embeddings.extend(null_embeddings)
        
        # Converti in array numpy
        embeddings_array = np.array(embeddings)
        
        # Crea DataFrame
        embedding_columns = [f'chemberta_{i}' for i in range(embeddings_array.shape[1])]
        berta_df = pd.DataFrame(
            embeddings_array, 
            columns=embedding_columns,
            index=idx.index
        )
        
        iberta_df = pd.concat([idx, berta_df], axis=1)
        
        iberta_df.to_csv(berta_dir, index=False)
        lprint(ll.SUCCESS, f"ChemBERTa Embedding Completed.")

    else:
        iberta_df = pd.read_csv(berta_dir, header=0)
        lprint(ll.INFO, f'A ChemBERTa file was found: UPLOADED')

    iberta_df['id'] = iberta_df['id'].astype(str)
    return iberta_df


