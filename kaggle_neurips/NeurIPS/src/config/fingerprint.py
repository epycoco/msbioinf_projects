import os
import pandas as pd
import numpy as np
from typing import List, Tuple, Optional
from src.logging import lprint, LoggingLevels as ll
from src.config.config import *

from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

def compute_morgan_maccs_fp (ismi_df: pd.DataFrame, radius: int = FP_RADIUS, n_bits: int = FP_BITS, phase: str = '') -> pd.DataFrame:
    try:
        idx, smiles = ismi_df['id'].astype(str), ismi_df['SMILES']
    except:
        idx, smiles = ismi_df.index, ismi_df['SMILES']
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]

    fp_dir = f'./data/{phase}_fp.csv'
    if not os.path.exists(fp_dir):
        morgan_generator = GetMorganGenerator(radius=radius, fpSize=n_bits)
        fingerprints = []
        valid_id = []
        for i, mol in enumerate(mols):
            try:
                morgan_fp = morgan_generator.GetFingerprint(mol)
                try:
                    maccs_fp = MACCSkeys.GenMACCSKeys(mol)
                    combined_fp = np.concatenate([np.array(morgan_fp), np.array(maccs_fp)])
                    valid_id.append(idx.to_list()[i])
                except Exception as e:
                    lprint(ll.WARN, f"Impossible to compute MACCSKeys for '{smiles.to_list()[i]}': {e}")
                    combined_fp = np.zeros(shape=(n_bits+167))

            except Exception as e:
                lprint(ll.WARN, f"Impossible to compute Morgan Fingerprint for '{smiles[i]}': {e}")
                combined_fp = np.zeros(shape=(n_bits+167))
            
            combined_fp = np.concatenate([np.array(morgan_fp), np.array(maccs_fp)])
            fingerprints.append(combined_fp)

        fp_df = pd.DataFrame(fingerprints, columns=[f"mmfp_{i+1}" for i in range (FP_BITS+167)], dtype=int)
        ifp_df = pd.concat([idx, fp_df], axis=1)
        ifp_df = ifp_df[ifp_df['id'].isin(valid_id)]
        ifp_df.to_csv(fp_dir, index=False)

    else:
        ifp_df = pd.read_csv(fp_dir)
        lprint(ll.INFO, f'A Fingerprint file was found: UPLOADED')
    
    ifp_df['id'] = ifp_df['id'].astype(str)
    return ifp_df