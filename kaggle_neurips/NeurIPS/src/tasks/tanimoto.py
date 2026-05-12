import os, sys
from typing import List, Tuple, Optional

import pandas as pd
import numpy as np

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


def compute_tanimoto(bits_a, bits_b):
    a = np.asarray(bits_a, dtype=bool)
    b = np.asarray(bits_b, dtype=bool)
    inter = np.logical_and(a, b).sum()
    union = np.logical_or(a, b).sum()
    return inter / union if union else 0.0


def tanimoto (fp_df1: pd.DataFrame, fp_df2: pd.DataFrame):
    df1_rows = fp_df1.shape[0]
    df2_rows = fp_df2.shape[0]
    top_sim = []

    for l2 in range(df2_rows):
        top = 0
        for l1 in range(df1_rows):
            similarity = compute_tanimoto(fp_df2.iloc[l2, :], fp_df1.iloc[l1, :])
            if similarity == 1.0:
                top = similarity if similarity > top else top
        top_sim.append(top)
    
    return top_sim



