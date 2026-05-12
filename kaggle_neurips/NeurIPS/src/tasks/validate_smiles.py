from typing import List, Tuple, Optional

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from src.logging import lprint, LoggingLevels as ll


def change_atom_holders(smiles:pd.Series, atom_holders = '*', substitution_group = 'C') -> pd.Series:
    """
        Take in input a series that contain a the SMILES
        and substituites the atom holders ('*' or another char that you can set)
        with an functional group or another molecule that you desire to insert.
    """
    return smiles.str.replace(atom_holders, substitution_group)


def get_canonical_smiles(smiles: str):
        """Convert SMILES to canonical form for consistency"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Chem.MolToSmiles(mol, canonical=True)
            else:
                return None
        except:
            return None


# Robust Data Loading with Complete R-Group Filtering
def validate_smiles_2d(smiles: str):
    """Completely clean and validate SMILES, removing all problematic patterns"""
    if not isinstance(smiles, str) or len(smiles) == 0:
        return None
    
    # List of all problematic patterns we've seen
    bad_patterns = [
        '[R]', '[R1]', '[R2]', '[R3]', '[R4]', '[R5]', 
        "[R']", '[R"]', 'R1', 'R2', 'R3', 'R4', 'R5',
        # Additional patterns that cause issues
        '([R])', '([R1])', '([R2])', 
    ]
    
    # Check for any bad patterns
    for pattern in bad_patterns:
        if pattern in smiles:
            return None
    
    # Additional check: if it contains ] followed by [ without valid atoms, likely polymer notation
    if '][' in smiles and any(x in smiles for x in ['[R', 'R]']):
        return None
    
    return get_canonical_smiles(smiles)
    

def validate_smiles_3d(smiles: str) -> str | None:
    """
    Converte un SMILES in una molecola RDKit con coordinate 3D ottimizzate.
    Ritorna None se la generazione fallisce.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # aggiunge idrogeni
        # genera conformero 3D
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            return None
        # ottimizza con MMFF (se disponibile, altrimenti UFF)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            AllChem.UFFOptimizeMolecule(mol)
        return smiles
    except:
        return None


def main(df1: pd.DataFrame, df2: Optional[pd.DataFrame] = None, not3d: bool = False):
    """
        Main function to validate SMILES:
            if only 'df1' is submitted -> it simply compute validation on SMILES column 
            if 'df1' and 'df2' are submitted -> the validation of SMILES 3D structures is based on df2 SMILES compatibility with RDKit 3D computations
    """
    if df2 is not None:
        df2['SMILES'] = change_atom_holders(df2['SMILES'])
        df2['SMILES'] = df2['SMILES'].apply(validate_smiles_2d)
        dup = df2["SMILES"].duplicated().to_list().count(True)
        if dup > 0:
            df2 = df2.drop_duplicates(subset='SMILES', keep='first') # Remove Duplicated SMILES
            lprint(ll.INFO, f"Removed {dup} duplicated SMILES from training data")
            
        for smi in df2['SMILES']:
            if smi is None:
                lprint(ll.WARN, f"It is impossibile to compute 3D descriptors because for one SMILES in test file is not possible to compute 3D structures")
                not3d = True
                break

    df1['SMILES'] = change_atom_holders(df1['SMILES'])
    df1['SMILES'] = df1['SMILES'].apply(validate_smiles_2d)
    
    dup = df1["SMILES"].duplicated().to_list().count(True)
    if dup > 0:
        df1 = df1.drop_duplicates(subset='SMILES', keep='first') # Remove Duplicated SMILES
        lprint(ll.INFO, f"Removed {dup} duplicated SMILES from training data")

    if not3d is False:
        if df2 is not None: 
            df2['SMILES'] = df2['SMILES'].apply(validate_smiles_3d)
            if any(df2['SMILES'].isnull()):
                lprint(ll.WARN, f"It is impossibile to compute 3D descriptors because of one SMILES in test file.")
                not3d = True
                return df1, df2, not3d
            else:
                lprint(ll.INFO, f"Validating 3D structures for df1...")
                df1['SMILES'] = df1['SMILES'].apply(validate_smiles_3d)
                invalid_train = df1['SMILES'].isnull().sum()
                df1 = df1[df1['SMILES'].notnull()].reset_index(drop=True)
                lprint(ll.INFO, f"Removing {invalid_train} invalid SMILES from training data")
        else:
            lprint(ll.INFO, f"Validating 3D structures for df1...")
            df1['SMILES'] = df1['SMILES'].apply(validate_smiles_3d)
            invalid_train = df1['SMILES'].isnull().sum()
            df1 = df1[df1['SMILES'].notnull()].reset_index(drop=True)
            lprint(ll.INFO, f"Removing {invalid_train} invalid SMILES from training data")

    lprint(ll.REPORT, f"Final df1 samples: {df1.shape[0]}")
    if df2 is not None:
        lprint(ll.REPORT, f"Final df2 samples: {df2.shape[0]}")
        return df1, df2, not3d

    return df1, None, not3d