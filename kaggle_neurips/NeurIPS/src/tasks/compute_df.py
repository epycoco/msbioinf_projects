from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.Chem.rdMolDescriptors import CalcTPSA, CalcNumRotatableBonds
from rdkit.Chem.Descriptors import MolWt, MolLogP
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator, GetAtomPairGenerator, GetTopologicalTorsionGenerator
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from sklearn.preprocessing import StandardScaler
import networkx as nx

from typing import List, Tuple, Optional
from src.logging import lprint, LoggingLevels as ll
from src.config.config import *


def smiles_to_fp (mol_list: List, radius: int = FP_RADIUS, n_bits: int = FP_BITS) -> pd.DataFrame:
    lprint(ll.INFO, "Computing molecular fingerprints...")
    generator = GetMorganGenerator(radius=radius, fpSize=n_bits)
    fingerprints = []

    for i, smiles in enumerate(mol_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            morgan_fp = generator.GetFingerprint(mol)
            maccs_fp = MACCSkeys.GenMACCSKeys(mol)
            combined_fp = np.concatenate([np.array(morgan_fp), np.array(maccs_fp)])
            fingerprints.append(combined_fp)
        else:
            fingerprints.append(np.zeros(n_bits + 167))

    fp_df = pd.DataFrame(fingerprints, columns=[f"fp_{i+1}" for i in range (FP_BITS+167)])
    return fp_df


def smiles_to_desc (mol_list: List, selected_descriptors: List[str]):
    descriptor_functions = {name: func for name, func in Descriptors.descList if name in selected_descriptors}
    descriptors = []
    descriptor_values = {}
    for mol in mol_list:
        for name, func in descriptor_functions.items():
            try:
                descriptor_values[name] = func(mol)
            except:
                descriptor_values[name] = None

        descriptor_values['MolWt'] = MolWt(mol)
        descriptor_values['LogP'] = MolLogP(mol)
        descriptor_values['TPSA'] = CalcTPSA(mol)
        descriptor_values['RotatableBonds'] = CalcNumRotatableBonds(mol)
        descriptor_values['NumAtoms'] = mol.GetNumAtoms()

        try:
            adj = rdmolops.GetAdjacencyMatrix(mol)
            G = nx.from_numpy_array(adj)
            if nx.is_connected(G):
                descriptor_values['graph_diameter'] = nx.diameter(G)
                descriptor_values['avg_shortest_path'] = nx.average_shortest_path_length(G)
            else:
                descriptor_values['graph_diameter'] = 0
                descriptor_values['avg_shortest_path'] = 0
            descriptor_values['num_cycles'] = len(list(nx.cycle_basis(G)))
        except:
            descriptor_values['graph_diameter'] = None
            descriptor_values['avg_shortest_path'] = None
            descriptor_values['num_cycles'] = None

        descriptors.append(descriptor_values)

    # Convert descriptors to DataFrame
    descriptors = [
        {k: v for k, v in d.items() if k != 'SMILES'}
        for d in descriptors if d is not None
    ]

    desc_df = pd.DataFrame(descriptors, columns=selected_descriptors)
    return desc_df



def smiles_to_combined_fingerprints_with_descriptors(smiles_list, selected_descriptors, scaler=None, radius=FP_RADIUS, n_bits=FP_BITS):
    lprint(ll.INFO, "Computing molecular descriptors and fingerprints...")
    generator = GetMorganGenerator(radius=radius, fpSize=n_bits)
    atom_pair_gen = GetAtomPairGenerator(fpSize=n_bits)
    torsion_gen = GetTopologicalTorsionGenerator(fpSize=n_bits)
    descriptor_functions = {name: func for name, func in Descriptors.descList if name in selected_descriptors}
    fingerprints = []
    descriptors = []
    valid_smiles = []
    invalid_indices = []

    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            morgan_fp = generator.GetFingerprint(mol)
            maccs_fp = MACCSkeys.GenMACCSKeys(mol)
            combined_fp = np.concatenate([np.array(morgan_fp), np.array(maccs_fp)])
            fingerprints.append(combined_fp)

            descriptor_values = {}
            for name, func in descriptor_functions.items():
                try:
                    descriptor_values[name] = func(mol)
                except:
                    descriptor_values[name] = None

            descriptor_values['MolWt'] = MolWt(mol)
            descriptor_values['LogP'] = MolLogP(mol)
            descriptor_values['TPSA'] = CalcTPSA(mol)
            descriptor_values['RotatableBonds'] = CalcNumRotatableBonds(mol)
            descriptor_values['NumAtoms'] = mol.GetNumAtoms()

            try:
                adj = rdmolops.GetAdjacencyMatrix(mol)
                G = nx.from_numpy_array(adj)
                if nx.is_connected(G):
                    descriptor_values['graph_diameter'] = nx.diameter(G)
                    descriptor_values['avg_shortest_path'] = nx.average_shortest_path_length(G)
                else:
                    descriptor_values['graph_diameter'] = 0
                    descriptor_values['avg_shortest_path'] = 0
                descriptor_values['num_cycles'] = len(list(nx.cycle_basis(G)))
            except:
                descriptor_values['graph_diameter'] = None
                descriptor_values['avg_shortest_path'] = None
                descriptor_values['num_cycles'] = None

            descriptors.append(descriptor_values)
            valid_smiles.append(smiles)
        else:
            fingerprints.append(np.zeros(n_bits + 167))
            descriptors.append(None)
            valid_smiles.append(None)
            invalid_indices.append(i)

    # Convert descriptors to DataFrame
    descriptors = [
        {k: v for k, v in d.items() if k != 'SMILES'}
        for d in descriptors if d is not None
    ]
    descriptors_df = pd.DataFrame(descriptors)
    
    # Filter selected descriptors
    descriptors_df = descriptors_df.filter(selected_descriptors)
    
    # Standardize the descriptors
    lprint(ll.INFO, f"Standardizing descriptors...")
    if not scaler:
        scaler = StandardScaler()

    desc_values = scaler.fit_transform(descriptors_df)
    descriptors_df = pd.DataFrame(desc_values, columns=descriptors_df.columns)
    lprint(ll.INFO, f"Descriptors standardized: mean={np.mean(desc_values):.2e}, std={np.std(desc_values):.2e}")

    lprint(ll.INFO, f"Computed descriptors shape: {descriptors_df.shape}")
    return descriptors_df, np.array(fingerprints), valid_smiles, invalid_indices, scaler


def compute(df:pd.DataFrame, selected_descriptors: List, scaler = None) -> pd.DataFrame:
    """
    Genera fingerprint molecolari e descrittori per i dataset di train e test separatamente.
    
    Parametri:
        train_df (pd.DataFrame): DataFrame di training con colonna 'SMILES' e altre colonne.
        test_df (pd.DataFrame): DataFrame di test con colonna 'SMILES' e altre colonne.
        selected_descriptors (list): Lista di descrittori molecolari da calcolare.
    
    Ritorna:
        tuple: (train_df, test_df) DataFrame aggiornati con fingerprint e descrittori.
    """
    
    
    lprint(ll.INFO, "Generating molecular fingerprints and descriptors...")
    
    """invalid_ids = []
    smiles_list = df["SMILES"].to_list()
    mol_list = []
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol_list.append(mol)
        else:
            invalid_ids.append(i)
    
    df = df.drop(invalid_ids)

    fp_df = smiles_to_fp(mol_list)
    desc_fp = smiles_to_desc(mol_list, selected_descriptors)

    mol_list = None"""


    descriptors, fingerprints, v_smiles, invalid_ids, scaler = smiles_to_combined_fingerprints_with_descriptors(df["SMILES"].to_list(), selected_descriptors, scaler=scaler)
    descriptors_df = pd.DataFrame(descriptors)
    descriptors_df.filter(selected_descriptors)
    fp_df = pd.DataFrame(fingerprints, columns=[f'FP_{i+1}' for i in range(fingerprints.shape[1])])
    final_df = pd.concat([descriptors_df, fp_df], axis=1)

    #final_df = pd.concat([df, desc_fp, fp_df], axis=1)
    lprint(ll.REPORT, f"Final df shape after feature generation: {final_df.shape}")

    return final_df.sort_index(axis=1), v_smiles, scaler


