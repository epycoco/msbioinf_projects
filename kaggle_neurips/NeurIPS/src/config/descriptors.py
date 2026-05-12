
import os
import pandas as pd
import numpy as np
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from typing import List, Tuple, Optional
from src.logging import lprint, LoggingLevels as ll
from src.config.config import *

# Descriptors
from mordred import Calculator, descriptors

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem import Descriptors, Descriptors3D
from rdkit.Chem.rdMolDescriptors import CalcTPSA, CalcNumRotatableBonds

# Graphs
import networkx as nx

# Descriptors Selection 
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# verbose
from tqdm import tqdm

base_rdkit_desc = {'graph_diameter','num_cycles','avg_shortest_path','MolWt', 'LogP', 'TPSA', 'RotatableBonds', 'NumAtoms'}
filters = {
    'Tg': list(set([
        "MolWt","HeavyAtomCount","NumAromaticRings","NumAliphaticRings","NumRotatableBonds","NumHAcceptors",
        "NumHDonors","FractionCSP3",
        # --- Topologici / connettività ---
        "Chi0", "Chi1", "Chi0n", "Chi1n", "Chi2n", "Chi3n", "Chi4n",
        "Kappa1", "Kappa2", "Kappa3", 
        "BalabanJ","BertzCT","HallKierAlpha","Ipc",
        # --- Fisico-chimici / polarità ---
        "MolLogP","MolMR","TPSA","NumValenceElectrons","MaxPartialCharge","MinPartialCharge","MaxAbsPartialCharge",

        # --- Elettrochimici / elettronici ---
        "FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3",  # vari livelli di Morgan density
        "EState_VSA1", "EState_VSA2", "EState_VSA3",  # esempi di descrittori EState VSA
        "VSA_EState1", "VSA_EState2", "VSA_EState3",  # complementari

        # --- Fingerprint (bitvector) ---
        "MorganFingerprint",  # da generare separatamente (non è un singolo descrittore)
        "TopologicalTorsionFingerprint",
        "AtomPairsFingerprint"
    ]).union(base_rdkit_desc)),

    'FFV': list(set([
        'AvgIpc','BalabanJ','BertzCT','Chi0','Chi0n','Chi0v','Chi1','Chi1n','Chi1v',
        'Chi2n','Chi2v','Chi3n','Chi3v','Chi4n','EState_VSA10','EState_VSA5',
        'EState_VSA7','EState_VSA8','EState_VSA9','ExactMolWt','FpDensityMorgan1',
        'FpDensityMorgan2','FpDensityMorgan3','FractionCSP3','HallKierAlpha',
        'HeavyAtomMolWt','Kappa1','Kappa2','Kappa3','MaxAbsEStateIndex',
        'MaxEStateIndex','MinEStateIndex','MolLogP','MolMR','MolWt','NHOHCount',
        'NOCount','NumAromaticHeterocycles','NumHAcceptors','NumHDonors',
        'NumHeterocycles','NumRotatableBonds','PEOE_VSA14','RingCount','SMR_VSA1',
        'SMR_VSA10','SMR_VSA3','SMR_VSA5','SMR_VSA6','SMR_VSA7','SMR_VSA9','SPS',
        'SlogP_VSA1','SlogP_VSA10','SlogP_VSA11','SlogP_VSA12','SlogP_VSA2',
        'SlogP_VSA3','SlogP_VSA4','SlogP_VSA5','SlogP_VSA6','SlogP_VSA7',
        'SlogP_VSA8','TPSA','VSA_EState1','VSA_EState10','VSA_EState2',
        'VSA_EState3','VSA_EState4','VSA_EState5','VSA_EState6','VSA_EState7',
        'VSA_EState8','VSA_EState9','fr_Ar_N','fr_C_O','fr_NH0','fr_NH1',
        'fr_aniline','fr_ether','fr_halogen','fr_thiophene'
    ]).union(base_rdkit_desc)),

    'Tc': list(set([
        'BalabanJ','BertzCT','Chi0','EState_VSA5','ExactMolWt','FpDensityMorgan1',
        'FpDensityMorgan2','FpDensityMorgan3','HeavyAtomMolWt','MinEStateIndex',
        'MolWt','NumAtomStereoCenters','NumRotatableBonds','NumValenceElectrons',
        'SMR_VSA10','SMR_VSA7','SPS','SlogP_VSA6','SlogP_VSA8','VSA_EState1',
        'VSA_EState7','fr_NH1','fr_ester','fr_halogen'
    ]).union(base_rdkit_desc)),

    'Density': list(set([
        'BalabanJ','Chi3n','Chi3v','Chi4n','EState_VSA1','ExactMolWt',
        'FractionCSP3','HallKierAlpha','Kappa2','MinEStateIndex','MolMR','MolWt',
        'NumAliphaticCarbocycles','NumHAcceptors','NumHeteroatoms',
        'NumRotatableBonds','SMR_VSA10','SMR_VSA5','SlogP_VSA12','SlogP_VSA5',
        'TPSA','VSA_EState10','VSA_EState7','VSA_EState8'
    ]).union(base_rdkit_desc)),

    'Rg': list(set([
        'AvgIpc','Chi0n','Chi1v','Chi2n','Chi3v','ExactMolWt','FpDensityMorgan1',
        'FpDensityMorgan2','FpDensityMorgan3','HallKierAlpha','HeavyAtomMolWt',
        'Kappa3','MaxAbsEStateIndex','MolWt','NOCount','NumRotatableBonds',
        'NumUnspecifiedAtomStereoCenters','NumValenceElectrons','PEOE_VSA14',
        'PEOE_VSA6','SMR_VSA1','SMR_VSA5','SPS','SlogP_VSA1','SlogP_VSA2',
        'SlogP_VSA7','SlogP_VSA8','VSA_EState1','VSA_EState8','fr_alkyl_halide',
        'fr_halogen'
    ]).union(base_rdkit_desc))
}



def is_numeric(series: pd.Series) -> bool:
    """
    Function that verify if the Series contains only numeric values
    """
    s = series.dropna().astype(str)
    is_all_numeric = pd.to_numeric(s, errors='coerce').notna().all()
    return is_all_numeric 


def smiles_to_3d(smiles: str):
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
        return mol
    except:
        return None
    

def base_filter_desc(idesc_df: pd.DataFrame) -> pd.DataFrame:
    idx = idesc_df['id']
    n_mols = len(idx)

    desc_df = idesc_df.drop(columns='id')
    desc_list = desc_df.columns.to_list()
    try:
        valid_desc = [des for des in desc_list if is_numeric(desc_df[des])]
        desc_df = desc_df[valid_desc]
        desc_df = desc_df.astype(np.float64)
    except Exception as e:
        lprint(ll.ERROR, f'An error occurred during the validation of descriptors values: {e}')
        return None

    idesc_df = pd.concat([idx, desc_df], axis=1)
    idt_rep_list = idesc_df.apply(lambda col: col.value_counts(normalize=True, dropna=True).max())

    desc_to_drop = []
    for i, desc in enumerate(valid_desc, 1):
        desc_values = desc_df[desc]
        desc_var = desc_values.var()
        itd_nan = idesc_df.loc[desc_values.isna(), 'id'].unique().tolist() # Id To Drop because of NaN values
        itd_rep = idt_rep_list[i]

        # filter on variance
        if desc_var <= VAR_THRESH:
            desc_to_drop.append(desc)

        # if the number of NaN > of the 2% of total molecules, select the descriptors in order to drop it
        elif len(itd_nan) > n_mols*NAN_THRESH:
            desc_to_drop.append(desc)
            # I called the variable "Id To Drop" because before I chose to drop them, but then I decided to do INFERENCE on Nan values for descriptors
            #   so I do NOT eliminate any rows

        # If the number of repeted values is greater than 90% of total molecules, select the descriptors in order to drop it
        elif itd_rep >= REP_THRESH:
            desc_to_drop.append(desc)
    
    # Drop not valid desc and save the index of not valid molecules
    if desc_to_drop is not None:
        idesc_df = idesc_df.drop(columns=desc_to_drop)
    
    return idesc_df


def compute_all_rdkit_desc(ismi_df: pd.DataFrame,
                           not3d:bool, 
                           desc_list: Optional[List[str]] = None, 
                           desc_dont: Optional[List[str]] = None,
                           phase: str = '') -> pd.DataFrame:
    """
    Compute 2D and 3D descriptors of RDKit from a pandas DataFrame containing SMILES.
    It requires that 3D strutures are generable from SMILES (using rdkit.Chem.AllChem.ETKDG + rdkit.Chem.AllChem.MMFF/UFF).
    
    Args:
        -   ismi_df (pd.DataFrame) : DataFrame with following columns ['id', 'SMILES'].
        -   rdkit_desc_list (Optional[List[str]]) : Optional list with the RDKit Descriptors name selected from all the descriptors of RDKit package
    
    Returns:
        -   desc_df (pd.DataFrame): DataFrame with 'id' and valid RDKit 2D/3D descriptors.
    """

    idx = ismi_df['id']
    smiles = ismi_df['SMILES']
    
    # 2D descriptors
    rdk2_dir = f'./data/{phase}_rdkit2.csv'
    if not os.path.exists(rdk2_dir):
        desc_values = {}
        results = []

        mols = [Chem.MolFromSmiles(smi) for smi in smiles]
        desc2D_list = Descriptors.descList
        if desc_dont:
            desc2D_list = [d for d in desc2D_list if d not in desc_dont]
        if desc_list:
            desc2D_list = [desc for desc in desc2D_list if desc[0] in desc_list]
        desc2D_names: List[str] = [d[0] for d in desc2D_list]
        desc2D_funcs = {name: func for name, func in desc2D_list}
        lprint(ll.INFO, f'Computing 2D RDKit {len(desc2D_names)} descriptors')

        for i, mol in enumerate(tqdm(mols, desc="2D Descriptors", unit="mol")):   
            desc_values = {}
            desc_values['id'] = idx[i]    
            for name, func in desc2D_funcs.items():
                try:
                    desc_values[name] = func(mol)
                except Exception as e:
                    desc_values[name] = None
                    lprint(ll.WARN, f"An error occur during '{name}' computation for the molecule {i}: {e}")
            results.append(desc_values)
        idesc2_df = pd.DataFrame(results)
        desc2D_names = idesc2_df.drop(columns='id').columns.to_list()
        lprint(ll.INFO, f'2D RDKit descriptors: COMPUTED.')
        idesc2_df.to_csv(rdk2_dir, index=False)

    else:
        idesc2_df = pd.read_csv(rdk2_dir, header=0)
        idesc2_df['id'] = idesc2_df['id'].astype(str)
        
        desc2D_names = idesc2_df.drop(columns='id').columns.to_list()
        lprint(ll.INFO, f'A 2D RDKit descriptors file was found: UPLOADED')

    # 3D descriptors
    if not3d is False:
        rdk3_dir = f'./data/{phase}_rdkit3.csv'
        if not os.path.exists(rdk3_dir):
            desc_values = {}
            results = []
            mols = [smiles_to_3d(smi) for smi in smiles]
            desc3D_list = Descriptors3D.descList
            if desc_dont:
                desc3D_list = [d for d in desc3D_list if d not in desc_dont]
            if desc_list:
                desc3D_list = [desc for desc in desc3D_list if desc[0] in desc_list]
            desc3D_names: List[str] = [d[0] for d in desc3D_list]
            desc3D_funcs = {name: func for name, func in desc3D_list}
            lprint(ll.INFO, f'Computing 3D RDKit {len(desc3D_names)} descriptors')
            
            for i, mol in enumerate(tqdm(mols, desc="3D Descriptors", unit="mol")):
                desc_values = {}
                desc_values['id'] = idx[i]
                for name, func in desc3D_funcs.items():
                    try:
                        desc_values[name] = func(mol)
                    except Exception as e:
                        desc_values[name] = None
                        lprint(ll.WARN, f"An error occur during '{name}' computation for the molecule {i}: {e}")
                results.append(desc_values)
            idesc3_df = pd.DataFrame(results)
            desc3D_names = idesc3_df.drop(columns='id').columns.to_list()
            lprint(ll.INFO, f'3D RDKit descriptors: COMPUTED.')
            idesc3_df.to_csv(rdk3_dir, index=False)
        else:
            idesc3_df = pd.read_csv(rdk3_dir, header=0)
            idesc3_df['id'] = idesc3_df['id'].astype(str)
            desc3D_names = idesc3_df.drop(columns='id').columns.to_list()
            lprint(ll.INFO, f'A 3D RDKit descriptors file was found: UPLOADED')

        try:
            idesc_df = pd.merge(idesc2_df, idesc3_df, on='id')
        except Exception as e:
            lprint(ll.ERROR, f"Merge between 2D and 3D dscriptors dataframe FAILED: {e}")
    
    else: 
        idesc_df = idesc2_df
    
    if phase != 'test':
        idesc_df = base_filter_desc(idesc_df)
        lprint(ll.INFO, f'RDKit descriptors dataset FILTERED.')

        if idesc_df.shape[0] == 0:
            lprint(ll.WARN, f"There is no valid molecules after RDKit Descriptors computation. The analysis will continue without RDKit Descriptors.")
            return None

    lprint(ll.REPORT, f"RDKit Descriptors DataFrame shape (include 'id' column): {idesc_df.shape}")
    idesc_df = idesc_df.reset_index()
    return idesc_df


def compute_all_mordred_desc (ismi_df: pd.DataFrame, 
                              not3d: bool, desc_list: Optional[List[str]] = None, 
                              desc_dont: Optional[List[str]] = None,
                              phase: str = '') -> pd.DataFrame:
    """
    Compute 2D and 3D descriptors of Mordred from a pandas DataFrame containing molecules 'id' and 'SMILES'.
    It requires that 3D strutures are generable from SMILES (using rdkit.Chem.AllChem.ETKDG + rdkit.Chem.AllChem.MMFF/UFF).
    
    Args:
        - ismi_df (pd.DataFrame) : DataFrame with following columns: ['id', 'SMILES'].
        - mordred_desc_list (Optional[List[str]]) : Optional list with the Mordred Descriptors name selected from all the descriptors of Mordred package
    
    Returns:
        - desc_df (pd.DataFrame) : DataFrame with 'id' and valid Mordred 2D/3D descriptors
    """
    idx, smiles = ismi_df['id'], ismi_df['SMILES']

    calculator = Calculator(descriptors, ignore_3D=not3d)
    
    if desc_list:
        all_mordred_desc_list = list(calculator.descriptors)
        desc_to_compute = [d for d in all_mordred_desc_list if str(d) in desc_list]
        calculator = Calculator(desc_to_compute, ignore_3D=not3d)
        
    elif desc_dont is not None:
        all_mordred_desc_list = list(calculator.descriptors)
        desc_to_compute = [d for d in all_mordred_desc_list if str(d) not in desc_dont]
        try:
            calculator = Calculator(desc_to_compute, ignore_3D=not3d)
        except Exception as e:
            lprint(ll.ERROR, f'Impossible to compute Mordred Descriptors: {e}')
            raise BrokenPipeError
    
    mordred_dir = f'./data/{phase}_mordred.csv'
    if not os.path.exists(mordred_dir):
        mols = [Chem.MolFromSmiles(smi) for smi in smiles]
        try:
            desc_df = calculator.pandas(mols)
            idesc_df = pd.concat([idx, desc_df], axis=1)
        except Exception as e:
            lprint(ll.WARN, f'Computation of Mordred Descriptors FAILED: {e}')
            return None

        idesc_df.to_csv(mordred_dir, index=False)

    else:
        idesc_df = pd.read_csv(mordred_dir)
        lprint(ll.INFO, f'A Mordred descriptors file was found: UPLOADED')
    
    if phase != 'test':
        idesc_df = base_filter_desc(idesc_df)
        lprint(ll.INFO, f'Mordred descriptors dataset FILTERED.')

        if idesc_df.shape[0] == 0:
            lprint(ll.WARN, f"There is no valid molecules after Mordred Descriptors filtering. The analysis will continue without Mordred Descriptors.")
            return None
            
    lprint(ll.REPORT, f"Mordred Descriptors DataFrame shape: {idesc_df.shape}")
    return idesc_df


def compute_all_desc (ismi_df: pd.DataFrame, 
                      not3d: bool, 
                      desc_list: Optional[List[str]] = None,
                      phase: str = '') -> pd.DataFrame:
    desc_dir = f"./data/{phase}_desc.csv"

    if not os.path.exists(desc_dir):
        ismi_df['id'] = ismi_df['id'].astype(str)
        
        lprint(ll.STEP, f'Computing RDKit Descriptors')
        idesc_rdk_df = compute_all_rdkit_desc(ismi_df= ismi_df, not3d= not3d, desc_list= desc_list, phase= phase)
        valid_rdk_id = idesc_rdk_df['id'].to_list()
        lprint(ll.SUCCESS, f'RDKit Descriptors Computed')

        desc_dont = idesc_rdk_df.drop(columns='id').columns.to_list() + ['ABC', 'ABCGG']
        ismi_df = ismi_df[ismi_df['id'].isin(valid_rdk_id)]
        
        lprint(ll.STEP, f'Computing Mordred Descriptors')
        idesc_mrd_df = compute_all_mordred_desc(ismi_df= ismi_df, not3d= not3d, desc_list= desc_list, desc_dont= desc_dont, phase= phase)
        valid_mrd_id = idesc_mrd_df['id'].to_list()
        lprint(ll.SUCCESS, f'Mordred Descriptors Computed')

        if idesc_rdk_df is not None:
            idesc_rdk_df['id'] = idesc_rdk_df['id'].astype(str)
            lprint (ll.REPORT, f"RDKit Descriptors DataFrame shape: {idesc_rdk_df.shape}")
        if idesc_mrd_df is not None:
            lprint (ll.REPORT, f"Mordred Descriptors DataFrame shape: {idesc_mrd_df.shape}")
            idesc_mrd_df['id'] = idesc_mrd_df['id'].astype(str)

        if (idesc_rdk_df is None) and (idesc_mrd_df is None):
            lprint(ll.WARN, f"There is no one valid molecules found during computation of Mordred/RDKit Descriptors. Descriptors feature would be omitted\n")
            return None
        
        elif phase == 'test':
            idesc_df = pd.merge(idesc_rdk_df, idesc_mrd_df, on="id")

        elif idesc_rdk_df is None:
            idesc_df = idesc_mrd_df[idesc_mrd_df['id'].isin(valid_mrd_id)] # selection of valid molecules

        elif idesc_mrd_df is None:
            idesc_df = idesc_rdk_df[idesc_rdk_df['id'].isin(valid_rdk_id)] # selection of valid molecules

        else:
            idesc_df = pd.merge(idesc_rdk_df, idesc_mrd_df, on="id")

        if phase != 'test':
            idesc_df = base_filter_desc(idesc_df= idesc_df)
        
        idesc_df.to_csv(desc_dir, index=False)
        lprint (ll.REPORT, f"All Descriptors DataFrame shape: {idesc_df.shape}")

    else:
        idesc_df = pd.read_csv(desc_dir, header=0)
        lprint(ll.INFO, f'A Descriptors file was found: UPLOADED')

    return idesc_df


def decrease_corr_min_yn(filter_corr: float):
        if filter_corr > 0.675:
            filter_corr = 0.675
        else:
            filter_corr -= 0.025

        yn = input(f"Do you want to decrease the minimun correlation threshold value at '{filter_corr}'? [Y/n]")
        if yn.strip()=='Y' or yn.strip()=='y' or yn.strip()==None:
            return filter_corr
        elif yn.strip()=='N' or yn.strip()=='n':
            return None
        else:
            lprint(ll.WARN, f"Digit only [Y] or [N]")
            decrease_corr_min_yn(filter_corr)



def filter_desc_pca_corr (desc_df: pd.DataFrame, 
                          target: str,
                          filter_by_loding: bool = True, 
                          filter_variance: float = 0.9,
                          filter_quantile: float = 0.70,
                          filter_correlation_max:float = 0.89,
                          filter_correlation_min:float = 0.70) -> pd.DataFrame:
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(desc_df)
    pca = PCA()
    X_pca = pca.fit_transform(X_scaled)
    cum_var = pca.explained_variance_ratio_.cumsum()

    loadings = pd.DataFrame(
        pca.components_.T,
        columns=[f'PC{i+1}' for i in range(pca.n_components_)],
        index = desc_df.columns
    )

    # Select PCs with the method of cumulative variance (filter_variance)    
    n_pc = len(cum_var)
    n_pc = next(i for i, v in enumerate(cum_var) if v > filter_variance)+1

    sel_PCs = [f"PC{i}" for i in range (n_pc)]
    scores_df = pd.DataFrame(X_pca[:, :n_pc], columns=sel_PCs)
    
    selected_descriptor = list()
    for i in range (n_pc):
        loading_ord = loadings.iloc[:, i].abs().sort_values(ascending=False) # descending order of descriptors in loading matrix feature
        
        # select just a bunch of descriptors per loadings features using a quantile filter threshold
        quantile_thres = loading_ord.quantile(filter_quantile) # I choose the loading value of the q-th quantile
        j=0
        while loading_ord.iloc[j] > quantile_thres:
            j+=1
        selected_descriptor += loading_ord[:j].index.tolist()
    
    # create a correlation matrix of descriptor selected from loadings 
    selected_descriptor = set(selected_descriptor)
    sel_desc_df = desc_df[list(selected_descriptor)].copy()
    
    corr_matrix = sel_desc_df.corr()

    fig = px.imshow(
        corr_matrix,
        text_auto=False,             # NON mostra i valori nelle celle
        aspect="equal",              # proporzioni corrette
        color_continuous_scale="RdBu_r",  # blu-rosso (positivo-negativo)
        title=f"{target} - Correlation between Descriptors"
        )
    fig.write_image(file=f"C:/Users/andre/Desktop/machine_learning_prjct/FINAL/ppt/{target}_desc_corr_desc", format='png')
    fig.show()

    # find the best value for each descriptors selected in all the loadings chose
    diz_best_descriptor_load = {}  
    for i in range(n_pc):
        series = loadings.iloc[:, i].abs().fillna(0.0)
        for feat, val in series.items():
            if val > diz_best_descriptor_load.get(feat, float("-inf")):
                diz_best_descriptor_load[feat] = (val)

    # descending order of descriptors selected
    best_descriptor_load = pd.Series(diz_best_descriptor_load).sort_values(ascending=False)
    print (best_descriptor_load)
    mean_loading = best_descriptor_load.mean()

    # verify if there is correlation between the descriptors selected
    descriptor_to_remove = set()
    for i in range (len(corr_matrix.columns)):
        for j in range (i+1, len(corr_matrix.columns)):

            idx_d1 = best_descriptor_load.index.get_loc(corr_matrix.columns[i])
            idx_d2 = best_descriptor_load.index.get_loc(corr_matrix.columns[j])

            # if correlation is over a threshold (filter_correlation_max), the descriptor to remove is that with the worst value in the descriptors list obtained
            if corr_matrix.iloc[i,j] > filter_correlation_max:
                descriptor_to_remove.add(best_descriptor_load.index[max(idx_d1, idx_d2)])

            # if correlation is around the min and max threshold, and at least one of two values is under the mean of desc_values
            # the descriptor to remove is that with the worst value in the descriptors list obtained
            if filter_by_loding and (filter_correlation_min < corr_matrix.iloc[i,j] < filter_correlation_max):
                if min(best_descriptor_load[corr_matrix.columns[i]], best_descriptor_load[corr_matrix.columns[j]]) < mean_loading:
                    descriptor_to_remove.add(best_descriptor_load.index[max(idx_d1, idx_d2)])

    descriptor_to_remove = list(descriptor_to_remove)
    final_descriptors = [col for col in selected_descriptor if col not in descriptor_to_remove]

    return desc_df[final_descriptors]



def filter_desc_corrTarget (desc_df: pd.DataFrame, target: pd.Series, t: str, corr_min: float = 0.60) -> pd.DataFrame | None:
    """
        Base the description selection on the correlation with the targets values
    """
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(desc_df))
    
    df = pd.concat([X_scaled, target], axis=1)
    df.columns = desc_df.columns.to_list() + [t]
    df.to_csv(f"./data/desc_corr_target/{t}_desc_target.csv", index=False)
    corr_matrix = df.corr(method='pearson')

    corr_target = corr_matrix[t].sort_values(ascending=False)
    corr_df = pd.DataFrame([corr_target.values], columns=corr_target.index)
    corr_df.to_csv(f"./data/{t}_desc_corr_with_target.csv", index=False)

    not_corr_target = [[i, abs(cor)] for i, cor in enumerate(corr_matrix[t]) if abs(cor)>0 and abs(cor)<corr_min]
    corr_target_ivv = [[i, abs(cor)] for i, cor in enumerate(corr_matrix[t]) if abs(cor)>=corr_min and abs(cor)<1]
    
    if len(not_corr_target) == corr_matrix.shape[0]-1:
        lprint(ll.REPORT, f'There are not Descriptors correlated.')
        return None
    
    elif corr_target_ivv == []:
        lprint(ll.REPORT, f'There are not Descriptors correlated with min correlation threshold = {corr_min}.')
        return None
    
    else:
        lprint(ll.REPORT, f'{len(corr_target_ivv)} Descriptors selected (min correlation threshold = {corr_min}.')
        desc_all_list = desc_df.columns.to_list()
        desc_list = [desc_all_list[i] for i, cor in corr_target_ivv]
        return desc_df[desc_list]



























def compute_graph_desc (mol_list: List):
    descriptors = []
    descriptor_values = {}
    for mol in mol_list:
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

    graph_df = pd.DataFrame(descriptors, columns=["graph_diameter", "avg_shortest_path", "num_cycles"])
    graph_df = base_filter_desc(graph_df)

    return graph_df
