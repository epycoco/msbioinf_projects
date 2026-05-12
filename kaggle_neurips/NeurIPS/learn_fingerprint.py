import os
from typing import Any, List, Optional, Tuple
import pandas as pd
import numpy as np

from rdkit.Chem import rdFingerprintGenerator

from rdkit import Chem
from rdkit.Chem import Mol, Conformer
from rdkit.Chem import AllChem, rdMolDescriptors, rdMolTransforms, rdMolAlign
from rdkit.ML.Cluster import Butina

from src.tasks import validate_smiles
from src.logging import lprint, LoggingLevels as ll
from src.config.config import *


# rdMolAlign.GetBestRMS(prbMol= mol, refMol= mol, prbId= cands[i], refId= kid, numThreads=CPU): between conformers of the same molecule
# rdMolAlign.GetBestRMS: between conformers of different molecules

E3_BITS = 1024
E3_LVL = 3
E3_RAD_MULT = 1.7
E3_TOPK = 10
RMS_THRES = 0.5
E3_MAX_CONFS = 200
E3_MMFF_MINIMIZATION_ITERS = 500
RMSD_THRESH = 1.0
CPU = os.cpu_count() - 2

IDX_WRONG_SMILES = pd.Series(['2147443170','2147443172','2147443173','2147443174','2147443175','2147443176','2147443177','2147443178','2147443179', 
                              '2147443180','2147443181','2147443182','2147443183','2147443184','2147443185','2147443186','2147443187','2147443188', 
                              '2147443189','2147443190','2147443192','2147443193','2147443194','2147443195','2147443197','2147443198','2147443199',  
                              '2147443200','2147443202','2147443203','2147443205','2147443206','2147443208','2147443209','2147443210','2147443213', 
                              '2147443214','2147443215','2147443216','2147443218','2147443219','2147443222','2147443223','2147443224','2147443225', 
                              '2147443226','2147443227','2147443228','2147443229','2147443230','2147443231','2147443232','2147443233','2147443234', 
                              '2147443235','2147443236','2147443238','2147443239','2147443240','2147443244','2147443245','2147443248','2147443249', 
                              '2147443250','2147443251','2147443254'], dtype=int, name='id')


def compute_torsional_fp (smiles: str, bits: int=128):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    tt_generator = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=bits)

    ttfp_list = [list(np.array(tt_generator.GetFingerprint(mol)), dtype=np.int8) for mol in mols]

    ttfp_df = pd.DataFrame(data= ttfp_list, columns= [f'ttfp_{i}' for i in range (bits)], dtype=np.int8)

    return ttfp_df


from rdkit.Chem.MolStandardize import rdMolStandardize
from e3fp.fingerprint.generate import fprints_dict_from_mol
from rdkit.Chem import rdMolTransforms, rdMolAlign

"""
def pruning_conformations (mol: Mol, energies: List, top_k: int):
    emin = energies[0][1]
    ewindow = 3.0  # kcal/mol
    # FILTER by energy 
    cands = [cid for cid, e in energies if (e - emin) <= ewindow]
    lprint(ll.REPORT, f'\tConformers Candidates: {cands}')

    drop = [cid for cid, e in energies if cid not in cands]
    for cid in drop:
        mol.RemoveConformer(cid)
    
    lprint(ll.REPORT, f'\tNum of conformers (after filtering by free energy) = {mol.GetNumConformers()}')  
    
    keep = [cands[0]]
    for i in range (1, len(cands)):
        is_far = True
        for kid in keep:
            rms = rdMolAlign.CalcRMS(prbMol= mol, refMol= mol, prbId= cands[i], refId= kid) # computing the similarity between conformers 
            if rms < RMSD_THRESH:
                is_far = False
                break
        if is_far:
            keep.append(cands[i])
            print(keep)
    keep = list(set(keep))

    drop = [cid for cid in cands if cid not in keep]
    for cid in drop:
        mol.RemoveConformer(id= cid)
    lprint(ll.REPORT, f'\tNum of conformers (after similarity analysis) = {mol.GetNumConformers()}')
    
    if mol.GetNumConformers() > 10:
        drop = keep[top_k:]
        for cid in drop:
            mol.RemoveConformer(id= cid)

    return mol
"""

def pruning_conformations (mol: Mol, energies: List, top_k: int):
    emin = energies[0][1]
    ewindow = 3.0  # kcal/mol
    # FILTER by energy 
    cands = [cid for cid, e in energies if (e - emin) <= ewindow]
    lprint(ll.REPORT, f'\tConformers Candidates: {cands}')

    drop = [cid for cid, e in energies if cid not in cands]
    for cid in drop:
        mol.RemoveConformer(cid)

    lprint(ll.REPORT, f'\tNum of conformers (after filtering by free energy) = {mol.GetNumConformers()}')  
    
    if mol.GetNumConformers() > 10:
        rms_matrix = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
        clusters = Butina.ClusterData(rms_matrix,mol.GetNumConformers(),distThresh = RMSD_THRESH, isDistData = True, reordering = True)

        energy_map = dict(energies)
        keep = []
        for cluster in clusters:
            best = min(cluster, key=lambda cid: energy_map[cid])
            keep.append(best)
        
        drop = [cid for cid in cands if cid not in keep]
        for cid in drop:
            mol.RemoveConformer(id= cid)

        lprint(ll.REPORT, f'\tNum of conformers (after similarity analysis) = {mol.GetNumConformers()}')
    
    if mol.GetNumConformers() > 10:
        drop = keep[top_k:]
        for cid in drop:
            mol.RemoveConformer(id= cid)

    return mol


def compute_e3fp_from_smiles(ids: List[str], smiles: List[str], tors_list: Optional[List]=None, bits: int=E3_BITS, level: int=E3_LVL, radius_multiplier: float=E3_RAD_MULT, top_k: int=E3_TOPK):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    e3fp_dict = {}

    for i, mol in enumerate(mols):
        lprint(ll.REPORT, f'{smiles[i]}')
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.pruneRmsThresh = RMS_THRES

        lprint(ll.INFO, f'\tComputing the optimal number of conformers to create')
        if tors_list is None:
            num_confs = min(E3_MAX_CONFS, E3_TOPK + 3*(rdMolDescriptors.CalcNumRotatableBonds(mol))) # the coeff '3' is so setted in order to literature settings
        else:
            num_confs = min(E3_MAX_CONFS, E3_TOPK + 3*tors_list[i])
        lprint(ll.REPORT, f'\tNumber of conformers = {num_confs}')
        
        lprint(ll.INFO, f'\tComputing the 3D conformers')
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

        mol_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')

        energies = []
        lprint(ll.INFO, f'\tComputing Free Energy of conformers')
        for cid in conf_ids:
            if mol_props is not None: 
                ff = AllChem.MMFFGetMoleculeForceField(mol, mol_props, confId=cid)
            else:
                ff = AllChem.UFFGetMoleculeForceField(mol)

            mol_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
            ff.Minimize(maxIts=E3_MMFF_MINIMIZATION_ITERS)
            energies.append((cid, ff.CalcEnergy()))
        energies.sort(key=lambda x: x[1])

        lprint(ll.INFO, f'\tPruning conformers...')
        mol_tokeep = pruning_conformations (mol= mol, energies= energies, top_k= top_k)
        lprint(ll.REPORT, f'\tBest conformers found: {mol_tokeep.GetNumConformers()}')

        fp_dict = fprints_dict_from_mol(
            mol, bits=bits, level=level, radius_multiplier=radius_multiplier,
            first=-1, stereo=True, counts=False, include_disconnected=True,
            rdkit_invariants=False, remove_duplicate_substructs=True
        )

        confs_fps = fp_dict[3]
        fp_union = set(confs_fps[0].indices)
        for cfp in confs_fps[1:]:
            fp = set(cfp.indices)
            fp_union.union(fp)
        fp_union = list(fp_union)

        bits_fp = np.zeros(E3_BITS, dtype=int)
        for id in fp_union:
            bits_fp[id] = 1
        e3fp_dict[ids[i]] = bits_fp

        
    e3fp_df = pd.DataFrame(e3fp_dict, columns=[f'e3fp_{j+1}' for j in range (E3_BITS)], dtype=np.int8)
        
    return e3fp_df



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



if __name__ == "__main__":
    train_df = pd.read_csv('./data/train_all.csv', header=0)
    train_df = train_df[train_df['id'].isin(IDX_WRONG_SMILES)]

    lprint(ll.STEP, f'STEP 2: Validating SMILES')
    train_df, _, _= validate_smiles.main(df1= train_df)
    lprint(ll.SUCCESS, f'STEP 2: Completed')
    train_df.to_csv('./data/clean.csv', index=False)
    
    print (train_df)

    lprint(ll.STEP, f'STEP 3: Computing Morgan + MACCKeys Fingerprint')
    mmfp_df = compute_morgan_maccs_fp (ismi_df= train_df, phase='aaaa')
    mmfp_df.to_csv('./data/aaaa_mmfp.csv', index=False)


    
    #lprint(ll.STEP, f'STEP 3: Computing E3 Fingerprint')
    #e3fp_df = compute_e3fp_from_smiles(ids= ids, smiles= smiles)
    #e3fp_df.to_csv('./data/train_e3fp.csv', index=False)
    #lprint(ll.SUCCESS, f'STEP 3: Completed')
       