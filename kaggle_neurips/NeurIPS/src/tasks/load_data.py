import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, Fragments, Lipinski
from rdkit.Chem import rdmolops

from src.logging import lprint, LoggingLevels as ll


def main():
    train_df = pd.read_csv('data/neurips-open-polymer-prediction-2025/train.csv')
    test_df = pd.read_csv('data/neurips-open-polymer-prediction-2025/test.csv')
    lprint(ll.INFO, f"Training samples: {len(train_df)}")
    lprint(ll.INFO, f"Test samples: {len(test_df)}")
    

    # train supplement load
    lprint(ll.INFO, "Loading supplementary training data")
    sup0_tc = pd.read_csv('data/neurips-open-polymer-prediction-2025/train_supplement/dataset1.csv').rename(columns={'TC_mean': 'Tc'})
    sup1_tg = pd.read_csv('data/neurips-open-polymer-prediction-2025/train_supplement/dataset3.csv')
    sup2_ffv = pd.read_csv('data/neurips-open-polymer-prediction-2025/train_supplement/dataset4.csv')


    ##External Data Load
    lprint(ll.INFO, "Loading external data")
    sup3_tg = pd.read_csv('data/train_supports/TgSS_enriched_cleaned.csv',usecols=['SMILES', 'Tg'])
    sup4_tg = pd.read_csv('data/train_supports/JCIM_sup_bigsmiles.csv', usecols=['SMILES', 'Tg (C)']).rename(columns={'Tg (C)': 'Tg'})
    sup5_tg3 = pd.read_excel('data/train_supports/data_tg3.xlsx').rename(columns={'Tg [K]': 'Tg'})
    sup5_tg3['Tg'] = sup5_tg3['Tg'] - 273.15
    sup6_poly = pd.read_csv('data/train_supports/Tg_SMILES_class_pid_polyinfo_median.csv', usecols=['SMILES','Tg'])

    sup7_dnst = pd.read_excel('data/train_supports/data_dnst1.xlsx').rename(columns={'density(g/cm3)': 'Density'})[['SMILES', 'Density']]
    sup7_dnst = sup7_dnst[(sup7_dnst['SMILES'].notnull())&(sup7_dnst['Density'].notnull())&(sup7_dnst['Density'] != 'nylon')]
    sup7_dnst['Density'] = sup7_dnst['Density'].astype('float64')
    sup7_dnst['Density'] -= 0.118

    sup8 = pd.read_csv('data/train_supports/PI1070.csv', usecols=['SMILES', 'density', 'Rg', 'thermal_conductivity']).rename(columns={"density": "Density", "thermal_conductivity": "Tc"})
    sup9_polymer = pd.read_csv('data/train_supports/extended_polymer_dataset.csv')

    data_extend_dict = {
        "Tg": [1, 3, 4, 5, 6, 9],
        "FFV":[2],
        "Tc": [0, 8, 9],
        "Density": [7, 8, 9],
        "Rg": [8]
    }

    main_df = (train_df, test_df)
    ext_dfs = (
        sup0_tc, sup1_tg, sup2_ffv, sup3_tg, sup4_tg, sup5_tg3, sup6_poly,
        sup7_dnst, sup8, sup9_polymer
    )
    return main_df, ext_dfs, data_extend_dict