import os, sys
import pandas as pd
import numpy as np

from src.config.config import *
from src.logging import lprint, LoggingLevels as ll

def is_numeric(series: pd.Series) -> bool:
    """
    Function that verify if the Series contains only numeric values
    """
    s = series.dropna().astype(str)
    is_all_numeric = pd.to_numeric(s, errors='coerce').notna().all()
    return is_all_numeric 

def base_filter_desc(idesc_df: pd.DataFrame) -> pd.DataFrame:
    idx = idesc_df['id']
    n_mols = len(idx)

    desc_df = idesc_df.drop(columns='id')
    valid_desc = desc_df.columns.to_list()
    try:
        valid_desc = [col for col in desc_df.columns if is_numeric(desc_df[col])]
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

        # If the number of repeted values is greater than 90% of total molecules, select the descriptors in order to drop it
        elif itd_rep >= REP_THRESH:
            desc_to_drop.append(desc)
    
    # Drop not valid desc and save the index of not valid molecules
    if desc_to_drop is not None:
        idesc_df = idesc_df.drop(columns=desc_to_drop)
    
    return idesc_df


def build_desc_df (desc_dir):
    mrd_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS./data/train_mordred.csv", header=0, low_memory=False)
    rdk2_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/train_rdkit2.csv", header=0)
    rdk3_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/train_rdkit3.csv", header=0)

    idesc_df = pd.merge(mrd_df, rdk2_df, on='id')
    idesc_df = pd.merge(idesc_df, rdk3_df, on='id')

    print(f'idesc_df shape before filtering: {idesc_df.shape}')

    idesc_df = base_filter_desc(idesc_df= idesc_df)
    idesc_df.to_csv(desc_dir, index=False)

    return idesc_df


def main():
    desc_dir = "C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/train_desc.csv"
    if not os.path.exists(desc_dir):
        idesc_df = build_desc_df(desc_dir)

    else:
        idesc_df = pd.read_csv(desc_dir)
        print('idesc dataframe found: UPLOADED')
    
    print (f'Final idesc DataFrame shape: {idesc_df.shape}')

    
main()