import pandas as pd
import numpy as np
from typing import List, Optional, Tuple

from src.tasks.compute_df import smiles_to_combined_fingerprints_with_descriptors
from src.tasks.io import load_best_model
from src.logging import lprint, LoggingLevels as ll

from src.config.fingerprint import compute_finger
from src.config.chemberta import compute_chemberta
from src.config.descriptors import compute_all_desc
from src.config.config import *

import joblib

def predict_from_smiles(test_df: pd.DataFrame, target: str, model_dir: str):
    """
    Generates predictions for a list of SMILES using pretrained models.

    Parameters:
        smiles_list (list): List of SMILES strings.
        target (str): Name of the target column to identify the model.
        model_dir (str): Directory where models are saved (default: 'model').

    Returns:
        tuple: (predictions, valid_smiles) Predictions for valid SMILES and list of valid SMILES.
    """
    lprint(ll.INFO, f"Generating predictions for {test_df.shape[0]} SMILES for target: {target}")
    
    idx = test_df['id']
    smiles = test_df['SMILES']

    """# Load model metadata and models
    metadata = load_best_model(target, model_dir)
    ensemble_method = metadata.get("ensemble_method")
    models = metadata.get("models")
    if not ensemble_method or not models:
        lprint(ll.ERROR, "Invalid metadata: ensemble_method or models missing")
        raise"""
    
    model = joblib.load(model_dir)
    target_feature = pd.read_csv('./data/target_train/Tg_train.csv', header=0).columns.to_list()

    fp_feat = [f"fp_{i+1}" for i in range (FP_BITS+167)]
    chemberta_feat = pd.read_csv('./data/train_chemberta.csv', header=0).columns.to_list()
    
    desc_feat = [d for d in target_feature if d not in (fp_feat + chemberta_feat + ['id', 'Tg'])]
    lprint(ll.WARN, f'Descriptors Feature: {desc_feat}\n')
    lprint(ll.INFO, "Converting SMILES to features...")
    if desc_feat == []:
        lprint(ll.INFO, f"No descriptors feature for Tg")
        # Convert SMILES to features
        
        fp_df = compute_finger(smiles.to_list(), phase='test')
        ifp_df = pd.concat([idx, fp_df], axis=1)

        berta_df = compute_chemberta(smiles, phase='test')
        iberta_df = pd.concat([idx, berta_df], axis=1)

        test_all_df = pd.merge(test_df, ifp_df, on='id')
        test_all_df = pd.merge(test_all_df, iberta_df, on='id')

    else:
        ## FORSE É MEGLIO CALCOLARE TUTTO SUL TEST E POI ANDARLO A PRENDERE
        fp_df = compute_finger(smiles.to_list(), phase='test')
        ifp_df = pd.concat([idx, fp_df], axis=1)

        berta_df = compute_chemberta(smiles, phase='test')
        iberta_df = pd.concat([idx, berta_df], axis=1)

        test_all_df = pd.merge(test_df, ifp_df, on='id')
        test_all_df = pd.merge(test_all_df, iberta_df, on='id')

    features = target_feature.copy()
    features.remove('Tg')
    X_test = test_all_df[features]
    print(X_test.columns.to_list()==features)

    print()
    lprint(ll.STEP, 'Prediction of Tg...')
    try:
        y_pred = model.predict(X_test)
        print ("Prediction for Tg:", y_pred)
    except Exception as e:
        lprint(ll.ERROR, f'{e}')
        
    return y_pred

    


if __name__ == "__main__":

    test_path = "./data/test_clean.csv"
    test_df = pd.read_csv(test_path, header=0)

    model_dir = "./model/Tg/MAE_f/model_LGB_Tg.pkl"
    
    y_pred = predict_from_smiles(test_df=test_df, target='Tg', model_dir=model_dir)
