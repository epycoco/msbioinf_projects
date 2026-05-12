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


####################################################
# CREA UNA LISTA GLOBALE CON DESCRITTORI DI RDKIT
# CREA UNA LISTA GLOBALE CON DESCRITTORI DI MORDRED
####################################################

def predict_from_smiles(test_df: pd.DataFrame, target: str, not3d: bool, model_dir="model"):
    """
    Generates predictions for a list of SMILES using pretrained models.

    Parameters:
        smiles_list (list): List of SMILES strings.
        target (str): Name of the target column to identify the model.
        model_dir (str): Directory where models are saved (default: 'model').

    Returns:
        tuple: (predictions, valid_smiles) Predictions for valid SMILES and list of valid SMILES.
    """
    lprint(ll.INFO, f"Generating predictions for {len(test_df.shape[0])} SMILES for target: {target}")
    
    idx = test_df['id']
    smiles = test_df['SMILES']

    """# Load model metadata and models
    metadata = load_best_model(target, model_dir)
    ensemble_method = metadata.get("ensemble_method")
    models = metadata.get("models")
    if not ensemble_method or not models:
        lprint(ll.ERROR, "Invalid metadata: ensemble_method or models missing")
        raise"""
    
    model_path = "./model/Tg/MAE/model_LGB_Tg.pkl"
    model = joblib.load(model_path)
    target_feature = list(model.feature_name())

    fp_feat = [f"fp_{i+1}" for i in range (FP_BITS+167)]
    chemberta_feat = pd.read_csv('./data/chemberta.csv', header=0).columns.to_list()
    
    desc_feat = [d for d in target_feature if d not in (fp_feat + chemberta_feat)]

    if desc_feat == []:
        # Convert SMILES to features
        lprint(ll.INFO, "Converting SMILES to features...")
        fp_df = compute_finger(smiles.to_list(), phase='test')
        ifp_df = pd.concat([idx, fp_df], axis=1)

        berta_df = compute_chemberta(smiles, phase='test')
        iberta_df = pd.concat([idx, berta_df], axis=1)

        test_all_df = pd.merge(test_df, ifp_df, on='id')
        test_all_df = pd.merge(test_all_df, iberta_df, on='id')

    else:
        return None

    # Create feature DataFrame
    X_test = test_all_df
    X_test.columns = X_test.columns.astype(str)
    X_test = X_test[sorted(X_test.columns)]

    # Load model metadata and models
    metadata = load_best_model(target, model_dir)
    ensemble_method = metadata.get("ensemble_method")
    models = metadata.get("models")
    if not ensemble_method or not models:
        lprint(ll.ERROR, "Invalid metadata: ensemble_method or models missing")
        raise

    # Generate predictions from individual models
    # Use top 2 models, consistent with train
    top_models = list(models.keys())
    lprint(ll.INFO, f"Using top 2 models for ensemble: {top_models}")
    test_list = [models[m].predict(X_test) for m in top_models]

    # Apply ensemble method
    if ensemble_method == "MeanEnsemble":
        lprint(ll.INFO, "Applying MeanEnsemble")
        predictions = np.mean(test_list, axis=0)
    elif ensemble_method == "WeightedEnsemble":
        lprint(ll.INFO, "Applying WeightedEnsemble")
        weights = metadata.get("weights", [1.0 / len(top_models)] * len(top_models))
        predictions = np.average(test_list, axis=0, weights=weights)
    elif ensemble_method == "StackingEnsemble":
        lprint(ll.INFO, "Applying StackingEnsemble")
        meta_model = metadata.get("meta_model")
        if meta_model is None:
            lprint(ll.ERROR, "No meta_model found in metadata for StackingEnsemble")
            raise ValueError("No meta_model found in metadata")
        X_meta_test = np.vstack(test_list).T
        predictions = meta_model.predict(X_meta_test)
    else:
        lprint(ll.ERROR, f"Unknown ensemble method: {ensemble_method}")
        raise

    lprint(ll.REPORT, f"Generated {len(predictions)} predictions for valid SMILES")
    return predictions, valid_smiles