import os
import joblib
from src.logging import lprint, LoggingLevels as ll
import pandas as pd

def save_best_model(models, model_name, target, weights=None, meta_model=None, scaler =None):
    """
    Saves individual models and ensemble metadata to the 'model' directory.

    Parameters:
        models (dict): Dictionary of trained models (e.g., {'XGB': xgb_model, ...}).
        model_name (str): Name of the model or ensemble method (e.g., 'MeanEnsemble').
        target (str): Name of the target column.
        weights (list, optional): Weights for WeightedEnsemble.
        meta_model (object, optional): Meta-model for StackingEnsemble.

    Returns:
        None
    """
    lprint(ll.INFO, f"Saving model: {model_name} for target: {target}")

    # Define the output directory
    output_dir = "model"
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, f"best_model_{target}_{model_name}.pkl")

    # Prepare metadata
    metadata = {
        "model_name": model_name,
        "target": target,
        "timestamp": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        "ensemble_method": model_name,
        "models": models,  # Save individual models
        "weights": weights,  # Save weights for WeightedEnsemble
        "meta_model": meta_model,  # Save meta-model for StackingEnsemble
        "scaler" : scaler
    }

    try:
        # Save metadata and models to disk
        joblib.dump(metadata, filename)
        lprint(ll.INFO, f"Successfully saved {model_name} for target {target} to {filename}")
    except Exception as e:
        lprint(ll.ERROR, f"Failed to save {model_name} for target {target}: {str(e)}")
        raise

def load_best_model(target, model_dir="model"):
    """
    Loads metadata and pretrained models from the specified directory.

    Parameters:
        target (str): Name of the target column to identify the model.
        model_dir (str): Directory where models are saved (default: 'model').

    Returns:
        dict: Metadata and loaded models (contains model_name, ensemble_method, models, weights, meta_model).
    """
    lprint(ll.INFO, f"Loading model for target: {target}")

    # Find the model file
    model_files = [f for f in os.listdir(model_dir) if f.startswith(f"best_model_{target}_") and f.endswith(".pkl")]
    if not model_files:
        lprint(ll.ERROR, f"No model file found for target {target} in {model_dir}")
        raise FileNotFoundError(f"No model found for target {target}")
    if len(model_files) > 1:
        lprint(ll.WARN, f"Multiple model files found for target {target}. Using the first: {model_files[0]}")

    model_file = os.path.join(model_dir, model_files[0])

    try:
        # Load metadata and models
        metadata = joblib.load(model_file)
        lprint(ll.INFO, f"Loaded metadata and models from {model_file}: {metadata['model_name']}")
        return metadata
    except Exception as e:
        lprint(ll.ERROR, f"Failed to load model {model_file}: {str(e)}")
        raise


import pickle
from glob import glob

def save_best_config(best_params, target, fold):
    output_dir = "best_params"
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, f"config_{target}_f{fold}.pkl")    
    with open(filename, 'wb') as f:
        pickle.dump(best_params, f, protocol=pickle.HIGHEST_PROTOCOL)

def load_all_best_configs(output_dir="best_params"):
    configs = {}
    files = glob(os.path.join(output_dir, "config_*.pkl"))
    
    for file in files:
        with open(file, 'rb') as f:
            configs[os.path.basename(file)] = pickle.load(f)
    return configs


def print_df_shape(df: pd.DataFrame, target: str):
    lprint(ll.REPORT, f'{target} DataFrame shape: {df.shape}')