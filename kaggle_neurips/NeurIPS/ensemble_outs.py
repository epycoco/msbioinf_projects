# General
import os
from pathlib import Path
from typing import Dict, Any, Optional, Callable, List, Tuple
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from src.logging import lprint, LoggingLevels as ll
from src.tasks.io import save_best_model, save_best_config
from src.tasks.train import train_with_optuna, get_model

# Train and Tuning
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold
from scipy.optimize import minimize
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# Model pool
from sklearn.linear_model import Ridge
from catboost import CatBoostRegressor
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
import lightgbm as lgb
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import ElasticNet, BayesianRidge, HuberRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import GradientBoostingRegressor

from sklearn.ensemble import StackingRegressor
from sklearn.base import clone

# Config and OptunaTuner
from src.config.config import *
import joblib

BERTA_FEAT = [f"chemberta_{i}" for i in range (BERTA_BITS)]
FP_FEAT = [f"mmfp_{i+1}" for i in range (FP_BITS+167)]
TARGETS = ['Density', 'Tg', 'FFV', 'Tc', 'Rg']

def loss_function(y_true, y_pred, mode: str, alpha: Optional[float] = None):
    if y_pred is None:
        return None

    if mode == "MSE":
        return mean_squared_error(y_true, y_pred)
    elif mode == "MAE":
        return mean_absolute_error(y_true, y_pred)
    elif mode == "CUSTOM":
        mse = mean_squared_error(y_true, y_pred)
        mae = mean_absolute_error(y_true, y_pred)
        custom = alpha * mse + (1 - alpha) * mae
        return custom
    else:
        raise ValueError("Loss non supportata: scegli tra 'MSE', 'MAE', 'CUSTOM'")
    

def weighted_ensemble(y_pred_df: pd.DataFrame, y_true: pd.Series):
    n_models = y_pred_df.shape[1]
    init_weights = np.ones(n_models) / n_models  # pesi iniziali uniformi

    def loss(weights):
        weighted_preds = np.average(y_pred_df.values, axis=1, weights=weights)
        return mean_absolute_error(y_true, weighted_preds)

    constraints = {'type': 'eq', 'fun': lambda w: np.sum(w) - 1}
    bounds = [(0, 1)] * n_models

    result = minimize(loss, init_weights, method='SLSQP', bounds=bounds, constraints=constraints) #Sequential Least Squares Programming: algoritmo di ottimizzazione numerica general-purpose per problemi vincolati e non lineari

    weighted_preds = np.average(y_pred_df.values, axis=1, weights=result.x)
    return weighted_preds, result.x


def stacking_ensemble_cv (y_pred_train_df: pd.DataFrame,
                          y_train: pd.Series,
                          y_pred_test_df: pd.DataFrame,
                          y_validation: pd.Series,
                          model_type: str) -> Tuple[np.ndarray, np.ndarray, Ridge]:
    
    best_params, best_value = train_with_optuna(model_type=model_type, X=y_pred_train_df, y=y_train, loss_mode='MAE', alpha=None)
    stacking_model = get_model(model_type=model_type, params=best_params)
    stacking_model.fit(y_pred_train_df, y_train)

    # Predizione finale sul test
    y_pred = stacking_model.predict(y_pred_test_df)
    print("R2 stacking:", r2_score(y_validation, y_pred))
    return stacking_model, y_pred
    


def build_ensemble (best_models: dict[str: (XGBRegressor | lgb.LGBMRegressor | CatBoostRegressor | 
                                     HistGradientBoostingRegressor | ExtraTreesRegressor | 
                                     RandomForestRegressor | MLPRegressor | ElasticNet | 
                                     BayesianRidge | HuberRegressor | KNeighborsRegressor | 
                                     GradientBoostingRegressor)],
                    X_train: pd.DataFrame,
                    y_train: pd.Series,
                    X_validation: pd.DataFrame,
                    y_validation: pd.Series,
                    ens_mean: bool = True,
                    ens_median: bool = True,
                    ens_weight: bool = True,
                    ens_stacking: bool = True) -> dict[str: Any]:

    # --- 2. predizioni dei top-2 su training e test ---
    # Predizioni sul validation set
    y_pred_train_df = pd.DataFrame({
        model_name: best_models[model_name].predict(X_train)
        for model_name in best_models
    })

    # Prediction on validation 
    if X_validation is not None:
        y_pred_valid_df = pd.DataFrame({
            model_name: best_models[model_name].predict(X_validation)
            for model_name in best_models
        })

    # --- Mean ---
    if ens_mean:
        mean_preds = y_pred_train_df.mean(axis=1)
    else:
        mean_preds = None

    # --- Median ---
    if ens_median:
        median_preds = y_pred_train_df.median(axis=1)

    # --- Weighted ---
    if ens_weight:
        weighted_preds, weights = weighted_ensemble(y_pred_train_df, y_train)
    else:
        weighted_preds, weights = None, None

    # --- Stacking ---
    if ens_stacking and (X_validation is not None):
        stacking_model, stacking_preds = stacking_ensemble_cv(
            y_pred_train_df = y_pred_train_df.values,
            y_train = y_train.values,
            y_pred_test_df = y_pred_valid_df.values,
            y_validation = y_validation.values,
            model_type = 'RFF'
        )
    else: 
        stacking_model, stacking_preds = None, None
    
    return {
        "mean_preds": mean_preds,
        "median_preds": median_preds,
        "weighted_preds": weighted_preds,
        "weights": weights,
        "stacking_model": stacking_model,
        "stacking_preds": stacking_preds
    }


def ensemble_preds(y_pred_test_df: pd.DataFrame, ens_type: str, weights: np.array, stack_model):
    if ens_type == "mean":
        return y_pred_test_df.mean(axis=1)
    
    elif ens_type == "median":
        return y_pred_test_df.median(axis=1)
    
    elif ens_type == "weighted":
        if weights is None:
            raise ValueError("weights must be provided for weighted ensemble")
        return (y_pred_test_df.values @ weights).flatten() # y = w1 * pred1 + w2 * pred2  (...)
    
    elif ens_type == "stacking":
        if stack_model is None:
            raise ValueError("stack_model must be provided for stacking")
        return stack_model.predict(y_pred_test_df.values)
    
    raise ValueError(f"Unknown ensemble type: {ens_type}")



def main ():
    test_y_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_clean.csv", header=0)
    test_fp_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_fp.csv", header=0)
    test_berta_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_chemberta.csv", header=0)
    test_desc_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_desc.csv", header=0)
    test_X_df = pd.merge(test_fp_df, test_berta_df, on='id')
    test_X_df = pd.merge(test_X_df, test_desc_df, on='id')
    X_validation = test_X_df.drop(columns='id')
    
    best_diz = {
        'Tg': ['ET', 'XGB'],
        'FFV': ['HGB', 'GBR'],
        'Tc': ['HGB', 'GBR'],
        'Density': ['HGB', 'GBR'],
        'Rg': ['ET', 'HGB']
    }

    for T in TARGETS:
        best1 = joblib.load(f"C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/model/{T}/MAE/model_{best_diz[T][0]}_{T}.pkl")
        best2 = joblib.load(f"C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/model/{T}/MAE/model_{best_diz[T][1]}_{T}.pkl")

        best_models= {
            f'{best_diz[T][0]}': best1,
            f'{best_diz[T][1]}': best2
        }

        target_df = pd.read_csv(f"C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/target_train/{T}_train.csv")
        X_train = target_df.drop(columns=['id', T])
        y_train = target_df[T]

        desc_list = X_train.drop(columns=[*FP_FEAT, *BERTA_FEAT]).columns.to_list()
        if desc_list is not None:
            X_validation = X_validation[[*FP_FEAT, *BERTA_FEAT, *desc_list]]

        y_test = test_y_df[T]
        y_pred_test_df = pd.DataFrame({
            model_name: best_models[model_name].predict(X_validation) 
            for model_name in best_models
        })
        
        ensamble_results = build_ensemble(best_models=best_models, X_train=X_train, y_train=y_train, X_validation=None, y_validation=None)

        # Evaluating Ensamble...
        ensemble_scores = {
            "mean": loss_function(y_train, ensamble_results["mean_preds"], mode='MAE'),
            "median": loss_function(y_train, ensamble_results["median_preds"], mode='MAE'),
            "weighted": loss_function(y_train, ensamble_results["weighted_preds"], mode='MAE'),
            "stacking": loss_function(y_train, ensamble_results["stacking_preds"], mode='MAE')
        }
        
        # ...and choose the best one
        valid_scores = {k: v for k, v in ensemble_scores.items() if v is not None}
        best_ensemble_type = min(valid_scores, key=valid_scores.get)
        lprint(ll.REPORT, f"Best ensemble: {best_ensemble_type} → {ensemble_scores[best_ensemble_type]:.4f}")

        final_preds = ensemble_preds(y_pred_test_df = y_pred_test_df,
                                     ens_type = best_ensemble_type, 
                                     weights = ensamble_results["weights"],
                                     stack_model=ensamble_results["stacking_model"])
        
        mae_val = mean_absolute_error(y_test, final_preds)

        ensemble_obj = {
            "type": best_ensemble_type,
            "MAE": mae_val,
            "preds": final_preds
        }
        
        lprint (ll.REPORT, f"{T} - {ensemble_obj}\n\n")


    
if __name__ == "__main__":
    main()
    