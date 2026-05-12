# General
import os
from pathlib import Path
from typing import Dict, Any, Optional, Callable, List
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from src.logging import lprint, LoggingLevels as ll
from src.tasks.io import save_best_model, save_best_config

# Train and Tuning
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold
from scipy.optimize import minimize
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# Model pool
from sklearn.linear_model import Ridge, Lasso
from catboost import CatBoostRegressor
from sklearn.ensemble import RandomForestRegressor
import xgboost as xgb
from xgboost import XGBRegressor
import lightgbm as lgb
from sklearn.svm import SVR
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import ElasticNet, BayesianRidge, HuberRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import GradientBoostingRegressor

# Config and OptunaTuner
from src.config.config import *
from src.tuning import OptunaTuner, create_regression_objective

import optuna, joblib
from optuna.samplers import TPESampler
from optuna.pruners import MedianPruner
optuna.logging.set_verbosity(optuna.logging.ERROR)  # Suppress Optuna's default logging


def suggest_params(model_type: str, trial: optuna.Trial, seed: int = 42) -> Dict[str, Any]:
    """Suggerimenti migliorati per i parametri dei modelli."""
    if model_type == "XGB":
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 100, 2000),  # Range più ragionevole
            'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3),
            'max_depth': trial.suggest_int('max_depth', 3, 12),  # Esteso per problemi complessi
            'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 100.0),  # Range esteso
            'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 100.0),
            'subsample': trial.suggest_float('subsample', 0.6, 1.0),  # Min più alto per stabilità
            'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
            'colsample_bylevel': trial.suggest_float('colsample_bylevel', 0.6, 1.0),
            'min_child_weight': trial.suggest_float('min_child_weight', 0.5, 20),  # Range esteso
            'gamma': trial.suggest_float('gamma', 0, 10),  # Range esteso
            'tree_method': 'hist',  # Fisso su hist per performance migliori
            'random_state': seed,  # Aggiunto per riproducibilità
            'n_jobs': -1
        }

    elif model_type == "LGB":
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 100, 3000),
            'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3),
            'max_depth': trial.suggest_int('max_depth', -1, 15),  # -1 per no limit
            'num_leaves': trial.suggest_int('num_leaves', 10, 300),  # Range esteso
            'min_child_samples': trial.suggest_int('min_child_samples', 5, 200),  # Min più basso
            'min_split_gain': trial.suggest_float('min_split_gain', 0, 1.0),  # Nome corretto
            'subsample': trial.suggest_float('subsample', 0.6, 1.0),
            'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 100.0),
            'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 100.0),
            'max_bin': trial.suggest_categorical('max_bin', [128, 255, 512]),  # Valori ottimali LGB
            'min_child_weight': trial.suggest_float('min_child_weight', 0.001, 50.0),
            'bagging_freq': trial.suggest_int('bagging_freq', 0, 10),
            'feature_fraction': trial.suggest_float('feature_fraction', 0.6, 1.0),  # Aggiunto
            'random_state': seed,
            'verbose': -1,  # Silenzia output
            'n_jobs': -1
        }
    elif model_type == "Cat":
        params = {
            'iterations': trial.suggest_int('iterations', 100, 2000),
            'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3),
            'depth': trial.suggest_int('depth', 4, 12),  # Range più ampio
            'l2_leaf_reg': trial.suggest_float('l2_leaf_reg', 1, 100.0),
            'bagging_temperature': trial.suggest_float('bagging_temperature', 0, 10),
            'random_strength': trial.suggest_float('random_strength', 0, 10),  # Aggiunto
            'border_count': trial.suggest_categorical('border_count', [32, 64, 128, 254]),  # Aggiunto
            'random_seed': seed,
            'verbose': False,
            'thread_count': -1,
            'task_type': 'GPU'
        }
    elif model_type == "RFF":
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 50, 1000),  # Range più ampio
            'max_depth': trial.suggest_int('max_depth', 3, 20),  # Più profondo
            'min_samples_split': trial.suggest_int('min_samples_split', 2, 20),
            'min_samples_leaf': trial.suggest_int('min_samples_leaf', 1, 10),
            'max_features': trial.suggest_categorical('max_features', ['sqrt', 'log2', None, 0.5, 0.8]),  # Più opzioni
            'bootstrap': trial.suggest_categorical('bootstrap', [True, False]),  # Aggiunto
            'min_impurity_decrease': trial.suggest_float('min_impurity_decrease', 0.0, 0.1),  # Aggiunto
            'random_state': seed
        }
    elif model_type == "ET":
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 50, 1000),
            'max_depth': trial.suggest_int('max_depth', 3, 25),
            'min_samples_split': trial.suggest_int('min_samples_split', 2, 20),
            'min_samples_leaf': trial.suggest_int('min_samples_leaf', 1, 15),
            'max_features': trial.suggest_categorical('max_features', ['sqrt', 'log2', None, 0.5, 0.8]),
            'bootstrap': trial.suggest_categorical('bootstrap', [True, False]),  # Aggiunto
            'min_impurity_decrease': trial.suggest_float('min_impurity_decrease', 0.0, 0.1),
            'random_state': seed
        }
    elif model_type == "HGB":
        params = {
            'max_iter': trial.suggest_int('max_iter', 100, 2000),
            'max_depth': trial.suggest_int('max_depth', 3, 20),  # Range esteso
            'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3),
            'min_samples_leaf': trial.suggest_int('min_samples_leaf', 5, 200),  # Range esteso
            'l2_regularization': trial.suggest_float('l2_regularization', 1e-8, 100.0),
            'max_leaf_nodes': trial.suggest_int('max_leaf_nodes', 10, 1000),  # Aggiunto
            'max_bins': trial.suggest_categorical('max_bins', [64, 128, 255]),  # Aggiunto
            'validation_fraction': trial.suggest_float('validation_fraction', 0.1, 0.3),  # Aggiunto
            'n_iter_no_change': trial.suggest_int('n_iter_no_change', 5, 20),  # Early stopping
            'random_state': seed
        }
    elif model_type == "MLP":
        params =  {
            'hidden_layer_sizes': trial.suggest_categorical('hidden_layer_sizes', 
                [(50,), (100,), (200,), (50, 50), (100, 50), (100, 100), (200, 100), (100, 50, 25)]),  # Più architetture
            'activation': trial.suggest_categorical('activation', ['relu', 'tanh']),  # Rimosso logistic (spesso problematico)
            'solver': trial.suggest_categorical('solver', ['adam', 'lbfgs']),  # lbfgs per dataset piccoli
            'alpha': trial.suggest_float('alpha', 1e-6, 1e-1),
            'learning_rate': trial.suggest_categorical('learning_rate', ['constant', 'adaptive']),  # Rimosso invscaling
            'learning_rate_init': trial.suggest_float('learning_rate_init', 1e-4, 1e-1),
            'batch_size': trial.suggest_categorical('batch_size', [32, 64, 128, 'auto']),  # Aggiunto
            'beta_1': trial.suggest_float('beta_1', 0.8, 0.999),  # Per Adam optimizer
            'beta_2': trial.suggest_float('beta_2', 0.9, 0.9999),  # Per Adam optimizer
            'early_stopping': True,  # Fisso
            'validation_fraction': 0.1,  # Fisso
            'random_state': seed
        }       
    elif model_type == "ElasticNet":
        params = {
            'alpha': trial.suggest_float('alpha', 1e-5, 100.0),  # Range esteso
            'l1_ratio': trial.suggest_float('l1_ratio', 0.01, 0.99),  # Evita estremi 0 e 1
            'fit_intercept': trial.suggest_categorical('fit_intercept', [True, False]),  # Aggiunto
            'max_iter': trial.suggest_int('max_iter', 1000, 10000),  # Aggiunto
            'selection': trial.suggest_categorical('selection', ['cyclic', 'random']),  # Aggiunto
            'random_state': seed
        }
    elif model_type == "BayesianRidge":
        params = {
            'alpha_1': trial.suggest_float('alpha_1', 1e-8, 1e-3),  # Range più specifico
            'alpha_2': trial.suggest_float('alpha_2', 1e-8, 1e-3),
            'lambda_1': trial.suggest_float('lambda_1', 1e-8, 1e-3),
            'lambda_2': trial.suggest_float('lambda_2', 1e-8, 1e-3),
            'fit_intercept': trial.suggest_categorical('fit_intercept', [True, False]),  # Aggiunto
        }
    elif model_type == "Huber":
        params =  {
            'epsilon': trial.suggest_float('epsilon', 1.1, 5.0),  # Range esteso
            'alpha': trial.suggest_float('alpha', 1e-6, 1.0),  # Range esteso
            'fit_intercept': trial.suggest_categorical('fit_intercept', [True, False]),  # Aggiunto
            'max_iter': trial.suggest_int('max_iter', 100, 2000)  # Aggiunto
        }
    elif model_type == "KNN":
        params = {
            'n_neighbors': trial.suggest_int('n_neighbors', 1, 30),  # Range ridotto (più ragionevole)
            'weights': trial.suggest_categorical('weights', ['uniform', 'distance']),
            'algorithm': trial.suggest_categorical('algorithm', ['auto', 'ball_tree', 'kd_tree', 'brute']),  # Aggiunto
            'leaf_size': trial.suggest_int('leaf_size', 10, 50),  # Aggiunto
            'p': trial.suggest_int('p', 1, 3),  # Esteso per includere Minkowski
            'metric': trial.suggest_categorical('metric', ['euclidean', 'manhattan', 'minkowski'])  # Aggiunto
        }
    elif model_type == "GBR":
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 50, 1500),
            'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3),
            'max_depth': trial.suggest_int('max_depth', 3, 12),
            'min_samples_split': trial.suggest_int('min_samples_split', 2, 20),
            'min_samples_leaf': trial.suggest_int('min_samples_leaf', 1, 15),
            'subsample': trial.suggest_float('subsample', 0.6, 1.0),
            'max_features': trial.suggest_categorical('max_features', ['sqrt', 'log2', None, 0.5, 0.8]),  # Aggiunto
            'min_impurity_decrease': trial.suggest_float('min_impurity_decrease', 1e-8, 0.1),  # Aggiunto
            'validation_fraction': trial.suggest_float('validation_fraction', 0.1, 0.3),  # Aggiunto
            'n_iter_no_change': trial.suggest_int('n_iter_no_change', 5, 20),  # Early stopping
            'random_state': seed
        }
    
    return params


def get_model(model_type: str, params: dict[str, Any]):
    """Return an instance of the specified model."""
    if model_type == "XGB":
        return XGBRegressor(**params)
    elif model_type == "LGB":
        return lgb.LGBMRegressor(**params)
    elif model_type == "Cat":
        return CatBoostRegressor(**params)  # Riduci verbosity
    elif model_type == "HGB":
        return HistGradientBoostingRegressor(**params)
    elif model_type == "ET":
        return ExtraTreesRegressor(**params)
    elif model_type == "RFF":
        return RandomForestRegressor(**params)
    elif model_type == "MLP":
        return MLPRegressor(**params, max_iter=MLP_ITER)
    elif model_type == "ElasticNet":
        return ElasticNet(**params)
    elif model_type == "BayesianRidge":
        return BayesianRidge(**params)
    elif model_type == "Huber":
        return HuberRegressor(**params)
    elif model_type == "KNN":
        return KNeighborsRegressor(**params)
    elif model_type == "GBR":
        return GradientBoostingRegressor(**params)
    else:
        raise ValueError(f"Model {model_type} not recognized.")
        

def loss_function(y_true, y_pred, mode: str, alpha: Optional[float] = None):
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


def suggest_best_params(model_type: str, trial: optuna.Trial, optim_params: dict[str, Any], seed: int = 42):
    if model_type == "XGB":
        optimized_params = {
            "n_estimators": trial.suggest_int("n_estimators", *optim_params["n_estimators"]),
            "learning_rate": trial.suggest_float("learning_rate", *optim_params["learning_rate"]),
            "max_depth": trial.suggest_int("max_depth", *optim_params["max_depth"]),
            "reg_lambda": trial.suggest_float("reg_lambda", *optim_params["reg_lambda"]),
            "reg_alpha": trial.suggest_float("reg_alpha", *optim_params["reg_alpha"]),
            "subsample": trial.suggest_float("subsample", *optim_params["subsample"]),
            "colsample_bytree": trial.suggest_float("colsample_bytree", *optim_params["colsample_bytree"]),
            "colsample_bylevel": trial.suggest_float("colsample_bylevel", *optim_params["colsample_bylevel"]),
            "min_child_weight": trial.suggest_float("min_child_weight", *optim_params["min_child_weight"]),
            "gamma": trial.suggest_float("gamma", *optim_params["gamma"]),
            "tree_method": "hist",
            "random_state": seed,
            "n_jobs": -1,
        }

    elif model_type == "LGB":
        optimized_params = {
            "n_estimators": trial.suggest_int("n_estimators", *optim_params["n_estimators"]),
            "learning_rate": trial.suggest_float("learning_rate", *optim_params["learning_rate"]),
            "max_depth": trial.suggest_int("max_depth", *optim_params["max_depth"]),
            "num_leaves": trial.suggest_int("num_leaves", *optim_params["num_leaves"]),
            "min_child_samples": trial.suggest_int("min_child_samples", *optim_params["min_child_samples"]),
            "min_split_gain": trial.suggest_float("min_split_gain", *optim_params["min_split_gain"]),
            "subsample": trial.suggest_float("subsample", *optim_params["subsample"]),
            "reg_lambda": trial.suggest_float("reg_lambda", *optim_params["reg_lambda"]),
            "reg_alpha": trial.suggest_float("reg_alpha", *optim_params["reg_alpha"]),
            "min_child_weight": trial.suggest_float("min_child_weight", *optim_params["min_child_weight"]),
            "bagging_freq": trial.suggest_int("bagging_freq", *optim_params["bagging_freq"]),
            "feature_fraction": trial.suggest_float("feature_fraction", *optim_params["feature_fraction"]),
            "max_bin": trial.suggest_categorical("max_bin", [128, 255, 512]),
            "random_state": seed,
            "verbose": -1,
            "n_jobs": -1,
        }

    elif model_type == "Cat":
        optimized_params = {
            "iterations": trial.suggest_int("iterations", *optim_params["iterations"]),
            "learning_rate": trial.suggest_float("learning_rate", *optim_params["learning_rate"]),
            "depth": trial.suggest_int("depth", *optim_params["depth"]),
            "l2_leaf_reg": trial.suggest_float("l2_leaf_reg", *optim_params["l2_leaf_reg"]),
            "bagging_temperature": trial.suggest_float("bagging_temperature", *optim_params["bagging_temperature"]),
            "border_count": trial.suggest_categorical("border_count", [32, 64, 128, 254]),
            "random_strength": trial.suggest_float("random_strength", *optim_params["random_strength"]),
            "random_seed": seed,
            "verbose": False,
            "thread_count": -1,
            'task_type': 'CPU'
        }

    elif model_type == "RFF":
        optimized_params = {
            "n_estimators": trial.suggest_int("n_estimators", *optim_params["n_estimators"]),
            "max_depth": trial.suggest_int("max_depth", *optim_params["max_depth"]),
            "min_samples_split": trial.suggest_int("min_samples_split", *optim_params["min_samples_split"]),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", *optim_params["min_samples_leaf"]),
            "max_features": trial.suggest_categorical("max_features", ["sqrt", "log2", None, 0.5, 0.8]),
            "bootstrap": trial.suggest_categorical("bootstrap", [True, False]),
            "min_impurity_decrease": trial.suggest_float("min_impurity_decrease", 0.0, 0.1),
            "random_state": seed,
        }

    elif model_type == "ET":
        optimized_params = {
            "n_estimators": trial.suggest_int("n_estimators", *optim_params["n_estimators"]),
            "max_depth": trial.suggest_int("max_depth", *optim_params["max_depth"]),
            "min_samples_split": trial.suggest_int("min_samples_split", *optim_params["min_samples_split"]),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", *optim_params["min_samples_leaf"]),
            "max_features": trial.suggest_categorical("max_features", ["sqrt", "log2", None, 0.5, 0.8]),
            "bootstrap": trial.suggest_categorical("bootstrap", [True, False]),
            "min_impurity_decrease": trial.suggest_float("min_impurity_decrease", 0.0, 0.1),
            "random_state": seed,
        }

    elif model_type == "HGB":
        optimized_params = {
            "max_iter": trial.suggest_int("max_iter", *optim_params["max_iter"]),
            "max_depth": trial.suggest_int("max_depth", *optim_params["max_depth"]),
            "learning_rate": trial.suggest_float("learning_rate", *optim_params["learning_rate"]),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", *optim_params["min_samples_leaf"]),
            "l2_regularization": trial.suggest_float("l2_regularization", *optim_params["l2_regularization"]),
            "max_leaf_nodes": trial.suggest_int("max_leaf_nodes", *optim_params["max_leaf_nodes"]),
            "max_bins": trial.suggest_categorical("max_bins", [64, 128, 255]),
            "validation_fraction": trial.suggest_float("validation_fraction", 0.1, 0.3),
            "n_iter_no_change": trial.suggest_int("n_iter_no_change", 5, 20),
            "random_state": seed,
        }

    elif model_type == "MLP":
        optimized_params = {
            "hidden_layer_sizes": optim_params["hidden_layer_sizes"],
            "activation": optim_params["activation"],
            "solver": optim_params["solver"],
            "learning_rate": optim_params["learning_rate"],
            "alpha": trial.suggest_float("alpha", *optim_params["alpha"]),
            "learning_rate_init": trial.suggest_float("learning_rate_init", *optim_params["learning_rate_init"]),
            "beta_1": trial.suggest_float("beta_1", *optim_params["beta_1"]),
            "beta_2": trial.suggest_float("beta_2", *optim_params["beta_2"]),
            "batch_size": trial.suggest_categorical("batch_size", [32, 64, 128, "auto"]),
            "early_stopping": True,
            "validation_fraction": 0.1,
            "random_state": seed,
        }

    elif model_type == "ElasticNet":
        optimized_params = {
            "alpha": trial.suggest_float("alpha", *optim_params["alpha"]),
            "l1_ratio": trial.suggest_float("l1_ratio", *optim_params["l1_ratio"]),
            "max_iter": trial.suggest_int("max_iter", *optim_params["max_iter"]),
            "selection": optim_params["selection"],
            "fit_intercept": trial.suggest_categorical("fit_intercept", [True, False]),
            "random_state": seed,
        }

    elif model_type == "BayesianRidge":
        optimized_params = {
            "alpha_1": trial.suggest_float("alpha_1", *optim_params["alpha_1"]),
            "alpha_2": trial.suggest_float("alpha_2", *optim_params["alpha_2"]),
            "lambda_1": trial.suggest_float("lambda_1", *optim_params["lambda_1"]),
            "lambda_2": trial.suggest_float("lambda_2", *optim_params["lambda_2"]),
            "fit_intercept": optim_params["fit_intercept"],
        }

    elif model_type == "Huber":
        optimized_params = {
            "epsilon": trial.suggest_float("epsilon", *optim_params["epsilon"]),
            "alpha": trial.suggest_float("alpha", *optim_params["alpha"]),
            "max_iter": trial.suggest_int("max_iter", *optim_params["max_iter"]),
            "fit_intercept": optim_params["fit_intercept"],
        }

    elif model_type == "KNN":
        optimized_params = {
            "n_neighbors": trial.suggest_int("n_neighbors", *optim_params["n_neighbors"]),
            "leaf_size": trial.suggest_int("leaf_size", *optim_params["leaf_size"]),
            "weights": optim_params["weights"],
            "algorithm": optim_params["algorithm"],
            "metric": optim_params["metric"],
            "p": trial.suggest_int("p", 1, 3),
        }

    elif model_type == "GBR":
        optimized_params = {
            "n_estimators": trial.suggest_int("n_estimators", *optim_params["n_estimators"]),
            "learning_rate": trial.suggest_float("learning_rate", *optim_params["learning_rate"]),
            "max_depth": trial.suggest_int("max_depth", *optim_params["max_depth"]),
            "min_samples_split": trial.suggest_int("min_samples_split", *optim_params["min_samples_split"]),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", *optim_params["min_samples_leaf"]),
            "subsample": trial.suggest_float("subsample", *optim_params["subsample"]),
            "min_impurity_decrease": trial.suggest_float("min_impurity_decrease", *optim_params["min_impurity_decrease"]),
            "n_iter_no_change": trial.suggest_int("n_iter_no_change", *optim_params["n_iter_no_change"]),
            "max_features": trial.suggest_categorical("max_features", ["sqrt", "log2", None, 0.5, 0.8]),
            "validation_fraction": trial.suggest_float("validation_fraction", 0.1, 0.3),
            "random_state": seed,
        }

    return optimized_params


KF = KFold(n_splits=5, shuffle=True, random_state=42)
def make_objective(model_type: str, X, y, loss_mode: str, alpha: float, optimized_params: Optional[dict] = None):
    def objective(trial: optuna.Trial):
        trial_seed = 42 + trial.number

        if not optimized_params:
            params = suggest_params(model_type=model_type, trial=trial, seed=trial_seed)
        else:
            params = suggest_best_params(model_type=model_type, trial=trial, optim_params=optimized_params, seed=trial_seed)

        fold_scores = []
        fold_r2_scores = []  # Lista per raccogliere R2 di ogni fold

        for fold, (train_idx, val_idx) in enumerate(KF.split(X), start=1):
            X_train, X_val = X[train_idx], X[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]

            model = get_model(model_type, params)
            model.fit(X_train, y_train)

            y_pred = model.predict(X_val)
            
            # Calcola loss per questo fold
            score = loss_function(y_val, y_pred, mode=loss_mode, alpha=alpha)
            fold_scores.append(score)
            
            # Calcola R2 per questo fold
            r2 = r2_score(y_val, y_pred)
            fold_r2_scores.append(r2)

        # Medie sui fold
        mean_score = np.mean(fold_scores)
        mean_r2 = np.mean(fold_r2_scores)

        lprint(ll.INFO, f'\t[Trial {trial.number+1}] hyperparameters = {params}')
        lprint(ll.REPORT, f"\t[Trial {trial.number+1} | seed={trial_seed}] Mean CV {loss_mode} = {mean_score:.4f} | R2 = {mean_r2:.4f}\n")
        
        return mean_score

    return objective



trial_results = []  # lista globale dove accumuliamo
def save_trials_callback(study, trial):
    trial_results.append({
        "number": trial.number,
        "params": trial.params,
        "value": trial.value
    })


def get_top_trials(trial_results, top_k=10):
    # Ordina in base a val_loss (trial["value"])
    sorted_trials = sorted(trial_results, key=lambda x: x["value"])
    return sorted_trials[:top_k]


def build_params_dict(trial_params: dict, model_type: str, seed: int = 42):
    if model_type == "XGB":
        params = {
            'n_estimators': trial_params.get('n_estimators'),
            'learning_rate': trial_params.get('learning_rate'),
            'max_depth': trial_params.get('max_depth'),
            'reg_lambda': trial_params.get('reg_lambda'),
            'reg_alpha': trial_params.get('reg_alpha'),
            'subsample': trial_params.get('subsample'),
            'colsample_bytree': trial_params.get('colsample_bytree'),
            'colsample_bylevel': trial_params.get('colsample_bylevel'),
            'min_child_weight': trial_params.get('min_child_weight'),
            'gamma': trial_params.get('gamma'),
            'tree_method': 'hist',
            'random_state': seed,
            'n_jobs': -1
        }

    elif model_type == "LGB":
        params = {
            'n_estimators': trial_params.get('n_estimators'),
            'learning_rate': trial_params.get('learning_rate'),
            'max_depth': trial_params.get('max_depth'),
            'num_leaves': trial_params.get('num_leaves'),
            'min_child_samples': trial_params.get('min_child_samples'),
            'min_split_gain': trial_params.get('min_split_gain'),
            'subsample': trial_params.get('subsample'),
            'reg_lambda': trial_params.get('reg_lambda'),
            'reg_alpha': trial_params.get('reg_alpha'),
            'max_bin': trial_params.get('max_bin'),
            'min_child_weight': trial_params.get('min_child_weight'),
            'bagging_freq': trial_params.get('bagging_freq'),
            'feature_fraction': trial_params.get('feature_fraction'),
            'random_state': seed,
            'verbose': -1,
            'n_jobs': -1
        }

    elif model_type == "Cat":
        params = {
            'iterations': trial_params.get('iterations'),
            'learning_rate': trial_params.get('learning_rate'),
            'depth': trial_params.get('depth'),
            'l2_leaf_reg': trial_params.get('l2_leaf_reg'),
            'border_count': trial_params.get('border_count'),
            'bagging_temperature': trial_params.get('bagging_temperature'),
            'random_strength': trial_params.get('random_strength'),
            'subsample': trial_params.get('subsample'),
            'random_seed': seed,
            'verbose': False,
            'thread_count': -1,
        }

    elif model_type == "RFF":
        params = {
            'n_estimators': trial_params.get('n_estimators'),
            'max_depth': trial_params.get('max_depth'),
            'min_samples_split': trial_params.get('min_samples_split'),
            'min_samples_leaf': trial_params.get('min_samples_leaf'),
            'max_features': trial_params.get('max_features'),
            'bootstrap': trial_params.get('bootstrap'),
            'min_impurity_decrease': trial_params.get('min_impurity_decrease'),
            'random_state': seed
        }

    elif model_type == "ET":
        params = {
            'n_estimators': trial_params.get('n_estimators'),
            'max_depth': trial_params.get('max_depth'),
            'min_samples_split': trial_params.get('min_samples_split'),
            'min_samples_leaf': trial_params.get('min_samples_leaf'),
            'max_features': trial_params.get('max_features'),
            'bootstrap': trial_params.get('bootstrap'),
            'min_impurity_decrease': trial_params.get('min_impurity_decrease'),
            'random_state': seed
        }

    elif model_type == "HGB":
        params = {
            'max_iter': trial_params.get('max_iter'),
            'max_depth': trial_params.get('max_depth'),
            'learning_rate': trial_params.get('learning_rate'),
            'min_samples_leaf': trial_params.get('min_samples_leaf'),
            'l2_regularization': trial_params.get('l2_regularization'),
            'max_leaf_nodes': trial_params.get('max_leaf_nodes'),
            'max_bins': trial_params.get('max_bins'),
            'validation_fraction': trial_params.get('validation_fraction'),
            'n_iter_no_change': trial_params.get('n_iter_no_change'),
            'random_state': seed
        }

    elif model_type == "MLP":
        params = {
            'hidden_layer_sizes': trial_params.get('hidden_layer_sizes'),
            'activation': trial_params.get('activation'),
            'solver': trial_params.get('solver'),
            'alpha': trial_params.get('alpha'),
            'learning_rate': trial_params.get('learning_rate'),
            'learning_rate_init': trial_params.get('learning_rate_init'),
            'batch_size': trial_params.get('batch_size'),
            'beta_1': trial_params.get('beta_1'),
            'beta_2': trial_params.get('beta_2'),
            'early_stopping': True,
            'validation_fraction': 0.1,
            'random_state': seed
        }

    elif model_type == "ElasticNet":
        params = {
            'alpha': trial_params.get('alpha'),
            'l1_ratio': trial_params.get('l1_ratio'),
            'fit_intercept': trial_params.get('fit_intercept'),
            'max_iter': trial_params.get('max_iter'),
            'selection': trial_params.get('selection'),
            'random_state': seed
        }

    elif model_type == "BayesianRidge":
        params = {
            'alpha_1': trial_params.get('alpha_1'),
            'alpha_2': trial_params.get('alpha_2'),
            'lambda_1': trial_params.get('lambda_1'),
            'lambda_2': trial_params.get('lambda_2'),
            'fit_intercept': trial_params.get('fit_intercept'),
        }

    elif model_type == "Huber":
        params = {
            'epsilon': trial_params.get('epsilon'),
            'alpha': trial_params.get('alpha'),
            'fit_intercept': trial_params.get('fit_intercept'),
            'max_iter': trial_params.get('max_iter')
        }

    elif model_type == "KNN":
        params = {
            'n_neighbors': trial_params.get('n_neighbors'),
            'weights': trial_params.get('weights'),
            'algorithm': trial_params.get('algorithm'),
            'leaf_size': trial_params.get('leaf_size'),
            'p': trial_params.get('p'),
            'metric': trial_params.get('metric')
        }

    elif model_type == "GBR":
        params = {
            'n_estimators': trial_params.get('n_estimators'),
            'learning_rate': trial_params.get('learning_rate'),
            'max_depth': trial_params.get('max_depth'),
            'min_samples_split': trial_params.get('min_samples_split'),
            'min_samples_leaf': trial_params.get('min_samples_leaf'),
            'subsample': trial_params.get('subsample'),
            'max_features': trial_params.get('max_features'),
            'min_impurity_decrease': trial_params.get('min_impurity_decrease'),
            'validation_fraction': trial_params.get('validation_fraction'),
            'n_iter_no_change': trial_params.get('n_iter_no_change'),
            'random_state': seed
        }

    else:
        raise ValueError(f"Model type {model_type} non supportato")

    return params


def find_optim_params (model_type: str, best_params_dict: list[dict[str, Any]]):
    ndict = len(best_params_dict)
    optim_params = {}
    if model_type == "XGB":
        params_modable = ["n_estimators","learning_rate","max_depth","reg_lambda","reg_alpha","subsample","colsample_bytree",
                          "colsample_bylevel","min_child_weight","gamma",]
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)
            
    elif model_type == "LGB":
        params_modable = ["n_estimators","learning_rate","max_depth","num_leaves","min_child_samples","min_split_gain",
                       "subsample","reg_lambda","reg_alpha","min_child_weight","bagging_freq","feature_fraction"]
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "Cat":
        params_modable = ["iterations","learning_rate","depth","l2_leaf_reg","bagging_temperature", "random_strength"]
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "RFF":
        params_modable = ['n_estimators', 'max_depth', 'min_samples_split', 'min_samples_leaf']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "ET":
        params_modable = ['n_estimators', 'max_depth', 'min_samples_split', 'min_samples_leaf']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "HGB":
        params_modable = ['max_iter', 'max_depth', 'learning_rate', 'min_samples_leaf', 'l2_regularization', 'max_leaf_nodes']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "MLP":
        optim_params['hidden_layer_sizes'] = best_params_dict[0]['hidden_layer_sizes']
        optim_params['activation'] = best_params_dict[0]['activation']
        optim_params['solver'] = best_params_dict[0]['solver']
        optim_params['learning_rate'] = best_params_dict[0]['learning_rate']
        params_modable = ['alpha', 'learning_rate_init', 'beta_1', 'beta_2']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "ElasticNet":
        optim_params['selection'] = best_params_dict[0]['selection']
        params_modable = ['alpha', 'l1_ratio', 'max_iter']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "BayesianRidge":
        optim_params['fit_intercept'] = best_params_dict[0]['fit_intercept']
        params_modable = ['alpha_1', 'alpha_2', 'lambda_1', 'lambda_2']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "Huber":
        optim_params['fit_intercept'] = best_params_dict[0]['fit_intercept']
        params_modable = ['epsilon', 'alpha', 'max_iter']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "KNN":
        optim_params['weights'] = best_params_dict[0]['weights']
        optim_params['algorithm'] = best_params_dict[0]['algorithm']
        optim_params['metric'] = best_params_dict[0]['metric']
        params_modable = ['n_neighbors', 'leaf_size']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    elif model_type == "GBR":
        params_modable = ['n_estimators', 'learning_rate', 'max_depth', 'min_samples_split', 'min_samples_leaf', 'subsample', 'min_impurity_decrease', 'n_iter_no_change']
        for p in params_modable:
            p_list = [best_params_dict[i][p] for i in range (ndict)]
            optim_params[p] = min(p_list), max(p_list)

    return optim_params


def train_with_optuna(model_type: str, X: pd.DataFrame, y: pd.Series, loss_mode: str, alpha: float, n_trials_phase1=20, n_trials_phase2=40):
    trial_results.clear()
    # Phase 1: Find best hyperparams
    lprint(ll.INFO, f'Optimization of hyperparameters for {model_type}...')
    study1 = optuna.create_study(direction="minimize", sampler=TPESampler(seed=42))
    study1.optimize(make_objective(model_type, X, y, loss_mode, alpha), n_trials=n_trials_phase1, show_progress_bar=False, callbacks=[save_trials_callback])
    
    lprint(ll.INFO, f'Choosing the best hyperparameters for {model_type}...')
    params_top10 = get_top_trials(trial_results, top_k=10)
    best_params_dict = [build_params_dict(trial["params"], model_type= model_type, seed=42) for trial in params_top10]
    optimized_params = find_optim_params (model_type= model_type, best_params_dict= best_params_dict)

    lprint(ll.REPORT, f'Optimized hyperparameters: {optimized_params}')

    # Phase 2: Find the best val_loss
    lprint(ll.INFO, f'Training with best hyperparameters...')
    study2 = optuna.create_study(direction="minimize", sampler=TPESampler(seed=42))
    study2.optimize(make_objective(model_type, X, y, loss_mode, alpha, optimized_params), n_trials=n_trials_phase2, show_progress_bar=False)
    
    return study2.best_params, study2.best_value


def train (X_train: pd.DataFrame, y_train: pd.Series, X_test: Optional[pd.DataFrame], target: str, loss_mode: str, alpha: Optional[float] = None) -> dict[str: Any]:
    # Train a model using a custom loss function (or MSE/MAE) but scoring the model on MAE
    best_models = {}
    best_values = {}

    X = X_train.values
    y = y_train.values

    if not alpha:
        model_dir = f"./model/{target}/{loss_mode}"
    else:
        model_dir = f"./model/{target}/{loss_mode}/alfa_{alpha}"
    
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    for m in MODELS:
        lprint(ll.INFO, f'{m}')
        if os.path.exists(f'{model_dir}/model_{m}_{target}.pkl'):
            model = joblib.load(f'{model_dir}/model_{m}_{target}.pkl')
            best_params = joblib.load(f'{model_dir}/params_{m}_{target}.pkl')
            lprint (ll.INFO, f'A model found per {m}: params_{m}_{target}.pkl and model_{m}_{target}.pkl Loaded')
            y_pred = model.predict(X)
            best_value = loss_function(y_true=y, y_pred=y_pred, mode=loss_mode)
            lprint (ll.REPORT, f"{loss_mode} = {best_value} | R2 = {r2_score(y_true=y, y_pred=y_pred)}")

        else:
            best_params, best_value = train_with_optuna(model_type=m, X=X, y=y, loss_mode=loss_mode, alpha=alpha)
            model = get_model(model_type= m, params= best_params)
            model.fit(X, y)
            joblib.dump(model, f"{model_dir}/model_{m}_{target}.pkl")
            joblib.dump(best_params, f"{model_dir}/params_{m}_{target}.pkl")

        best_models[m] = model
        best_values[m] = best_value
        print('\n')
    
    # Compute Ensamble
    ensamble_results = run_ensembles(best_models, best_values, X_train=X, y_train=y, X_test=X_test)

    # Evaluating Ensamble
    ensemble_scores = {
        "mean": loss_function(y, ensamble_results["mean_preds"], mode='MAE', alpha= alpha),
        "weighted": loss_function(y, ensamble_results["weighted_preds"], mode='MAE', alpha= alpha)
    }
    if ensamble_results["stacking_oof"]:
        ensemble_scores["stacking"] = loss_function(y, ensamble_results["stacking_oof"], mode='MAE', alpha= alpha)
    
    #   and choose the best one
    best_ensemble_type = min(ensemble_scores, key=ensemble_scores.get)
    lprint(ll.REPORT, f"Best ensemble: {best_ensemble_type} → {ensemble_scores[best_ensemble_type]:.4f}")

    if best_ensemble_type == "mean":
        final_preds = ensamble_results["mean_preds"]
    elif best_ensemble_type == "weighted":
        final_preds = ensamble_results["weighted_preds"]
    elif best_ensemble_type == "stacking":
        final_preds = ensamble_results["stacking_oof"]
    else:
        raise ValueError("Tipo ensemble non riconosciuto")
    
    mae_val = mean_absolute_error(y, final_preds)

    ensemble_obj = {
        "type": best_ensemble_type,
        "top2_models": ensamble_results["top2_names"],
        "mae": mae_val                                   
    }
    
    # --- Salvataggio solo ensemble migliore ---
    if best_ensemble_type == "mean":
        ensemble_obj["preds"] = ensamble_results["mean_preds"]

    elif best_ensemble_type == "weighted":
        ensemble_obj["preds"] = ensamble_results["weighted_preds"]
        ensemble_obj["weights"] = ensamble_results["weights"]

    elif best_ensemble_type == "stacking":
        ensemble_obj["model"] = ensamble_results["stacking_model"]
        ensemble_obj["oof_preds"] = ensamble_results["stacking_oof"]
        if ensamble_results["stacking_test"] is not None:
            ensemble_obj["test_preds"] = ensamble_results["stacking_test"]

    # unico file .pkl con tutto dentro
    joblib.dump(ensemble_obj, f"{model_dir}/best_ensemble_{target}.pkl")

    return ensemble_obj

        
def mean_ensemble(preds_list):
    preds_array = np.vstack(preds_list)  # shape: (n_models, n_samples)
    return np.mean(preds_array, axis=0)


def weighted_ensemble(preds_list, y_true):
    init_weights = [1 / len(preds_list)] * len(preds_list)
    def loss(weights):
        weighted_preds = np.average(preds_list, axis=0, weights=weights)
        return mean_absolute_error(y_true, weighted_preds)
    constraints = {'type': 'eq', 'fun': lambda w: 1 - sum(w)}
    bounds = [(0, 1)] * len(preds_list)
    result = minimize(loss, init_weights, method='SLSQP', bounds=bounds, constraints=constraints)
    return np.average(preds_list, axis=0, weights=result.x), result.x


def stacking_ensemble(preds_list, y_true, test_list):
    X_meta_train = np.vstack(preds_list).T
    X_meta_test = np.vstack(test_list).T
    meta_model = Ridge(alpha=1.0)
    meta_model.fit(X_meta_train, y_true)
    oof_preds = meta_model.predict(X_meta_train)
    test_preds = meta_model.predict(X_meta_test)
    return oof_preds, test_preds, meta_model


def run_ensembles(best_models: dict[str: (XGBRegressor | lgb.LGBMRegressor | CatBoostRegressor | 
                                     HistGradientBoostingRegressor | ExtraTreesRegressor | 
                                     RandomForestRegressor | MLPRegressor | ElasticNet | 
                                     BayesianRidge | HuberRegressor | KNeighborsRegressor | 
                                     GradientBoostingRegressor)], 
                    best_values: dict[str: float], 
                    X_train: pd.DataFrame, 
                    y_train: pd.DataFrame, 
                    X_test: Optional[pd.DataFrame] = None) -> dict[str: Any]:
    """
    Esegue ensemble (mean, weighted, stacking) sui 2 migliori modelli.

    best_models : dict {model_name: fitted_model}
    best_values : dict {model_name: score}
    X_train, y_train : dati per l'addestramento meta-modello / pesi
    X_test : opzionale, per avere predizioni finali sui dati di test
    """

    # --- 1. ordina modelli in base alle performance ---
    sorted_models = sorted(best_values.items(), key=lambda x: x[1])
    top2_names = [m for m, _ in sorted_models[:2]]
    top2_models = [best_models[m] for m in top2_names]

    lprint(ll.REPORT, f"Top 2 models: {top2_names}")

    # --- 2. predizioni dei top-2 su training ---
    preds_list = [model.predict(X_train) for model in top2_models]

    # --- 3. ensemble mean ---
    mean_preds = mean_ensemble(preds_list)

    # --- 4. ensemble pesato (ottimizzato su MAE) ---
    weighted_preds, weights = weighted_ensemble(preds_list, y_train)

    # --- 5. stacking (Ridge come meta-modello) ---
    if X_test is not None:
        test_list = [model.predict(X_test) for model in top2_models]
        oof_preds, test_preds, meta_model = stacking_ensemble(preds_list, y_train, test_list)
    else:
        oof_preds, test_preds, meta_model = None, None, None

    return {
        "mean_preds": mean_preds, 
        "weighted_preds": weighted_preds,
        "weights": weights,
        "stacking_test": test_preds,
        "stacking_oof": oof_preds,
        "stacking_model": meta_model,
        "top2_names": top2_names
    }

