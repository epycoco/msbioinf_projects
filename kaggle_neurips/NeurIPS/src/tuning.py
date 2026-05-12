import optuna
from optuna.samplers import TPESampler
from optuna.pruners import MedianPruner
import pandas as pd
import numpy as np
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import mean_absolute_error
from typing import Dict, Any, Optional, Callable, List
import pickle
from datetime import datetime
import warnings
import optuna.logging

warnings.filterwarnings('ignore')
optuna.logging.set_verbosity(optuna.logging.ERROR)  # Suppress Optuna's default logging

from src.config.config import *
from src.logging import lprint, LoggingLevels as ll

# models
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

class OptunaTuner:
    def __init__(self, sampler_type: str = 'tpe', pruner_type: str = 'median', study_name: Optional[str] = None):
        self.sampler_type = sampler_type
        self.pruner_type = pruner_type
        self.sampler = TPESampler(seed=SEED, n_startup_trials=10)
        self.pruner = MedianPruner(n_startup_trials=5, n_warmup_steps=10, interval_steps=1)
        self.study_name = study_name or f"study_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}"
        self.study = None
        self.best_params = None
        self.best_score = None
        self.optimization_history = []

    def create_study(self, direction: str = 'minimize'):
        self.study = optuna.create_study(
            direction=direction,
            sampler=self.sampler,
            pruner=self.pruner,
            study_name=self.study_name
        )
        lprint(ll.INFO, f"Study Created: {self.study_name}")
        lprint(ll.INFO, f"Sampler: {self.sampler_type}, Pruner: {self.pruner_type}")
        return self.study
    
    def get_model_instance(self, model_name, params):
        """Return an instance of the specified model."""
        if model_name == "XGB":
            return XGBRegressor(**params)
        elif model_name == "LGB":
            return lgb.LGBMRegressor(**params)
        elif model_name == "Cat":
            return CatBoostRegressor(**params, verbose=False)  # Riduci verbosity
        elif model_name == "HGB":
            return HistGradientBoostingRegressor(**params)
        elif model_name == "ET":
            return ExtraTreesRegressor(**params)
        elif model_name == "RFF":
            return RandomForestRegressor(**params)
        elif model_name == "MLP":
            return MLPRegressor(**params, max_iter=MLP_ITER)
        elif model_name == "ElasticNet":
            return ElasticNet(**params)
        elif model_name == "BayesianRidge":
            return BayesianRidge(**params)
        elif model_name == "Huber":
            return HuberRegressor(**params)
        elif model_name == "KNN":
            return KNeighborsRegressor(**params)
        elif model_name == "GBR":
            return GradientBoostingRegressor(**params)
        else:
            raise ValueError(f"Model {model_name} not recognized.")
        
    def suggest_model_params(self, model_type: str, trial: optuna.Trial) -> Dict[str, Any]:
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
                'random_state': 42,  # Aggiunto per riproducibilità
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
                'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
                'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 100.0),
                'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 100.0),
                'max_bin': trial.suggest_categorical('max_bin', [128, 255, 512]),  # Valori ottimali LGB
                'min_child_weight': trial.suggest_float('min_child_weight', 0.001, 50.0),
                'bagging_freq': trial.suggest_int('bagging_freq', 0, 10),
                'feature_fraction': trial.suggest_float('feature_fraction', 0.6, 1.0),  # Aggiunto
                'random_state': 42,
                'verbose': -1,  # Silenzia output
                'n_jobs': -1
            }
        elif model_type == "Cat":
            params ={
                "hidden_layer_sizes": trial.suggest_categorical("hidden_layer_sizes",
                        [(128,), (256,), (384,), (512,),(256,128), (256,256), (384,192),(512,256), (512,256,128)]),
                "activation": trial.suggest_categorical("activation", ["relu", "tanh"]),
                "solver": trial.suggest_categorical("solver", ["adam", "sgd"]),
                "alpha": trial.suggest_float("alpha", 1e-6, 1e-2), # L2 reg più stringente con più dati
                "learning_rate": trial.suggest_categorical("learning_rate", ["constant", "adaptive"]),
                "learning_rate_init": trial.suggest_float("learning_rate_init", 1e-4, 5e-3),
                "batch_size": trial.suggest_categorical("batch_size", [64, 128, 256, "auto"]),
                "beta_1": trial.suggest_float("beta_1", 0.85, 0.99),
                "beta_2": trial.suggest_float("beta_2", 0.95, 0.9999),
                "momentum": trial.suggest_float("momentum", 0.8, 0.99),
                "nesterovs_momentum": trial.suggest_categorical("nesterovs_momentum", [True, False]),
                "early_stopping": True,
                "validation_fraction": 0.1,
                "max_iter": trial.suggest_int("max_iter", 200, 600),
                "n_iter_no_change": trial.suggest_int("n_iter_no_change", 10, 30),
                "tol": trial.suggest_float("tol", 1e-5, 1e-3),
                "random_state": 42,
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
                'random_state': 42
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
                'random_state': 42
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
                'random_state': 42
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
                'random_state': 42
            }       
        elif model_type == "ElasticNet":
            params = {
                'alpha': trial.suggest_float('alpha', 1e-5, 100.0),  # Range esteso
                'l1_ratio': trial.suggest_float('l1_ratio', 0.01, 0.99),  # Evita estremi 0 e 1
                'fit_intercept': trial.suggest_categorical('fit_intercept', [True, False]),  # Aggiunto
                'normalize': trial.suggest_categorical('normalize', [True, False]),  # Aggiunto (deprecato in sklearn >= 1.2)
                'max_iter': trial.suggest_int('max_iter', 1000, 10000),  # Aggiunto
                'selection': trial.suggest_categorical('selection', ['cyclic', 'random']),  # Aggiunto
                'random_state': 42
            }
        elif model_type == "BayesianRidge":
            params = {
                'alpha_1': trial.suggest_float('alpha_1', 1e-8, 1e-3),  # Range più specifico
                'alpha_2': trial.suggest_float('alpha_2', 1e-8, 1e-3),
                'lambda_1': trial.suggest_float('lambda_1', 1e-8, 1e-3),
                'lambda_2': trial.suggest_float('lambda_2', 1e-8, 1e-3),
                'n_iter': trial.suggest_int('n_iter', 100, 1000),  # Aggiunto
                'fit_intercept': trial.suggest_categorical('fit_intercept', [True, False]),  # Aggiunto
                'normalize': trial.suggest_categorical('normalize', [True, False])  # Aggiunto
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
                'random_state': 42
            }
      
        return params


    def optimize(self, objective_func: Callable, n_trials: int = 100, timeout: Optional[int] = None, show_progress: bool = False):
        if not self.study:
            self.create_study()

        lprint(ll.INFO, f"Starting Optimization: {n_trials} trials")

        def logging_callback(study, trial):
            self.optimization_history.append({
                'trial': len(study.trials),
                'value': trial.value,
                'params': trial.params,
                'state': trial.state.name,
                'duration': trial.duration.total_seconds() if trial.duration else None
            })
            if trial.value is not None:
                lprint(ll.REPORT, f"Trial {len(study.trials)}: {trial.value:.4f} | Best: {study.best_value:.4f}")

        try:
            self.study.optimize(
                objective_func,
                n_trials=n_trials,
                timeout=timeout,
                show_progress_bar=show_progress,
                callbacks=[logging_callback]
            )
            if not self.study.trials or all(t.state != optuna.trial.TrialState.COMPLETE for t in self.study.trials):
                lprint(ll.ERROR, "No trials completed successfully")
                raise ValueError("No trials completed successfully")
            self.best_params = self.study.best_params
            self.best_score = self.study.best_value
            lprint(ll.SUCCESS, f"Optimization completed!")
            lprint(ll.REPORT, f"Best MAE: {-self.best_score:.4f}")
            lprint(ll.REPORT, f"Best params: {self.best_params}")
        except Exception as e:
            lprint(ll.ERROR, f"Optimization failed: {str(e)}")
            raise
        return self.study

    def cross_validate_with_pruning(self, model, X, y, cv_folds: int = 5, scoring: str = 'neg_mean_absolute_error'):
        kf = KFold(n_splits=cv_folds, shuffle=True, random_state=SEED)
        try:
            scores = cross_val_score(model, X, y, cv=kf, scoring=scoring, n_jobs=-1)
            return scores.mean(), scores.std()
        except Exception as e:
            lprint(ll.ERROR, f"Cross-validation failed: {str(e)}")
            raise

    def get_optimization_summary(self) -> Dict[str, Any]:
        if not self.study:
            lprint(ll.ERROR, "No study available")
            return {}
        trials_df = self.study.trials_dataframe()
        completed_trials = trials_df[trials_df['state'] == 'COMPLETE']
        summary = {
            'study_name': self.study_name,
            'n_trials': len(self.study.trials),
            'n_completed_trials': len(completed_trials),
            'n_pruned_trials': len(trials_df[trials_df['state'] == 'PRUNED']),
            'best_value': self.study.best_value if completed_trials.shape[0] > 0 else None,
            'best_params': self.study.best_params if completed_trials.shape[0] > 0 else None,
            'best_trial_number': self.study.best_trial.number if completed_trials.shape[0] > 0 else None,
            'optimization_time_total': sum([t.duration.total_seconds() for t in self.study.trials if t.duration]),
            'avg_trial_time': np.mean([t.duration.total_seconds() for t in self.study.trials if t.duration]) if self.study.trials else 0,
            'sampler': self.sampler_type,
            'pruner': self.pruner_type
        }
        return summary

    def save_study(self, filepath: str):
        if not self.study:
            lprint(ll.ERROR, "No study to save")
            return
        study_data = {
            'study': self.study,
            'best_params': self.best_params,
            'best_score': self.best_score,
            'optimization_history': self.optimization_history,
            'summary': self.get_optimization_summary()
        }
        with open(filepath, 'wb') as f:
            pickle.dump(study_data, f)
        lprint(ll.INFO, f"Study saved: {filepath}")

    def load_study(self, filepath: str):
        with open(filepath, 'rb') as f:
            study_data = pickle.load(f)
        self.study = study_data['study']
        self.best_params = study_data['best_params']
        self.best_score = study_data['best_score']
        self.optimization_history = study_data['optimization_history']
        lprint(ll.INFO, f"Study loaded: {filepath}")
        return study_data

    def export_trials_to_csv(self, filepath: str):
        if not self.study:
            lprint(ll.ERROR, "No study available")
            return
        trials_df = self.study.trials_dataframe()
        trials_df.to_csv(filepath, index=False)
        lprint(ll.INFO, f"Trials exported: {filepath}")

def create_regression_objective(
    tuner: OptunaTuner,
    model_type: str,
    X,
    y,
    param_suggest_func: Callable,
    cv_folds: int = 5,
    scoring: str = 'neg_mean_absolute_error'
):
    def objective(trial):
        params = param_suggest_func(tuner, model_type, trial)
        model = tuner.get_model_instance(model_type, params)
        score_mean, score_std = tuner.cross_validate_with_pruning(
            model, X, y, cv_folds, scoring
        )
        mae = -score_mean 
        
        trial.report(mae, 0)
        if trial.should_prune():
            raise optuna.TrialPruned()
        return mae  
    
    return objective