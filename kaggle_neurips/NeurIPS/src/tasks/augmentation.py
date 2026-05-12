import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, MACCSkeys, AllChem, rdmolops
from rdkit.Chem.rdMolDescriptors import CalcTPSA, CalcNumRotatableBonds
from rdkit.Chem.Descriptors import MolWt, MolLogP
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator, GetAtomPairGenerator, GetTopologicalTorsionGenerator
import networkx as nx
from sklearn.mixture import GaussianMixture
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from src.logging import lprint, LoggingLevels as ll
from src.config.config import *

import os

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from datetime import datetime
import os

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

import matplotlib.pyplot as plt
import os

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from src.logging import lprint, LoggingLevels as ll
from src.config.config import SEED

def evaluate_gmm_components(X, target, max_components=20, random_state=SEED, save_path="gmm_evaluation"):
    """
    Evaluates the optimal number of components for the Gaussian Mixture Model using BIC and AIC.
    Standardizes input data, includes checks for data validity, and evaluates covariance types.

    Parameters:
        X (pd.DataFrame or np.ndarray): Input data for the GMM (standardized).
        target (str): Name of the target for identifying the plot.
        max_components (int): Maximum number of components to test.
        random_state (int): Random seed for reproducibility.
        save_path (str): Directory to save the plots (created if it doesn't exist).

    Returns:
        dict: Dictionary containing BIC and AIC scores for each number of components and covariance types.
    """
    os.makedirs(save_path, exist_ok=True)
    lprint(ll.INFO, f"Save directory: {save_path}")

    # Check validity of input data
    if isinstance(X, pd.DataFrame):
        X = X.values  # Convert DataFrame to NumPy array
    if not isinstance(X, np.ndarray):
        lprint(ll.ERROR, f"Input X must be a NumPy array or DataFrame, received: {type(X)}")
        raise ValueError("Invalid input X")
    if np.any(np.isnan(X)) or np.any(np.isinf(X)):
        lprint(ll.ERROR, f"Input X contains NaN or infinite values for target {target}")
        raise ValueError("Invalid data: NaN or infinite values found")
    lprint(ll.INFO, f"Dataset size for {target}: {X.shape}")

    # Evaluate GMM for components 1 to max_components (original functionality)
    bic_scores = []
    aic_scores = []
    log_likelihoods = []
    components_range = list(range(1, max_components + 1))

    for n_components in components_range:
        try:
            gmm = GaussianMixture(n_components=n_components, random_state=random_state)
            gmm.fit(X)
            bic = gmm.bic(X)
            aic = gmm.aic(X)
            log_likelihood = gmm.score(X) * X.shape[0]  # Total log-likelihood
            bic_scores.append(bic)
            aic_scores.append(aic)
            log_likelihoods.append(log_likelihood)
            lprint(ll.DEBUG, f"Components: {n_components}, BIC: {bic:.2f}, AIC: {aic:.2f}, Log-likelihood: {log_likelihood:.2f}")
        except Exception as e:
            lprint(ll.ERROR, f"Error during GMM fit for {n_components} components: {str(e)}")
            raise

    # Calculate elbow point for BIC
    def find_elbow(scores):
        if len(scores) < 3:
            return components_range[np.argmin(scores)]
        relative_changes = [(scores[i-1] - scores[i]) / abs(scores[i-1]) for i in range(1, len(scores))]
        threshold = 0.05  # Threshold for relative change
        for i, change in enumerate(relative_changes):
            if change < threshold:
                return components_range[i+1]
        return components_range[np.argmin(scores)]

    optimal_components_bic = components_range[np.argmin(bic_scores)]
    optimal_components_aic = components_range[np.argmin(aic_scores)]
    elbow_components = find_elbow(bic_scores)

    # Plot original BIC/AIC results
    plt.figure(figsize=(10, 6))
    plt.plot(components_range, bic_scores, label='BIC', marker='o')
    plt.plot(components_range, aic_scores, label='AIC', marker='s')
    plt.axvline(x=elbow_components, color='r', linestyle='--', label='Elbow point')
    plt.xlabel('Number of Components')
    plt.ylabel('Score')
    plt.title(f'Evaluation of GMM Components - Target: {target} (Standardized Data)')
    plt.legend()
    plt.grid(True)

    filename = f"{save_path}/gmm_components_{target}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    lprint(ll.INFO, f"Plot saved to: {filename}")

    lprint(ll.REPORT, f"Optimal number of components (BIC): {optimal_components_bic} (BIC: {min(bic_scores):.2f})")
    lprint(ll.REPORT, f"Optimal number of components (AIC): {optimal_components_aic} (AIC: {min(aic_scores):.2f})")
    lprint(ll.REPORT, f"Elbow point components: {elbow_components}")

    # Warning if BIC is negative
    if min(bic_scores) < 0:
        lprint(ll.WARN, f"Negative BIC detected ({min(bic_scores):.2f}). Verify data standardization and model validity.")

    # Evaluate GMM for covariance types (1 to 6 components)
    covariance_types = ['full', 'tied', 'diag', 'spherical']
    bic_scores_cov = {cov_type: [] for cov_type in covariance_types}
    aic_scores_cov = {cov_type: [] for cov_type in covariance_types}
    components_range_cov = list(range(1, 7))

    lprint(ll.INFO, f"Evaluating covariance types for target: {target} (components 1 to 6)")
    best_bic = float('inf')
    best_cov_type = None
    best_n_components = None

    for cov_type in covariance_types:
        for n_components in components_range_cov:
            try:
                gmm = GaussianMixture(n_components=n_components, covariance_type=cov_type, random_state=random_state)
                gmm.fit(X)
                bic = gmm.bic(X)
                aic = gmm.aic(X)
                bic_scores_cov[cov_type].append(bic)
                aic_scores_cov[cov_type].append(aic)
                lprint(ll.DEBUG, f"Covariance: {cov_type}, Components: {n_components}, BIC: {bic:.2f}, AIC: {aic:.2f}")

                # Track the best model
                if bic < best_bic:
                    best_bic = bic
                    best_cov_type = cov_type
                    best_n_components = n_components
            except Exception as e:
                lprint(ll.ERROR, f"Error during GMM fit for covariance {cov_type}, {n_components} components: {str(e)}")
                raise

    # Plot covariance type results
    plt.figure(figsize=(12, 8))
    for cov_type in covariance_types:
        plt.plot(components_range_cov, bic_scores_cov[cov_type], label=f'{cov_type} (BIC)', marker='o')
    plt.axvline(x=best_n_components, color='r', linestyle='--', label=f'Best Model: {best_cov_type}, {best_n_components} components (BIC: {best_bic:.2f})')
    plt.xlabel('Number of Components')
    plt.ylabel('BIC Score')
    plt.title(f'GMM Covariance Type Comparison - Target: {target} (Standardized Data)')
    plt.legend()
    plt.grid(True)

    cov_filename = f"{save_path}/gmm_covariance_types_{target}.png"
    plt.savefig(cov_filename, dpi=300, bbox_inches='tight')
    plt.close()
    lprint(ll.INFO, f"Covariance type plot saved to: {cov_filename}")

    lprint(ll.REPORT, f"Best model: {best_cov_type} covariance, {best_n_components} components (BIC: {best_bic:.2f})")

    return {
        'components': components_range,
        'bic_scores': bic_scores,
        'aic_scores': aic_scores,
        'log_likelihoods': log_likelihoods,
        'optimal_components_bic': optimal_components_bic,
        'optimal_components_aic': optimal_components_aic,
        'elbow_components': elbow_components,
        'covariance_types': covariance_types,
        'components_cov': components_range_cov,
        'bic_scores_cov': bic_scores_cov,
        'aic_scores_cov': aic_scores_cov,
        'best_cov_type': best_cov_type,
        'best_n_components_cov': best_n_components
    }


def augment_smiles(smiles_list, targets, num_augments):
    augmented_smiles = []
    augmented_targets = []
    for smiles, target in zip(smiles_list, targets):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        augmented_smiles.append(smiles)
        augmented_targets.append(target)
        for _ in range(num_augments):
            rand_smiles = Chem.MolToSmiles(mol, doRandom=True)
            augmented_smiles.append(rand_smiles)
            augmented_targets.append(target)
    return augmented_smiles, np.array(augmented_targets)


def augment_descriptors(X, y, target, n_samples, n_components, scaler=None, random_state=SEED, cov_type='full', descale_fp=False, descale_target=False):
    """
    Augments the dataset by generating synthetic samples using a Gaussian Mixture Model.
    Assumes X contains standardized non-FP features and binary FP features (0/1), and y is not scaled.
    Non-FP features remain standardized. FP features and target can optionally be inverse-transformed manually.

    Parameters:
        X (pd.DataFrame or np.ndarray): Feature matrix (standardized non-FP features and binary FP features).
        y (pd.Series or np.ndarray): Target values (not scaled).
        target (str): Name of the target column.
        n_samples (int): Number of synthetic samples to generate.
        n_components (int): Number of GMM components.
        scaler (StandardScaler, optional): Pre-fitted scaler for non-FP features (and optionally FP/target).
        random_state (int): Random seed for reproducibility.
        cov_type (str): Covariance type for GMM ('full', 'tied', 'diag', 'spherical').
        descale_fp (bool): If True, apply manual inverse transform to synthetic FP features (not recommended for binary data).
        descale_target (bool): If True, apply manual inverse transform to synthetic target (not recommended if y is not scaled).

    Returns:
        pd.DataFrame: Augmented DataFrame containing original and synthetic data (non-FP and FP features).
    """

    if isinstance(X, np.ndarray):
        X = pd.DataFrame(X)
    elif not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a pandas DataFrame or a NumPy array")

    X.columns = X.columns.astype(str)

    if isinstance(y, np.ndarray):
        y = pd.Series(y)
    elif not isinstance(y, pd.Series):
        raise ValueError("y must be a pandas Series or a NumPy array")

    if np.any(np.isnan(X)) or np.any(np.isinf(X)):
        lprint(ll.ERROR, f"Input X contains NaN or infinite values for target {target}")
        raise ValueError("Invalid data: NaN or infinite values found")
    if np.any(np.isnan(y)) or np.any(np.isinf(y)):
        lprint(ll.ERROR, f"Input y contains NaN or infinite values for target {target}")
        raise ValueError("Invalid data: NaN or infinite values found in y")

    desc_columns = [col for col in X.columns.to_list() if not col.startswith('FP')]
    fp_columns = [col for col in X.columns.to_list() if col.startswith('FP')]
    lprint(ll.INFO, f"Descriptor features for {target}: {len(desc_columns)} columns")
    lprint(ll.INFO, f"FP features for {target}: {len(fp_columns)} columns")

    df = pd.DataFrame(X)
    df[target] = y.values

    gmm_input = df[desc_columns + [target]]
    gmm = GaussianMixture(n_components=n_components, covariance_type=cov_type, random_state=random_state)
    gmm.fit(gmm_input)
    synthetic_data, _ = gmm.sample(n_samples)
    synthetic_df = pd.DataFrame(synthetic_data, columns=desc_columns + [target])

    synthetic_df[fp_columns] = df[fp_columns].iloc[:n_samples].values

    if scaler is not None:
        columns_to_descale = []
        if descale_fp:
            columns_to_descale.extend(fp_columns)
        if descale_target:
            columns_to_descale.append(target)

        if columns_to_descale:
            if not hasattr(scaler, 'mean_') or not hasattr(scaler, 'scale_'):
                lprint(ll.ERROR, f"Scaler trained for {target}")
                raise 
            try:
                feature_names = scaler.feature_names_in_ if hasattr(scaler, 'feature_names_in_') else desc_columns
                feature_to_index = {name: idx for idx, name in enumerate(feature_names)}
            except AttributeError as e:
                lprint(ll.WARN, f"{e}")
                feature_to_index = {name: idx for idx, name in enumerate(desc_columns)}

            synthetic_fp_target = synthetic_df[columns_to_descale].copy()

            for col in columns_to_descale:
                if col in feature_to_index:
                    idx = feature_to_index[col]
                    mean = scaler.mean_[idx]
                    scale = scaler.scale_[idx]
                    synthetic_fp_target[col] = synthetic_fp_target[col] * scale + mean
                    lprint(ll.DEBUG, f"Descaled {col} with mean={mean:.4f}, scale={scale:.4f}")
                else:
                    lprint(ll.WARN, f"Feature {col} not found in scaler, skipping scaling")
                    continue

            if descale_fp:
                synthetic_fp_target[fp_columns] = (synthetic_fp_target[fp_columns] > 0.5).astype(int)

            synthetic_df[columns_to_descale] = synthetic_fp_target

    augmented_df = pd.concat([df, synthetic_df], ignore_index=True)
    lprint(ll.INFO, f"Generated {n_samples} synthetic samples for {target} with {cov_type} covariance")
    lprint(ll.INFO, f"Augmented DataFrame shape for {target}: {augmented_df.shape}")

    return augmented_df.sort_index(axis=1)