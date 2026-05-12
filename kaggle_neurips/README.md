# NeurIPS Open Polymer Prediction 2025

A machine learning pipeline for predicting polymer physicochemical properties from SMILES representations, developed for the [NeurIPS 2025 Open Polymer Prediction Challenge](https://www.kaggle.com/competitions/neurips-open-polymer-prediction-2025).

## Overview

This project predicts five key polymer properties from molecular SMILES strings:

| Target | Description |
|--------|-------------|
| `Tg` | Glass Transition Temperature (K) |
| `FFV` | Fractional Free Volume |
| `Tc` | Crystallization Temperature (K) |
| `Density` | Polymer Density (g/cm³) |
| `Rg` | Radius of Gyration (Å) |

The pipeline combines multiple molecular featurization strategies with an ensemble of gradient-boosted tree models, tuned via Optuna.

## Project Structure

```
NeurIPS/
├── main.py                    # Main training pipeline
├── ensemble_outs.py           # Ensemble model definition and training
├── build_desc_df.py           # Descriptor DataFrame builder
├── learn_fingerprint.py       # Fingerprint learning utilities
├── prova_predict.py           # Inference / prediction script
├── src/
│   ├── config/
│   │   ├── config.py          # Global constants and paths
│   │   ├── models.py          # Per-target model hyperparameters
│   │   ├── chemberta.py       # ChemBERTa embedding computation
│   │   ├── descriptors.py     # Molecular descriptor computation and filtering
│   │   └── fingerprint.py     # Morgan + MACCS fingerprint computation
│   ├── tasks/
│   │   ├── augmentation.py    # GMM-based data augmentation
│   │   ├── berta.py           # ChemBERTa fine-tuning utilities
│   │   ├── compute_df.py      # Feature DataFrame assembly
│   │   ├── hashing.py         # SMILES hashing and deduplication
│   │   ├── io.py              # Model/config save & load
│   │   ├── load_data.py       # Data loading (official + supplementary)
│   │   ├── merge_data.py      # Merging official and external datasets
│   │   ├── predict.py         # Prediction utilities
│   │   ├── split_data.py      # Train/validation splitting
│   │   ├── tanimoto.py        # Tanimoto similarity between train and test
│   │   ├── train.py           # Training loop with Optuna tuning
│   │   └── validate_smiles.py # SMILES validation and 3D readiness check
│   └── logging.py             # Structured logging utilities
├── data/
│   ├── neurips-open-polymer-prediction-2025/   # Official competition data
│   │   ├── train.csv
│   │   ├── test.csv
│   │   └── sample_submission.csv
│   ├── train_supports/        # External/supplementary datasets
│   ├── best_params/           # Saved Optuna best configs (.pkl)
│   └── fitting_results.txt    # Per-model MAE and R² on training data
├── ChemBERTa-77M-MTR/         # Pre-trained ChemBERTa model (HuggingFace)
└── best_params/               # Best hyperparameter configs per target & fold
```

## Pipeline

The pipeline runs in sequential steps:

1. **Data Loading & Merging** — loads the competition dataset and merges it with supplementary polymer databases (`train_supports/`)
2. **SMILES Validation** — validates and normalizes SMILES strings; checks 3D conformer feasibility
3. **Feature Embedding** — computes three feature types for each molecule:
   - **Morgan + MACCS fingerprints** (256-bit Morgan radius-2 + 167-bit MACCS keys)
   - **Molecular descriptors** via RDKit and Mordred (2D/3D, filtered by variance, NaN rate, and redundancy)
   - **ChemBERTa embeddings** (384-dim, from `seyonec/ChemBERTa-zinc-base-v1` fine-tuned on MTR)
4. **EDA per target**:
   - Outlier removal (IQR method)
   - GMM-based feature inference for missing descriptor values
   - Descriptor filtering by correlation with target and PCA/inter-descriptor correlation
5. **Training & Tuning** — per target, trains an ensemble of models with Optuna hyperparameter optimization (K=5 cross-validation, MAE loss):
   - XGBoost, RandomForest, ExtraTrees, HistGradientBoosting, GradientBoosting
6. **Ensemble** — combines models with optimized stacking weights

## Installation

```bash
pip install -r requirements.txt
```

Key dependencies:
- `rdkit`, `mordred` — molecular featurization
- `transformers`, `torch` — ChemBERTa embeddings
- `xgboost`, `lightgbm`, `catboost`, `scikit-learn` — model training
- `optuna` — hyperparameter tuning
- `pandas`, `numpy`, `plotly`

## Usage

### Full training pipeline

```bash
python main.py
```

The script caches intermediate results (merged data, clean SMILES, features) so subsequent runs skip already-computed steps.

### Prediction on new data

```bash
python prova_predict.py
```

Expects a CSV with columns `id` and `SMILES`. Outputs predictions for all five targets.

### Hyperparameter tuning only

```bash
python ensemble_outs.py
```

Runs Optuna tuning per target and saves best configs to `best_params/`.

## Configuration

Key parameters in `src/config/config.py`:

```python
TARGETS       = ['Tg', 'FFV', 'Tc', 'Density', 'Rg']
MODELS        = ['XGB', 'RFF', 'ET', 'HGB', 'GBR']
K_FOLD        = 5
OPTIM_TRIALS  = 20
FP_RADIUS     = 2
FP_BITS       = 256
BERTA_BITS    = 384
SEED          = 42
DEVICE        = "cuda"
```

## Results

Best validation performance (MAE / R²) on training data:

| Target | Best Model | MAE | R² |
|--------|-----------|-----|----|
| Tg | RandomForest | 55.59 | 0.648 |
| FFV | HistGradientBoosting | 0.00075 | 0.9979 |
| Tc | HistGradientBoosting | 0.0198 | 0.883 |
| Density | HistGradientBoosting | 0.00096 | 0.9999 |
| Rg | HistGradientBoosting | 0.028 | 0.9999 |

## External Datasets Used

The following supplementary datasets are merged during training:

- `PI1070.csv` — polyimide dataset
- `Tg_SMILES_class_pid_polyinfo_median.csv` — glass transition temperatures
- `TgSS_enriched_cleaned.csv` — enriched Tg dataset
- `Tc_SMILES.csv` — crystallization temperatures
- `JCIM_sup_bigsmiles.csv` — BigSMILES polymer dataset
- `extended_polymer_dataset.csv` — general polymer properties
- `data_dnst1.xlsx`, `data_tg3.xlsx` — density and Tg supplementary data

## Notes

- The `ChemBERTa-77M-MTR/` directory must be present locally (tracked via Git LFS).
- 3D descriptor computation is automatically skipped if any test SMILES cannot generate a 3D conformer.
- Feature caching: processed CSVs (descriptors, fingerprints, embeddings) are saved to `data/` and reused on subsequent runs to save computation time.
