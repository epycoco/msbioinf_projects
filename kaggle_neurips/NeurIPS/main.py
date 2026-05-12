import os, sys
from pathlib import Path
parent_folder = Path(__file__).resolve().parent
if str(parent_folder) not in sys.path:
    sys.path.insert(0, str(parent_folder))

from typing import List, Tuple, Optional
import pandas as pd
import numpy as np
from src.tasks import load_data, merge_data, validate_smiles, compute_df, split_data, train, tanimoto

from src.config.fingerprint import compute_morgan_maccs_fp
from src.config.chemberta import compute_chemberta
from src.config.descriptors import compute_all_desc, filter_desc_corrTarget, filter_desc_pca_corr

from src.logging import lprint, LoggingLevels as ll
from src.config.config import *

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.mixture import GaussianMixture
import plotly.express as px


def feature_embadding (ismi_df: pd.DataFrame, phase: str,
                       bdesc: bool=True, not3d: bool=False, desc_list: Optional[List[str]] = None,
                       be3fp: bool=False, 
                       bberta: bool=True,
                       step: int=1) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
        ismi_df must contain the following columns = ['id', 'SMILES']
    """
    idesc_df = None
    immfp_df = None
    i3fp_df = None
    iberta_df = None
    step_count = 1

    ismi_df['id'] = ismi_df['id'].astype(str)

    if bdesc is True:
        lprint(ll.STEP, f"STEP 2.{step}.{step_count}: Computing Descriptors")
        step_count += 1
        idesc_df = compute_all_desc(ismi_df, not3d, desc_list=desc_list, phase= phase)
        lprint(ll.SUCCESS, "Descriptors Computed.\n")
        idesc_df['id'] = idesc_df["id"].astype(str)        

    lprint(ll.STEP, f"STEP 2.{step}.{step_count}: Computing Morgan Fingerprint + MACCSKey")
    step_count += 1
    immfp_df = compute_morgan_maccs_fp(ismi_df= ismi_df, phase= phase)
    immfp_df['id'] = immfp_df["id"].astype(str)  
    lprint(ll.SUCCESS, "Morgan Fingerprint + MACCSKey Computed.\n")

    if bberta is True:
        lprint(ll.STEP, f"STEP 2.{step}.{step_count}: Computing ChemBERTa Embeddings")
        step_count += 1
        iberta_df = compute_chemberta(ismi_df, phase= phase)
        iberta_df['id'] = iberta_df["id"].astype(str)  
        lprint(ll.SUCCESS, "ChemBERTa Embeddings Computed.")
    
    if be3fp is True:
        lprint(ll.STEP, f"STEP 2.{step}.{step_count}: Computing 3D Fingerprint")
        ## from learn_fingerprint
        lprint(ll.SUCCESS, "3D Fingerprint Computed")
    
    all_df = pd.merge(ismi_df, immfp_df, on='id')
    if idesc_df is not None:
        all_df = pd.merge(all_df, idesc_df, on='id')

    if iberta_df is not None:
        all_df = pd.merge(all_df, iberta_df, on='id')

    if i3fp_df is not None:
        all_df = pd.merge(all_df, i3fp_df, on='id')
    
    if all_df.shape == ismi_df.shape:
        lprint(ll.ERROR, f'There is no valid feature embaddings')
        raise ValueError

    return all_df, idesc_df, immfp_df, i3fp_df, iberta_df 


def best_gmm_ncomponents(X, max_components=50, n_splits=5, random_state=42, verbose=False):
    scores = {}
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    scount = 0
    for n in range(1, max_components + 1):
        fold_scores = []
        for train_idx, test_idx in kf.split(X):
            gmm = GaussianMixture(n_components=n, random_state=random_state)
            gmm.fit(X[train_idx])
            fold_scores.append(gmm.score(X[test_idx]))
        mean_score = np.mean(fold_scores)
        scores[n] = mean_score
        if mean_score < max(scores, key=scores.get):
            scount += 1
            if scount > 5:
                break

    best_n = max(scores, key=scores.get)
    return best_n


def feature_inference (inan_df: pd.DataFrame, ifp_df: pd.DataFrame):
    for i, desc in enumerate(inan_df.columns.to_list(), 1):
        if inan_df[desc].hasnans:
            lprint(ll.INFO, f"[{i}] Processing descriptor: {desc}")

            valid_df = inan_df[['id', desc]].dropna()
            valid_ifp_df = ifp_df.loc[valid_df['id'].to_list(), :]
            X_valid = valid_ifp_df.drop(columns='id').to_numpy(dtype=np.int8)
            y_valid = valid_df[desc].to_numpy(dtype=np.float32)

            bestk = best_gmm_ncomponents(X_valid)

            # addestramento GMM sulle fingerprint valide
            gmm = GaussianMixture(n_components=bestk, random_state=42)
            gmm.fit(X_valid)

            # cluster assignment e media descrittore per cluster
            cluster_labels = gmm.predict(X_valid)
            cluster_means = {
                c: y_valid[cluster_labels == c].mean()
                for c in np.unique(cluster_labels)
            }

            nan_df  = inan_df[inan_df[desc].isna()]
            nan_ifp_df = ifp_df[ifp_df['id'].isin(nan_df['id'])]
            X_nan = nan_ifp_df.drop(columns='id').to_numpy(dtype=np.int8)

            probs = gmm.predict_proba(X_nan)

            # imputazione: media ponderata delle medie dei cluster
            inferred_values = []
            for p in probs:
                weighted_mean = np.sum([p[c] * cluster_means[c] for c in range(bestk)])
                inferred_values.append(weighted_mean)

            # aggiorna il DataFrame originale
            inan_df.loc[inan_df[desc].isna(), desc] = inferred_values

    return inan_df
                

# Find Outliers to remove in target column in order to remove noise 
def remove_outliers_iqr(df: pd.DataFrame, target_col: str, id_col: str = "id") -> pd.DataFrame:
    """
        Return a pd.DataFrame without the outliers
    """
    Q1 = df[target_col].quantile(0.25)
    Q3 = df[target_col].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    df[id_col] = df[id_col].astype(str)

    # Trova outlier
    outliers = df[(df[target_col] < lower_bound) | (df[target_col] > upper_bound)]
    outlier_ids = outliers[id_col].tolist()

    # Rimuovi outlier
    df_clean = df[~df[id_col].isin(outlier_ids)].reset_index(drop=True)

    return df_clean




def main():
    #--------------#
    # Data loading #
    #--------------#
    if not os.path.exists(TRAIN_PATH): 
        lprint(ll.STEP, "STEP 0.1: Loading data")
        try:
            (train_df, test_df), ext_dfs, data_extend_dict = load_data.main()
            lprint(ll.INFO, f'Train shape: {train_df.shape}, Test shape: {test_df.shape}')
            lprint(ll.INFO, f"External datasets: {len(ext_dfs)}")
        except Exception as e:
            lprint(ll.ERROR, f"Failed to load data: {str(e)}")
            raise
        lprint(ll.SUCCESS, "STEP 0.1: Completed")
        
        #------------#
        # Data Merge #
        #------------#
        lprint(ll.STEP, "STEP 0.2: Merging official and support data")
        try:
            train_df = merge_data.main(train_df, ext_dfs, data_extend_dict)
            train_df.to_csv(TRAIN_PATH, index=0)
            lprint(ll.INFO, f"Merged train shape: {train_df.shape}")
        except Exception as e:
            lprint(ll.ERROR, f"Failed to merge data: {str(e)}")
            raise
        lprint(ll.SUCCESS, "STEP 0.2: Completed")
    
    else:
        #--------------#--------------#
        # Already Merged Data Loading #
        #--------------#--------------#
        train_df = pd.read_csv(TRAIN_PATH, header=0)
        test_df = pd.read_csv(TEST_PATH, header=0)


    #----------------#---------------#
    # SMILES Validation and Cleaning #
    #----------------#---------------#
    lprint(ll.STEP, "STEP 1: SMILES validation and cleaning")
    if not os.path.exists(TRAIN_CLEAN_PATH) or not os.path.exists(TEST_CLEAN_PATH):
        try:
            train_df, test_df, not3d = validate_smiles.main(train_df, test_df)
            lprint(ll.INFO, f"Validated Train shape: {train_df.shape}, Test shape: {test_df.shape}")
        except Exception as e:
            lprint(ll.ERROR, f"Failed to validate SMILES: {str(e)}")

        train_df.to_csv(TRAIN_CLEAN_PATH, index=False)
        lprint(ll.INFO, "Clean Train dataset file: SAVED")
        test_df.to_csv(TEST_CLEAN_PATH, index=False)
        lprint(ll.INFO, "Clean Test dataset file: SAVED")
    
    else:
        train_df = pd.read_csv(TRAIN_CLEAN_PATH, header=0)
        lprint(ll.INFO, "A clean Train dataset file was found: UPLOADED")
        test_df = pd.read_csv(TEST_CLEAN_PATH, header=0)
        lprint(ll.INFO, "A clean Test dataset file was found: UPLOADED")

        not3d = False
        test_df['SMILES'] = test_df['SMILES'].apply(validate_smiles.validate_smiles_3d)
        for smi in test_df['SMILES']:
            if smi is None:
                lprint(ll.WARN, f"It is impossibile to compute 3D descriptors because of one SMILES in test file.")
                not3d = True
                break
    lprint(ll.SUCCESS, "STEP 1: Completed\n")


    #-------------------#--------------#
    # Feature Embedding for Train/Test #
    #-------------------#--------------#
    lprint(ll.STEP, "STEP 2.1: Feature Embedding for Train Set")
    train_all_df, idesc_df, immfp_df, i3fp_df, iberta_df = feature_embadding(ismi_df=train_df, phase='train', bdesc=True, not3d=not3d, bberta=True, step=1)
    if idesc_df is not None:
        lprint(ll.REPORT, f'Full idesc_df set shape: {idesc_df.shape}')
    if immfp_df is not None:
        lprint(ll.REPORT, f'Full immfp_df set shape: {immfp_df.shape}')
    if i3fp_df is not None:
        lprint(ll.REPORT, f'Full immfp_df set shape: {i3fp_df.shape}')
    if iberta_df is not None:
        lprint(ll.REPORT, f'Full iberta_df set shape: {iberta_df.shape}')
    lprint(ll.REPORT, f'Full Train set shape: {train_all_df.shape}')
    lprint(ll.SUCCESS, "STEP 2.1: Completed.\n")

    lprint(ll.STEP, "STEP 2.2: Feature Embedding for Test Set")
    desc_list = None
    if idesc_df is not None:
        desc_list = idesc_df.drop(columns='id').columns.to_list()
    test_all_df, idesc_test_df, immfp_test_df, i3fp_test_df, iberta_test_df = feature_embadding(ismi_df=test_df, phase='test', bdesc=True, not3d=not3d, desc_list=desc_list, bberta=True, step=2)
    lprint(ll.SUCCESS, "STEP 2.2: Completed.\n")


    #-------------------#--------------#
    # SMILES Similarity by Fingerprint #
    #-------------------#--------------#
    lprint(ll.STEP, "STEP 3: Compute SMILES Similarity on Fingerprint between Train and Test")
    mmfp_df = immfp_df.drop(columns='id')
    mmfp_test_df = immfp_test_df.drop(columns='id')
    tanimoto_sim = tanimoto.tanimoto(mmfp_df, mmfp_test_df)
    lprint(ll.SUCCESS, f'Tanimoto similarity for each molecule on test set: {tanimoto_sim}\n')


    #----------------------#
    # TRAIN AND VALIDATION #
    #----------------------#
    ftt_df = None
    lprint(ll.STEP, f'STEP 4: Training and Validation')
    for i, target in enumerate(TARGETS):
        #-------------------------#
        # EDA ON TARGET DATAFRAME #
        #-------------------------#
        lprint(ll.STEP, f'STEP 4.{i}.0: EDA on {target} subset')

        mmfp_feature = immfp_df.drop(columns='id').columns.to_list()
        target_features = ['id', target] + mmfp_feature

        if iberta_df is not None:
            berta_feature = iberta_df.drop(columns='id').columns.to_list()
            target_features += berta_feature

        if i3fp_df is not None:
            e3fp_feature = i3fp_df
            target_features += e3fp_feature

        if idesc_df is not None:
            t_feat_no_desc = target_features.copy()
            desc_feature = idesc_df.drop(columns=['id']).columns.to_list()
            target_features += desc_feature

            target_train_df = train_all_df[target_features]
            target_train_df.dropna(subset= target, inplace=True)  # create the subset of target

            lprint(ll.REPORT, f'{target} Train Subset shape: {train_all_df.shape}')

            lprint(ll.STEP, f"STEP 4.{i}.0.a: Remove Outliers")
            target_train_df = remove_outliers_iqr(df= target_train_df, target_col= target) # remove outliers

            lprint(ll.STEP, f"STEP 4.{i}.0.b: Feature Inference on Descriptors")
            target_train_df[['id', *desc_feature]] = feature_inference(inan_df= target_train_df[['id', *desc_feature]], ifp_df= target_train_df[['id', *mmfp_feature]])

            lprint(ll.STEP, f"STEP 4.{i}.0.c: Filtering Descriptors by their Correlation with {target} values")
            desc_filter_df = filter_desc_corrTarget(desc_df= target_train_df[desc_feature], target= target_train_df[target], t= target) # remove descriptor not correlated
            
            
            # if there are not correlated descriptors with target
            if desc_filter_df is None:
                ftt_df = target_train_df[t_feat_no_desc] # Final Target Train DataFrame
                lprint(ll.REPORT, f"Final Target Train Subset shape: {ftt_df.shape}")

            # if there are descriptors correlated with target -> PCA and correlation between descriptors
            else:
                desc_feature = desc_filter_df.columns.to_list()
                
                # feature selection with PCA and Correlation 
                lprint(ll.STEP, f"STEP 4.{i}.0.d: Filtering Descriptors by PCA and Correlation Matrix")
                desc_filter_df = filter_desc_pca_corr(desc_df= target_train_df[desc_feature], target= target)

                if desc_filter_df is None:
                    ftt_df = target_train_df[t_feat_no_desc]
                
                else:
                    desc_feature = desc_filter_df.columns.to_list()
                    target_features = t_feat_no_desc + desc_feature
                    ftt_df = target_train_df[target_features]
        
        else:
            target_train_df = train_all_df[target_features].dropna(subset= target)
            target_train_df = remove_outliers_iqr(df= target_train_df, target_col= target) # remove outliers
            ftt_df = target_train_df # because in this 'else' condition 'idesc_df' is None, and there are not modifications in 'target_train_df' columns (drop only rows)

        lprint(ll.REPORT, f"Final Target Train Subset shape: {ftt_df.shape}")
        lprint(ll.STEP, f'STEP 4.{i}.0: Completed\n')

        lprint(ll.STEP, f'STEP 4.{i}.1: Training on {target} with MAE as Loss Function')
        if ftt_df is not None:
            target_train_dir = f'./data/target_train/'
            target_train_file = f'{target_train_dir}/{target}_train.csv'
            if not os.path.exists(target_train_dir):
                os.makedirs(target_train_dir)
            ftt_df.to_csv(target_train_file, index=False)
        else:
            raise ValueError

        X = ftt_df.drop(columns=['id', target])
        y = ftt_df[target]     

        info_val_loss = {
                        #a_idx, mae
            #'rmse':    [ None, 10000],
            'MAE':      [ None, 10000],
            #'alfa_m':  [ None, 10000],
            #'alfa_l':  [ None, 10000],
            #'alfa_u':  [ None, 10000]
        }

        # MAE
        lprint(ll.STEP, f'STEP 4.{i}.1: Build Ensamble Model for {target}')
        ens_mae = train.train(X_train=X, y_train=y, X_test=None, target=target, loss_mode="MAE")
        lprint (ll.WARN, ens_mae)
        if ens_mae["mae"] < info_val_loss["MAE"][1]:
            info_val_loss["MAE"] = [ens_mae, ens_mae["type"], ens_mae["mae"]]
        lprint(ll.SUCCESS, f'STEP 4.{i}.1: Completed\n')


        """#RMSE
        lprint(ll.STEP, f'STEP 3.{i}.2: Training on {target} with MSE as Loss Function')
        ens_mse = train.train(X_train=X, y_train=y, X_test=None, target=target, loss_mode="MSE")
        rmse = ens_mse["val_loss"] ** 0.5
        if rmse < info_val_loss["rmse"][1]:
            info_val_loss["rmse"] = [ens_mse, rmse, ens_mse['mae']]
        lprint(ll.SUCCESS, f'STEP 3.{i}.2: Completed\n')
        
        # CUSTOM alpha=0.5
        lprint(ll.STEP, f'STEP 3.{i}.3: Training on {target} with Custom Loss Function (alpha = 0.5)')
        ens_custom = train.train(X_train=X, y_train=y, X_test=None, target=target, loss_mode="CUSTOM", alpha=0.5)
        if ens_custom["val_loss"] < info_val_loss["alfa_m"][1]:
            info_val_loss["alfa_m"] = [ens_custom, ens_custom["val_loss"]]
        lprint(ll.SUCCESS, f'STEP 3.{i}.3: Completed\n')

        for i, alfa in enumerate(ALFA_LOW):
            # --- ALFA LOW ---
            if not stop_low:
                lprint(ll.STEP, f'STEP 3.{i}.{i+3}.L: Training on {target} with Custom Loss Function (alpha = {alfa})')
                ens_custom_low = train.train(X_train=X, y_train=y, X_test=None, target=target, loss_mode="CUSTOM", alpha=alfa)

                if ens_custom_low["val_loss"] < info_val_loss["alfa_l"][1]:
                    info_val_loss["alfa_l"] = [ens_custom_low, ens_custom_low["val_loss"]]
                    no_improve_low = 0  # reset contatore se migliora
                else:
                    no_improve_low += 1
                    if no_improve_low >= 2:
                        stop_low = True
                        lprint(ll.WARN, f"Early stopping ALFA_LOW at index {i}, last val_loss={ens_custom_low['val_loss']:.4f}")
                lprint(ll.SUCCESS, f'STEP 3.{i}.{i+3}.L: Completed\n')


            # --- ALFA UP ---
            if not stop_up:
                lprint(ll.STEP, f'STEP 3.{i}.{i+3}.U: Training on {target} with Custom Loss Function (alpha = {ALFA_UP[i]})')
                ens_custom_up = train.train(X_train=X, y_train=y, X_test=None, target=target, loss_mode="CUSTOM", alpha=ALFA_UP[i])

                if ens_custom_up["val_loss"] < info_val_loss["alfa_u"][1]:
                    info_val_loss["alfa_u"] = [ens_custom_up, ens_custom_up["val_loss"]]
                    no_improve_up = 0
                else:
                    no_improve_up += 1
                    if no_improve_up >= 2:
                        stop_up = True
                        lprint(ll.WARN, f"Early stopping ALFA_UP at index {i}, last val_loss={ens_custom_up['val_loss']:.4f}")

                lprint(ll.SUCCESS, f'STEP 3.{i}.{i+3}.U: Completed\n')

            # if both flag are False, break
            if stop_low and stop_up:
                break"""


if __name__ == "__main__":
    main()