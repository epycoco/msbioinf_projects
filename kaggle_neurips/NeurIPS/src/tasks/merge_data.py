import pandas as pd
import numpy as np

from src.logging import lprint, LoggingLevels as ll

def update_or_add_smiles_entries(main_df, ext_df, target_column, 
                                 add_new_smiles=True, fill_value=False):
    """
    Merge polymer datasets with robust handling of new SMILES addition and missing value filling.
    
    Parameters:
    - main_df: Main DataFrame (e.g., train.csv)
    - ext_df: External DataFrame
    - target_column: Column to merge (Tg, FFV, Tc, Density, Rg)
    - add_new_smiles: Whether to add new SMILES with valid target values
    - fill_value: Whether to fill missing values in main_df using ext_df
    
    Returns:
    - Updated DataFrame
    """
    
    # Valid property columns
    valid_columns = {'Tg', 'FFV', 'Tc', 'Density', 'Rg'}
    
    # Validate target column
    if target_column not in valid_columns:
        raise ValueError(f"Target column must be one of: {', '.join(valid_columns)}")
    
    # Check required columns
    required_columns = {'SMILES', target_column}
    if not required_columns.issubset(ext_df.columns):
        missing = required_columns - set(ext_df.columns)
        raise ValueError(f"External dataset missing columns: {', '.join(missing)}")
    
    if 'SMILES' not in main_df.columns:
        lprint(ll.ERROR, "Main dataset must contain 'SMILES' column")
        raise
    
    # Initialize counters
    filled_count = 0
    new_smiles_count = 0
    
    # === [1] Fill missing values if requested ===
    if fill_value:
        merged = main_df.merge(
            ext_df[['SMILES', target_column]].dropna(subset=[target_column]),
            on='SMILES', 
            how='left',
            suffixes=('', '_ext')
        )
        
        # Fill only if main has NA and external has value
        mask = main_df[target_column].isna() & merged[f'{target_column}_ext'].notna()
        filled_count = mask.sum()
        
        if filled_count > 0:
            main_df.loc[mask, target_column] = merged.loc[mask, f'{target_column}_ext']
            lprint(ll.INFO, f"Filled {filled_count} missing values for {target_column}")
        else:
            lprint(ll.INFO, f"No missing values filled for {target_column}")
    
    # === [2] Process new SMILES if requested ===
    if add_new_smiles:
        existing_smiles = set(main_df['SMILES'])
        
        # Get new SMILES with valid target values
        new_smiles = ext_df[
            ~ext_df['SMILES'].isin(existing_smiles) & 
            ext_df[target_column].notna()
        ].copy()
        
        if not new_smiles.empty:
            # Create template with all columns from main_df
            new_rows = pd.DataFrame(columns=main_df.columns)
            
            # Fill SMILES and target column
            new_rows['SMILES'] = new_smiles['SMILES']
            new_rows[target_column] = new_smiles[target_column]
            
            # Set other columns to NA (except ID if exists)
            for col in new_rows.columns:
                if col not in ['SMILES', target_column, 'id']:
                    new_rows[col] = np.nan
            
            # Generate new IDs if column exists
            if 'id' in main_df.columns:
                max_id = main_df['id'].max()
                new_rows['id'] = range(max_id + 1, max_id + 1 + len(new_smiles))
            
            # Safe concatenation by ensuring matching columns
            main_df = pd.concat([
                main_df,
                new_rows[main_df.columns]  # Ensure column order matches
            ], ignore_index=True)
            
            new_smiles_count = len(new_smiles)
            lprint(ll.INFO, f"Added {new_smiles_count} new SMILES entries with valid {target_column} values")
        else:
            lprint(ll.INFO, "No new SMILES with valid target values found to add")
    
    return main_df


def main(df_train, ext_dfs, data_extend_dict) -> pd.DataFrame:
    for target, df_idx in data_extend_dict.items():
        for idx in df_idx:
            ext_df = ext_dfs[idx]

            df_train = update_or_add_smiles_entries(df_train, ext_df, target, add_new_smiles=True, fill_value=True)
    lprint(ll.REPORT, f"Post-merge training samples: {len(df_train)}")
    lprint(ll.REPORT, f"Post-merge missing values:\n{df_train.isna().sum()}")
    lprint(ll.REPORT, f"Post-merge missing values (%):\n{(df_train.isna().mean()*100).round(2)}%")
    lprint(ll.REPORT, f"Post-merge colums:{df_train.columns.to_list()}")
    return df_train

