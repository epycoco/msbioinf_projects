import pandas as pd

test_mordred_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_mordred.csv", header=0)
test_rdk2_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_rdkit2.csv", header=0)
test_rdk3_df = pd.read_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_rdkit3.csv")

test_desc_df = pd.merge(test_mordred_df, test_rdk2_df, on='id')
test_desc_df = pd.merge(test_desc_df, test_rdk3_df, on='id')

test_desc_df.to_csv("C:/Users/andre/Desktop/machine_learning_prjct/FINAL/NeurIPS/data/test_desc.csv")