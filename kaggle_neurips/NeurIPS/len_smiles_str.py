import pandas as pd

df = pd.read_csv("./data/train_all.csv")

smiles = df['SMILES'].to_list()

print(len(smiles))
str_smiles = []

for s in smiles:
    str_smiles.append(len(s))


print (max(str_smiles))

