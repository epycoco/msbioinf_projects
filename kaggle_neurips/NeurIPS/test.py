import pandas as pd
from src.config.descriptors import filters

from src.tasks.predict import predict_from_smiles

TARGETS = ['Tg', 'FFV', 'Tc', 'Density', 'Rg']

def main():
    df_test = pd.read_csv("data/neurips-open-polymer-prediction-2025/test.csv", sep=",")
    smiles = df_test.loc[:,"SMILES"].to_list()

    pred_all = {}
    for t in TARGETS:
        predictions, valid_smiles = predict_from_smiles(smiles, t, model_dir=f"./model/{t}")
        print("Valid SMILES:", valid_smiles)
        print("Predictions:", predictions)
        pred_all[t] = predictions

    df_pred = pd.DataFrame(pred_all)
    df_submit = pd.concat([df_test["id"], df_pred], axis=1)

    df_submit.to_csv("submission.csv", index=False)


if __name__ == "__main__":

    main()