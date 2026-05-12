import pandas as pd
import math
from mordred import Calculator, descriptors

"""idf = pd.read_csv('./data/rdkit2.csv')

idx = idf['id']
df = idf.drop(columns='id')

idf = pd.concat([idx, df], axis=1)
# print (idf.columns.to_list()) # quando crei una serie da un dataset, viene salvata anche il nome della feature
id_to_drop = set()
for d in df.columns.to_list():
    d_feat = df[d]
    itd = set(idf.loc[d_feat.isna(), 'id'].unique().tolist()) # questo comando permette di individuare i NaN ed inserire in una lista i valori di un'altra feature associati al Nan (in questo caso gli id)
    id_to_drop.update(itd)

    print(type(id_to_drop), len(id_to_drop), id_to_drop) 

    """
"""
tg_df = pd.read_csv('./data/target_train/Tg_train.csv', header=0)
variance = tg_df["Tg"].var()
print(f'Variance: {variance}')


from scipy import stats
# supponiamo di avere un DataFrame df con colonna "feature"
data = tg_df["Tg"].dropna()

# 1) Shapiro-Wilk
stat, p = stats.shapiro(data)
print("Shapiro-Wilk: statistic = %.4f, p-value = %.4f" % (stat, p))
if p > 0.05:
    print("Non si rifiuta H0 → i dati possono essere considerati normali")
else:
    print("Rifiuto H0 → i dati non sono normali")

# 2) D’Agostino K²
stat2, p2 = stats.normaltest(data)
print("D’Agostino K²: statistic = %.4f, p-value = %.4f" % (stat2, p2))



import plotly.express as px
def plot_feature_histogram(df: pd.DataFrame, feature_col:str, nbins=None, title=None):
    # Assicuriamoci di eliminare i NaN
    data = df[feature_col].dropna()

    fig = px.histogram(
        data_frame = df,        # o puoi passare data direttamente come serie
        x = feature_col,        
        nbins = nbins,           # numero di bin (opzionale)
        title = title or f"Histogram of {feature_col}"
    )
    fig.update_layout(
        xaxis_title = feature_col,
        yaxis_title = "Count",
        bargap = 0.1
    )
    fig.show()

# Esempio d’uso:
plot_feature_histogram(tg_df, "Tg", nbins=30, title="Distribuzione di Tg")"""



