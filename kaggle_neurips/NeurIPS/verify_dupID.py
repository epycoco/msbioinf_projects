import pandas as pd


ifeat_df = pd.read_csv('./data/train_fp.csv', header=0)
print (ifeat_df.shape)

ifeat_df = ifeat_df.drop_duplicates(subset='id')
print (ifeat_df.shape)