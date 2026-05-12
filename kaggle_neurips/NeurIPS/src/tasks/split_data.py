from src.logging import lprint, LoggingLevels as ll
from src.config.config import *

def separate_subtables(train_df):
	labels = TARGETS
	subtables = {}
	for label in labels:
		subtables[label] = train_df[['SMILES', label]][train_df[label].notna()]
	return subtables


def main(train_df):
    subtables = separate_subtables(train_df)
    return subtables

