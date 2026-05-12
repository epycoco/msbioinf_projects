# FILE PATH
TRAIN_PATH = './data/train_all.csv'
TEST_PATH = './data/neurips-open-polymer-prediction-2025/test.csv'
TRAIN_CLEAN_PATH = './data/train_clean.csv'
TEST_CLEAN_PATH = './data/test_clean.csv'

# DATASET
FP_RADIUS = 2
FP_BITS = 256
BERTA_BITS = 384
VAR_THRESH = 0.01
NAN_THRESH = 0.02
REP_THRESH = 0.9

TARGETS = ['Tg', 'FFV', 'Tc', 'Density', 'Rg']

# TRAIN
ALFA_LOW = [0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10]
ALFA_UP  = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]
#  "Cat", 'LGB', "MLP", "KNN",
MODELS = ['XGB', "RFF", "ET", "HGB", "GBR"]
CUSTOM_MODELS = ['XGB', 'LGB']
SEED = 42
OPTIM_TRIALS = 20
K_FOLD = 5
TOP_K = 2
EARLY_STOPPING = 30
MLP_ITER = 500

DEVICE = "cuda"