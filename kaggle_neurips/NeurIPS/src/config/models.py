from src.config.config import SEED

model_config = {
    "Tg": {
        "xgb": {
            'n_estimators': 2173,
            'learning_rate': 0.0672418745539774,
            'max_depth': 6,
            "reg_lambda": 5.545520219149715,
            'random_state': SEED
        },
        "lgb": {
            'n_estimators': 46820,
            'max_depth': 8,
            'num_leaves': 146,
            'min_child_samples': 38,
            'min_gain_to_split': 0.19026687141998727,
            'learning_rate': 0.023300052655922467,
            'subsample': 0.917793716129109,
            'colsample_bytree': 0.6754104432113612,
            'reg_lambda': 2.8777986886787335,
            'reg_alpha': 4.066237690162047e-06,
            'verbose': -1
        },
        "cat": {
            'iterations': 2000,
            'learning_rate': 0.05,
            'depth': 6,
            'l2_leaf_reg': 3.0,
            'random_seed': SEED,
            'verbose': 0
        },
        "rff": {
            'n_estimators': 1000,
            'max_depth': 6,
            'min_samples_split': 5,
            'min_samples_leaf': 2,
            'random_state': SEED,
            'n_jobs': -1
        }
    },
    "FFV": {
        "xgb": {
            'n_estimators': 2202,
            'learning_rate': 0.07220580588586338,
            'max_depth': 4,
            "reg_lambda": 2.8872976032666493,
            'random_state': SEED
        },
        "lgb": {
            'n_estimators': 2055,
            'max_depth': 12,
            'num_leaves': 113,
            'min_child_samples': 86,
            'min_gain_to_split': 6.877341773175484e-08,
            'learning_rate': 0.020883673346212925,
            'subsample': 0.662301859231842,
            'colsample_bytree': 0.74050226305929,
            'reg_lambda': 9.214041793767644e-08,
            'reg_alpha': 3.668184812446153e-08,
            'verbose': -1
        },
        "cat": {
            'iterations': 1500,
            'learning_rate': 0.07,
            'depth': 4,
            'l2_leaf_reg': 3.0,
            'random_seed': SEED,
            'verbose': 0
        },
        "rff": {
            'n_estimators': 800,
            'max_depth': 4,
            'min_samples_split': 5,
            'min_samples_leaf': 2,
            'random_state': SEED,
            'n_jobs': -1
        }
    },
    "Tc": {
        "xgb": {
            'n_estimators': 1488,
            'learning_rate': 0.010456188013762864,
            'max_depth': 5,
            "reg_lambda": 9.970345982204618,
            'random_state': SEED
        },
        "lgb": {
            'n_estimators': 2594,
            'max_depth': 6,
            'num_leaves': 91,
            'min_child_samples': 49,
            'min_gain_to_split': 0.0012110933106124974,
            'learning_rate': 0.08638922129126234,
            'subsample': 0.820046754327763,
            'colsample_bytree': 0.8307665860595937,
            'reg_lambda': 6.487495732828185e-08,
            'reg_alpha': 0.05455567427687097,
            'verbose': -1
        },
        "cat": {
            'iterations': 1800,
            'learning_rate': 0.03,
            'depth': 5,
            'l2_leaf_reg': 5.0,
            'random_seed': SEED,
            'verbose': 0
        },
        "rff": {
            'n_estimators': 900,
            'max_depth': 5,
            'min_samples_split': 5,
            'min_samples_leaf': 2,
            'random_state': SEED,
            'n_jobs': -1
        }
    },
    "Density": {
        "xgb": {
            'n_estimators': 1958,
            'learning_rate': 0.10955287548172478,
            'max_depth': 5,
            "reg_lambda": 3.074470087965767,
            'random_state': 4
        },
        "lgb": {
            'n_estimators': 2055,
            'max_depth': 12,
            'num_leaves': 113,
            'min_child_samples': 86,
            'min_gain_to_split': 6.877341773175484e-08,
            'learning_rate': 0.020883673346212925,
            'subsample': 0.662301859231842,
            'colsample_bytree': 0.74050226305929,
            'reg_lambda': 9.214041793767644e-08,
            'reg_alpha': 3.668184812446153e-08,
            'verbose': -1
        },
        "cat": {
            'iterations': 1500,
            'learning_rate': 0.07,
            'depth': 5,
            'l2_leaf_reg': 3.0,
            'random_seed': 4,
            'verbose': 0
        },
        "rff": {
            'n_estimators': 800,
            'max_depth': 5,
            'min_samples_split': 5,
            'min_samples_leaf': 2,
            'random_state': SEED,
            'n_jobs': -1
        }
    },
    "Rg": {
        "xgb": {
            'n_estimators': 1958,
            'learning_rate': 0.10955287548172478,
            'max_depth': 5,
            "reg_lambda": 3.074470087965767,
            'random_state': 4
        },
        "lgb": {
            'n_estimators': 69978,
            'max_depth': 7,
            'num_leaves': 46,
            'min_child_samples': 51,
            'min_gain_to_split': 0.44253308570110284,
            'learning_rate': 0.003043998797424937,
            'subsample': 0.9146909819391108,
            'colsample_bytree': 0.7492210909652985,
            'reg_lambda': 1.5492503391864046e-08,
            'reg_alpha': 0.0003550612286494759,
            'verbose': -1
        },
        "cat": {
            'iterations': 2000,
            'learning_rate': 0.05,
            'depth': 5,
            'l2_leaf_reg': 3.0,
            'random_seed': 4,
            'verbose': 0
        },
        "rff": {
            'n_estimators': 1000,
            'max_depth': 5,
            'min_samples_split': 5,
            'min_samples_leaf': 2,
            'random_state': SEED,
            'n_jobs': -1
        }
    }
}