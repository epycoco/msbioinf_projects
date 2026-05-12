import os
import json
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from transformers import AutoTokenizer, AutoModel
from sklearn.model_selection import train_test_split
import pandas as pd
from tqdm import tqdm

# ========================================
# Configuration
# ========================================
hf_token = ""
MODEL_NAME = "DeepChem/ChemBERTa-77M-MLM"
tokenizer = AutoTokenizer.from_pretrained("DeepChem/ChemBERTa-77M-MLM", token=hf_token)

EPOCHS = 200
BATCH_SIZE = 16
LR = 1e-4
PATIENCE = 20
MAX_LEN = 128
TARGET_COLUMN = 'Tg'  # singolo target

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
SAVE_DIR = "saved_model"

# ========================================
# Utility function: save everything
# ========================================
def save_model(model, tokenizer, save_dir=SAVE_DIR, config=None):
    os.makedirs(save_dir, exist_ok=True)
    
    # Model
    model_path = os.path.join(save_dir, "chemberta_best.pt")
    torch.save(model.state_dict(), model_path)
    print(f"Model saved at: {model_path}")
    
    # Tokenizer
    tokenizer.save_pretrained(save_dir)
    print(f"Tokenizer saved at: {save_dir}")
    
    # Config
    if config is None:
        config = {
            "max_len": MAX_LEN,
            "batch_size": BATCH_SIZE,
            "target_column": TARGET_COLUMN,
            "model_name": MODEL_NAME
        }
    config_path = os.path.join(save_dir, "config.json")
    with open(config_path, "w") as f:
        json.dump(config, f, indent=4)
    print(f"Configuration saved at: {config_path}")

# ========================================
# Dataset
# ========================================
class ChemDataset(Dataset):
    def __init__(self, smiles, targets, tokenizer, max_len=128):
        self.smiles = smiles.values
        self.targets = targets.fillna(0).values.astype(float)
        self.mask = ~targets.isna().values
        self.tokenizer = tokenizer
        self.max_len = max_len

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, idx):
        sm = str(self.smiles[idx])
        target = torch.tensor(self.targets[idx], dtype=torch.float)
        mask = torch.tensor(self.mask[idx], dtype=torch.bool)
        encoding = self.tokenizer(sm, padding='max_length', truncation=True,
                                  max_length=self.max_len, return_tensors='pt')
        return {
            'input_ids': encoding['input_ids'].squeeze(),
            'attention_mask': encoding['attention_mask'].squeeze(),
            'targets': target,
            'mask': mask
        }

# ========================================
# Single-target Regressor
# ========================================
class ChemBERTaRegressor(nn.Module):
    def __init__(self, base_model):
        super().__init__()
        self.bert = base_model
        self.dropout = nn.Dropout(0.3)
        self.regressor = nn.Linear(self.bert.config.hidden_size, 1)

    def forward(self, input_ids, attention_mask):
        outputs = self.bert(input_ids=input_ids, attention_mask=attention_mask)
        cls_emb = outputs.last_hidden_state[:, 0, :]
        x = self.dropout(cls_emb)
        return self.regressor(x)

# ========================================
# Early Stopping
# ========================================
class EarlyStopping:
    def __init__(self, patience=5, min_delta=0.0):
        self.patience = patience
        self.min_delta = min_delta
        self.best_loss = float('inf')
        self.counter = 0
        self.early_stop = False

    def __call__(self, val_loss, model):
        if val_loss < self.best_loss - self.min_delta:
            self.best_loss = val_loss
            self.counter = 0
            torch.save(model.state_dict(), "best_model.pt")
            print(f"Saved best model with validation loss: {val_loss:.4f}")
        else:
            self.counter += 1
            print(f"No improvement, counter: {self.counter}/{self.patience}")
            if self.counter >= self.patience:
                self.early_stop = True
                print("Early stopping triggered")

# ========================================
# Load Data
# ========================================
smiles_df = pd.read_csv("output/train.ismiles", sep="\t")[[ 'smiles' ]]
targets_df = pd.read_csv("output/train.tt", sep="\t")[[TARGET_COLUMN]]
targets_df = targets_df.apply(pd.to_numeric, errors='coerce').replace(9999, pd.NA)

train_idx, val_idx = train_test_split(range(len(smiles_df)), test_size=0.2, random_state=42)

# ========================================
# Load ChemBERTa
# ========================================
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME, token=hf_token)
base_model = AutoModel.from_pretrained(MODEL_NAME, token=hf_token)

# ========================================
# Initialize Model
# ========================================
model = ChemBERTaRegressor(base_model).to(DEVICE)
optimizer = torch.optim.Adam(model.parameters(), lr=LR)
loss_fn = nn.MSELoss(reduction='none')
early_stopping = EarlyStopping(patience=PATIENCE)

train_dataset = ChemDataset(smiles_df.iloc[train_idx], targets_df.iloc[train_idx], tokenizer)
val_dataset = ChemDataset(smiles_df.iloc[val_idx], targets_df.iloc[val_idx], tokenizer)

train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE)

# ========================================
# Training Loop
# ========================================
for epoch in range(EPOCHS):
    model.train()
    train_loss = 0
    for batch in tqdm(train_loader, desc=f"Epoch {epoch+1}/{EPOCHS}"):
        optimizer.zero_grad()
        input_ids = batch['input_ids'].to(DEVICE)
        attention_mask = batch['attention_mask'].to(DEVICE)
        targets = batch['targets'].to(DEVICE)
        mask = batch['mask'].to(DEVICE)

        outputs = model(input_ids, attention_mask)
        loss = loss_fn(outputs, targets)
        loss = (loss * mask.float()).sum() / mask.float().sum()
        loss.backward()
        optimizer.step()
        train_loss += loss.item()

    avg_train_loss = train_loss / len(train_loader)

    model.eval()
    val_loss_total = 0
    with torch.no_grad():
        for batch in val_loader:
            input_ids = batch['input_ids'].to(DEVICE)
            attention_mask = batch['attention_mask'].to(DEVICE)
            targets = batch['targets'].to(DEVICE)
            mask = batch['mask'].to(DEVICE)

            outputs = model(input_ids, attention_mask)
            loss = loss_fn(outputs, targets)
            loss = (loss * mask.float()).sum() / mask.float().sum()
            val_loss_total += loss.item()

    avg_val_loss = val_loss_total / len(val_loader)
    print(f"Epoch {epoch+1}: Train Loss={avg_train_loss:.4f}, Val Loss={avg_val_loss:.4f}")

    early_stopping(avg_val_loss, model)
    if early_stopping.early_stop:
        print("Stopping training early.")
        break

# ========================================
# Save final model + tokenizer + config
# ========================================
save_model(model, tokenizer)