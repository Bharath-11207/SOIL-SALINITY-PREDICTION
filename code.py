
import pandas as pd
import numpy as np
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from torchdiffeq import odeint_adjoint as odeint
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix, ConfusionMatrixDisplay
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_class_weight
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from imblearn.over_sampling import SMOTE

csv_path = '/content/drive/MyDrive/Soil_Salinity/DL/Salinity_Samples_Merged.csv'
df = pd.read_csv(csv_path)
df['salinity'] = df['salinity'] - 1
df['salinity'] = df['salinity'].apply(lambda x: 0 if x <= 1 else 1)

features = ['SI1', 'SI2', 'SI3', 'SI4', 'SI5', 'SI6', 'SI7', 'SI8', 'SI9', 'NDSI', 'CRSI', 'NDVI']
X = df[features].values.astype(np.float32)
y = df['salinity'].values.astype(int)

# Feature Scaling
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

rf = RandomForestClassifier(random_state=42)
rf.fit(X_scaled, y)
importances = rf.feature_importances_
important_indices = np.argsort(importances)[-8:]
X_scaled = X_scaled[:, important_indices]
selected_features = [features[i] for i in important_indices]

smote = SMOTE(random_state=42)
X_resampled, y_resampled = smote.fit_resample(X_scaled, y)

train_X, test_X, train_y, test_y = train_test_split(
    X_resampled, y_resampled, test_size=0.3, random_state=42, stratify=y_resampled)

class_weights = compute_class_weight('balanced', classes=np.unique(y_resampled), y=y_resampled)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
weights_tensor = torch.tensor(class_weights, dtype=torch.float32).to(device)

class SalinityDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.long)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

train_ds = SalinityDataset(train_X, train_y)
test_ds = SalinityDataset(test_X, test_y)
train_loader = DataLoader(train_ds, batch_size=64, shuffle=True)
test_loader = DataLoader(test_ds, batch_size=64)

class ODEF(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.linear1 = nn.Linear(dim, dim * 2)
        self.act = nn.ReLU()
        self.linear2 = nn.Linear(dim * 2, dim)

    def forward(self, t, x):
        return self.linear2(self.act(self.linear1(x)))

class LiquidNN(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super().__init__()
        self.fc_in = nn.Linear(input_dim, hidden_dim)
        self.bn = nn.BatchNorm1d(hidden_dim)
        self.dropout = nn.Dropout(0.4)
        self.odefunc = ODEF(hidden_dim)
        self.integration_time = torch.tensor([0, 1]).float()
        self.fc_out = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        x = torch.relu(self.bn(self.fc_in(x)))
        x = self.dropout(x)
        x = odeint(self.odefunc, x, self.integration_time.to(x.device))[1]
        return self.fc_out(x)

model = LiquidNN(input_dim=len(selected_features), hidden_dim=256, output_dim=2).to(device)
loss_fn = nn.CrossEntropyLoss(weight=weights_tensor)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

for i, w in enumerate(class_weights):
    print(f"Class {i} weight: {w:.4f}")

train_loss_hist, val_loss_hist = [], []
epochs = 100

patience = 10
min_delta = 0.001
best_val_loss = float('inf')
counter = 0

for epoch in range(epochs):
    model.train()
    total_loss = 0
    for xb, yb in train_loader:
        xb, yb = xb.to(device), yb.to(device)
        pred = model(xb)
        loss = loss_fn(pred, yb)
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        total_loss += loss.item()
    train_loss_hist.append(total_loss / len(train_loader))

    model.eval()
    with torch.no_grad():
        val_loss = sum(loss_fn(model(xb.to(device)), yb.to(device)).item() for xb, yb in test_loader) / len(test_loader)
        val_loss_hist.append(val_loss)

    print(f"Epoch {epoch+1}: Train Loss={train_loss_hist[-1]:.4f} | Val Loss={val_loss_hist[-1]:.4f}")
    if val_loss < best_val_loss - min_delta:
        best_val_loss = val_loss
        counter = 0
    else:
        counter += 1
        if counter >= patience:
            print(f"Early stopping at epoch {epoch+1}")
            break

plt.plot(train_loss_hist, label='Train Loss')
plt.plot(val_loss_hist, label='Validation Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.show()

model.eval()
y_true, y_pred = [], []

with torch.no_grad():
    for xb, yb in test_loader:
        preds = model(xb.to(device)).argmax(1).cpu().numpy()
        y_true.extend(yb.numpy())
        y_pred.extend(preds)

print(classification_report(y_true, y_pred, target_names=['Class 0', 'Class 1']))
cm = confusion_matrix(y_true, y_pred)
print(cm)
ConfusionMatrixDisplay(cm, display_labels=['0', '1']).plot()
plt.show()
