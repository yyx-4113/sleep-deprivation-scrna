"""
PepDesign-Active ConvVAE Training v2 — ModelScope Workspace
============================================================
Fixed architecture: stride-2 convs + ConvTranspose upsampling (proper mirror).
Lower lr (1e-4), gradient clipping, warmup beta schedule.
"""
import torch, torch.nn as nn, torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import numpy as np, pandas as pd, time, json, random, os

print(f"Device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}")
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# -------------------------------------------------------
# 1. Load
# -------------------------------------------------------
csv_path = '/mnt/workspace/merged_peptide_library.csv'
if not os.path.exists(csv_path):
    print(f"ERROR: {csv_path} not found. Upload merged_peptide_library.csv first!")
    exit(1)

df = pd.read_csv(csv_path)
unique = df.drop_duplicates(subset=['peptide'])
unique = unique[unique['length'] >= 2]
peptides_all = unique['peptide'].tolist()
random.seed(42)
peptides = random.sample(peptides_all, min(52517, len(peptides_all)))
print(f"{len(peptides)} peptides, mean len={np.mean([len(p) for p in peptides]):.1f}")

# -------------------------------------------------------
# 2. One-hot
# -------------------------------------------------------
AA = 'ACDEFGHIKLMNPQRSTVWY'
AA2I = {aa: i for i, aa in enumerate(AA)}
MAX_L, N_AA = 20, 20

def to_onehot(seq):
    m = np.zeros((MAX_L, N_AA), dtype=np.float32)
    for i, aa in enumerate(seq[:MAX_L]):
        if aa in AA2I: m[i, AA2I[aa]] = 1.0
    return m

print("One-hot encoding...")
oh = np.stack([to_onehot(p) for p in peptides])
oh = np.transpose(oh, (0, 2, 1))  # (N, 20, 20)
print(f"Shape: {oh.shape}")

# -------------------------------------------------------
# 3. ConvVAE (proper mirror architecture)
# -------------------------------------------------------
class ConvVAE(nn.Module):
    """Encoder: Conv1d(s=2) -> Conv1d(s=2) -> Pool(1) -> FC -> latent.
       Decoder: FC -> reshape -> ConvTranspose1d(s=2) -> ConvTranspose1d(s=2).
       Mirrors exactly."""
    def __init__(self, latent_dim=64):
        super().__init__()
        # Encoder: 20->10->5->1
        self.enc_conv1 = nn.Conv1d(N_AA, 32, 3, stride=2, padding=1)
        self.enc_conv2 = nn.Conv1d(32, 64, 3, stride=2, padding=1)
        self.enc_bn1 = nn.BatchNorm1d(32)
        self.enc_bn2 = nn.BatchNorm1d(64)
        self.enc_pool = nn.AdaptiveAvgPool1d(1)
        self.enc_fc1 = nn.Linear(64, 128)
        self.enc_fc_mu = nn.Linear(128, latent_dim)
        self.enc_fc_logvar = nn.Linear(128, latent_dim)
        # Decoder: 1->5->10->20
        self.dec_fc1 = nn.Linear(latent_dim, 128)
        self.dec_fc2 = nn.Linear(128, 64 * 5)
        self.dec_dc1 = nn.ConvTranspose1d(64, 32, 3, stride=2, padding=1, output_padding=1)
        self.dec_bn1 = nn.BatchNorm1d(32)
        self.dec_dc2 = nn.ConvTranspose1d(32, N_AA, 3, stride=2, padding=1, output_padding=1)

    def encode(self, x):
        h = F.relu(self.enc_bn1(self.enc_conv1(x)))   # (B,32,10)
        h = F.relu(self.enc_bn2(self.enc_conv2(h)))    # (B,64,5)
        h = self.enc_pool(h).view(h.size(0), -1)       # (B,64)
        h = F.relu(self.enc_fc1(h))                     # (B,128)
        return self.enc_fc_mu(h), self.enc_fc_logvar(h)

    def reparameterize(self, mu, lv):
        return mu + torch.exp(0.5 * lv) * torch.randn_like(mu)

    def decode(self, z):
        h = F.relu(self.dec_fc1(z))                     # (B,128)
        h = F.relu(self.dec_fc2(h)).view(h.size(0), 64, 5)  # (B,64,5)
        h = F.relu(self.dec_bn1(self.dec_dc1(h)))       # (B,32,11) w/ output_padding=1
        h = self.dec_dc2(h)                              # (B,20,21) w/ output_padding=1
        return h[:, :, :MAX_L]                           # trim to (B,20,20)

    def forward(self, x):
        mu, lv = self.encode(x)
        z = self.reparameterize(mu, lv)
        return self.decode(z), mu, lv

model = ConvVAE().to(device)
print(f"Params: {sum(p.numel() for p in model.parameters()):,}")

# -------------------------------------------------------
# 4. Train (lr=1e-4, grad clip, KL warmup)
# -------------------------------------------------------
def vae_loss(recon, x, mu, lv, beta=0.1):
    """Masked CE + beta * KL."""
    B, C, L = x.shape
    tgt = x.argmax(dim=1)
    mask = (x.sum(dim=1) > 0).float()
    r = recon.permute(0, 2, 1).reshape(-1, C)
    ce = F.cross_entropy(r, tgt.reshape(-1), reduction='none')
    bce = (ce * mask.reshape(-1)).sum() / (mask.sum() + 1e-8)
    kl = -0.5 * torch.sum(1 + lv - mu.pow(2) - lv.exp()) / B
    return bce + beta * kl

X = torch.tensor(oh, dtype=torch.float32)
loader = DataLoader(TensorDataset(X), batch_size=128, shuffle=True, num_workers=2, pin_memory=True)
opt = torch.optim.Adam(model.parameters(), lr=1e-4)

print(f"Training: {len(peptides)} peptides, lr=1e-4, grad_clip=1.0")
losses, t0 = [], time.time()

for ep in range(80):
    model.train()
    # KL warmup: beta rises slowly from 0.0 to 0.3 over 40 epochs, stays at 0.3
    beta = min(0.3, 0.3 * ep / 40)
    ep_loss = 0.0
    for (batch,) in loader:
        batch = batch.to(device)
        opt.zero_grad()
        recon, mu, lv = model(batch)
        loss = vae_loss(recon, batch, mu, lv, beta)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)  # gradient clipping
        opt.step()
        ep_loss += loss.item()
    losses.append(ep_loss / len(loader))
    if (ep + 1) % 10 == 0:
        print(f"  Epoch {ep+1:3d}/80  loss={losses[-1]:.4f}  beta={beta:.3f}  {time.time()-t0:.0f}s")

elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s. Loss: {losses[0]:.2f} -> {losses[-1]:.2f}")

# -------------------------------------------------------
# 5. Reconstruction accuracy
# -------------------------------------------------------
model.eval()
test_n = min(5000, len(peptides))
test_t = torch.tensor(oh[:test_n], dtype=torch.float32).to(device)
with torch.no_grad():
    recon_idx = model(test_t)[0].argmax(dim=1).cpu().numpy()

ids = []
for i in range(test_n):
    orig = peptides[i]
    rseq = ''.join(AA[j] for j in recon_idx[i] if j < len(AA))
    if len(orig) > 0:
        ids.append(sum(1 for k in range(min(len(orig), len(rseq))) if orig[k] == rseq[k]) / len(orig))
acc = np.mean(ids)
print(f"Reconstruction accuracy: {acc:.3f} ({acc*100:.1f}%)")

# -------------------------------------------------------
# 6. Latent vectors
# -------------------------------------------------------
print("Extracting latent vectors...")
model.eval()
lats = []
with torch.no_grad():
    for i in range(0, len(oh), 1024):
        batch = torch.tensor(oh[i:i+1024], dtype=torch.float32).to(device)
        mu, _ = model.encode(batch)
        lats.append(mu.cpu().numpy())
lats = np.vstack(lats).astype(np.float32)
print(f"Latents: {lats.shape}")

from sklearn.decomposition import PCA
pca_l = PCA(n_components=64).fit(lats)
var = pca_l.explained_variance_ratio_.sum()
print(f"PCA variance: {var:.3f}")

# -------------------------------------------------------
# 7. Save
# -------------------------------------------------------
out = '/mnt/workspace'
torch.save(model.state_dict(), f'{out}/conv_vae_weights.pt')
np.savez_compressed(f'{out}/vae_latents.npz', latents=lats, peptides=np.array(peptides))
pd.DataFrame({'epoch': range(1, len(losses)+1), 'loss': losses}).to_csv(f'{out}/vae_loss.csv', index=False)

summary = {
    'model': 'ConvVAE_v2', 'latent_dim': 64, 'n_peptides': len(peptides),
    'epochs': 80, 'initial_loss': float(losses[0]), 'final_loss': float(losses[-1]),
    'reconstruction_accuracy': float(acc), 'latent_pca_variance': float(var),
    'training_time_s': float(elapsed), 'device': str(device), 'lr': 1e-4,
}
with open(f'{out}/training_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(range(1, len(losses)+1), losses, 'b-', lw=1.5)
ax.axvline(40, color='gray', ls='--', alpha=0.5, label='KL warmup end')
ax.set_xlabel('Epoch'); ax.set_ylabel('Loss')
ax.set_title(f'ConvVAE v2 Loss (acc={acc*100:.1f}%)'); ax.legend()
plt.tight_layout(); plt.savefig(f'{out}/fig_s1_vae_loss.png', dpi=300)

n_v = min(3000, len(lats))
iv = np.random.choice(len(lats), n_v, replace=False)
l2d = PCA(n_components=2).fit_transform(lats[iv])
KD = {'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3}
colors = [np.mean([KD.get(aa,0) for aa in peptides[i]]) for i in iv]
fig, ax = plt.subplots(figsize=(8, 7))
sc = ax.scatter(l2d[:,0], l2d[:,1], c=colors, cmap='RdYlBu', s=3, alpha=0.5)
plt.colorbar(sc, label='Hydrophobicity')
ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
ax.set_title(f'VAE Latent PCA (v2, var={var:.3f})')
plt.tight_layout(); plt.savefig(f'{out}/fig_s2_latent_pca.png', dpi=300)

print("\n" + "=" * 50)
print("TRAINING COMPLETE")
print("=" * 50)
for k, v in summary.items():
    print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")
print(f"\nFiles in {out}/: conv_vae_weights.pt, vae_latents.npz, vae_loss.csv, training_summary.json, fig_s1/s2")
