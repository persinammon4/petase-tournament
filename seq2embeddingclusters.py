# THIS IS NOT BEING USED, USING PCA AND K-MEANS IN IPYNB NOTEBOOK

# ==========================================
# PETase Sequence Clustering Pipeline (ProtBERT, unsupervised)
# ==========================================

import pandas as pd
import numpy as np
import torch
from transformers import BertTokenizer, BertModel
import umap
import hdbscan
import matplotlib.pyplot as plt

# -------------------------------
# Step 1: Load sequences
# -------------------------------
# CSV format: columns = ['sequence_id','sequence']
df = pd.read_csv("data/petase-sequences.txt")  # replace with your file
sequences = df['sequence'].tolist()

# -------------------------------
# Step 2: Load ProtBERT
# -------------------------------
tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)
model = BertModel.from_pretrained("Rostlab/prot_bert")
model.eval()

def get_protbert_embedding(seq):
    """Return the sequence-level embedding for a protein."""
    seq = " ".join(seq)  # ProtBERT expects space-separated amino acids
    tokens = tokenizer(seq, return_tensors="pt")
    with torch.no_grad():
        outputs = model(**tokens)
    return outputs.pooler_output.squeeze().numpy()  # 1024-dim vector

# -------------------------------
# Step 3: Generate embeddings
# -------------------------------
print("Generating ProtBERT embeddings for all sequences...")
embeddings = np.array([get_protbert_embedding(seq) for seq in sequences])
print("Embeddings generated:", embeddings.shape)

# -------------------------------
# Step 4: Dimensionality reduction (UMAP)
# -------------------------------
print("Reducing dimensionality with UMAP...")
reducer = umap.UMAP(n_neighbors=5, min_dist=0.1, metric='cosine', random_state=42)
emb_2d = reducer.fit_transform(embeddings)
print("UMAP reduction complete:", emb_2d.shape)

# -------------------------------
# Step 5: Clustering (HDBSCAN)
# -------------------------------
print("Clustering sequences with HDBSCAN...")
clusterer = hdbscan.HDBSCAN(min_cluster_size=2, metric='euclidean')
cluster_labels = clusterer.fit_predict(emb_2d)
df['cluster_id'] = cluster_labels
print("Clusters found:", np.unique(cluster_labels))

# -------------------------------
# Step 6: Save results
# -------------------------------
df.to_csv("data/petase_clusters_unsupervised.csv", index=False)
print("Pipeline complete! Output saved to 'data/petase_clusters_unsupervised.csv'")

# -------------------------------
# Step 7: Optional visualization
# -------------------------------
plt.figure(figsize=(8,6))
for cluster_id in np.unique(cluster_labels):
    idx = cluster_labels == cluster_id
    label = f"Cluster {cluster_id}" if cluster_id != -1 else "Noise"
    plt.scatter(emb_2d[idx,0], emb_2d[idx,1], label=label)
plt.legend()
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("PETase sequences clustered by ProtBERT embeddings (unsupervised)")
plt.show()
