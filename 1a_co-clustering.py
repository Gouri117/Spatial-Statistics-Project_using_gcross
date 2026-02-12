import pandas as pd
import numpy as np
from sklearn.cluster import SpectralCoclustering
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------- Load Data ----------------
df = pd.read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/cn_sample_matrix_k20_with_outcome.csv")

# Save metadata and isolate numeric matrix
metadata = df[["image", "primary_outcome"]].copy()
metadata["sample_cluster"] = model.row_labels_

cn_matrix = df.drop(columns=["image", "primary_outcome"])
cn_matrix.index = metadata["image"]

# ---------------- Fit Co-Clustering Model ----------------
model = SpectralCoclustering(n_clusters=4, random_state=42)
model.fit(cn_matrix)

# ---------------- Add Cluster Assignments ----------------
# For samples (rows)
metadata["sample_cluster"] = model.row_labels_

# For CNs (columns)
cn_clusters = pd.DataFrame({
    "CN": cn_matrix.columns,
    "cn_cluster": model.column_labels_
})

# ---------------- Save Results ----------------
metadata.to_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/sample_co_clusters.csv", index=False)
cn_clusters.to_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/cn_co_clusters.csv", index=False)

# ---------------- Optional: Visualize Clustered Matrix ----------------
# Rearrange matrix by clustering
fit_data = cn_matrix.iloc[np.argsort(model.row_labels_)]
fit_data = fit_data.iloc[:, np.argsort(model.column_labels_)]

plt.figure(figsize=(12, 8))
sns.heatmap(fit_data, cmap="viridis")
plt.title("Co-clustered CN x Sample matrix")
plt.xlabel("CN")
plt.ylabel("Sample")
plt.tight_layout()
plt.savefig("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/coclustering_heatmap.png")

