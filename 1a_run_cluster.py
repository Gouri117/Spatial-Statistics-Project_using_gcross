import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

# -------------------- Load Data --------------------
df = pd.read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_merged_with_xy_and_celltype.csv")

# -------------------- Extract AUC Features --------------------
auc_data = df.drop(columns=["cell_id", "image", "window_center_id", "X", "Y", "CELL_TYPE"])
scaled_auc = auc_data.apply(pd.to_numeric, errors='coerce').fillna(0)

# -------------------- Hyperparameter Sweep for K --------------------
k_values = list(range(10, 45, 5))
silhouette_scores = []

for k in k_values:
    print(f"Running MiniBatchKMeans for K = {k}")
    mbk = MiniBatchKMeans(n_clusters=k, batch_size=1000, n_init=10, max_iter=300, random_state=42)
    cluster_labels = mbk.fit_predict(scaled_auc)
    
    score = silhouette_score(scaled_auc, cluster_labels)
    silhouette_scores.append(score)

# -------------------- Plot Silhouette Scores --------------------
plt.figure(figsize=(8, 6))
plt.plot(k_values, silhouette_scores, marker='o', linewidth=2)
plt.title("Silhouette Score vs Number of Clusters (K)")
plt.xlabel("Number of Clusters (K)")
plt.ylabel("Average Silhouette Score")
plt.xticks(k_values)
plt.grid(True)
plt.tight_layout()
plt.savefig("/nfs/turbo/umms-ukarvind/shared_data/gcross_data/example/results/silhouette_k_sweep.png", dpi=300)
plt.close()

# -------------------- Print Best K --------------------
best_k = k_values[np.argmax(silhouette_scores)]
print(f"Best K based on silhouette score: {best_k}")

