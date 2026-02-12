import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

# -------------------- Load Data --------------------
df = pd.read_csv("/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_merged.csv")

# -------------------- Split Metadata and Features --------------------
meta_cols = ["cell_id", "image", "window_center_id", "primary_outcome", "x", "y", "CELL_TYPE"]
feature_cols = [col for col in df.columns if col not in meta_cols]

X = df[feature_cols].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# -------------------- Subsample for Silhouette --------------------
np.random.seed(42)
subset_indices = np.random.choice(X_scaled.shape[0], size=5000, replace=False)
X_subset = X_scaled[subset_indices]

# -------------------- Sweep Over K --------------------
k_values = range(10, 45, 5)
silhouette_scores = []

for k in k_values:

    print(f"Running MiniBatchKMeans with K={k}...")
    model = MiniBatchKMeans(n_clusters=k, batch_size=100, random_state=42, n_init = 10)
    model.fit(X_scaled)
    labels_subset = model.predict(X_subset)
    score = silhouette_score(X_subset, labels_subset)
    silhouette_scores.append(score)
    print(f"K={k} → Silhouette Score = {score:.4f}")

# -------------------- Best K Based on Silhouette --------------------
best_k = k_values[np.argmax(silhouette_scores)]
print(f"\nBest K based on silhouette score: {best_k}")

# -------------------- Cluster Using Best K --------------------
final_model_best = MiniBatchKMeans(n_clusters=best_k, batch_size=100, random_state=42, n_init = 20)
df["cluster_k{}".format(best_k)] = final_model_best.fit_predict(X_scaled)

# Save clustered file using best K
output_best = f"/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_k{best_k}.csv"
df.to_csv(output_best, index=False)
print(f"\nClustered data saved to: {output_best}")

