import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import MiniBatchKMeans

# -------------------- Load Data --------------------
file_path = "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_merged.csv"
df = pd.read_csv(file_path)

# -------------------- Split Metadata and Features --------------------
meta_cols = ["cell_id", "image", "window_center_id", "primary_outcome", "x", "y", "CELL_TYPE"]
feature_cols = [col for col in df.columns if col not in meta_cols]

X = df[feature_cols].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# -------------------- Define K values --------------------
k_values_to_save = [20, 25, 30, 40, 45]

# -------------------- Run MiniBatchKMeans for each K and assign labels --------------------
for k in k_values_to_save:
    print(f"Clustering with K={k}")
    model = MiniBatchKMeans(n_clusters=k, batch_size=100, random_state=42, n_init=20)
    df[f"cluster_k{k}"] = model.fit_predict(X_scaled)

# -------------------- Save the updated DataFrame --------------------
output_path = "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_clustered_multiple_K.csv"
df.to_csv(output_path, index=False)
print(f"\n Saved clustered file: {output_path}")
