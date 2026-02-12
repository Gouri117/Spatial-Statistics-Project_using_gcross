
import pandas as pd
from sklearn.cluster import MiniBatchKMeans

# -------------------- Load Data --------------------
file_path = "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_full_matrix_by_window_all.csv"
auc = pd.read_csv(file_path)

# -------------------- Verify cell_id --------------------
# If cell_id is missing, create it
if "cell_id" not in auc.columns:
    auc["cell_id"] = range(len(auc))

# -------------------- Prepare Feature Matrix --------------------
# Exclude metadata columns
feature_cols = [col for col in auc.columns if col not in ["cell_id", "image", "window_center_id"]]
features = auc[feature_cols]

# -------------------- Run MiniBatch KMeans --------------------
best_k = 20
mbkm = MiniBatchKMeans(n_clusters=best_k, batch_size=1000, n_init=5, max_iter=100, random_state=42)
cluster_labels = mbkm.fit_predict(features)

# -------------------- Add Cluster Assignments --------------------
auc["neighborhood_cluster"] = cluster_labels

# -------------------- Save Output --------------------
out_path = "/nfs/turbo/umms-ukarvind/shared_data/mfIHC_collab_data/stanford_crc/results/g_cross_auc_clustered_k20_all.csv"
auc.to_csv(out_path, index=False)

print(f"Clustered AUC matrix saved to: {out_path}")
