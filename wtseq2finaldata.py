import pandas as pd

df = pd.read_csv("data/predictive-pet-zero-shot-test-2025.csv")
df1 = pd.read_csv("data/petase_families_pca.csv")
df2 = pd.DataFrame()

sequence2cluster = dict()

cluster2rankings = {"1": 1.0, "0":0.9, "2":0.8, "4":0.7, "3":0.6}
max_activity_1 = 0.7
max_activity_2 = 1.4
max_expression = 15

for _, row in df1.iterrows():
    sequence = row["sequence"]
    cluster = str(row["cluster_id"])
    sequence2cluster[sequence] = cluster

df2["sequence"] = df["sequence"]

# Assign scaled attributes
df2["activity_1"] = df["sequence"].apply(
    lambda seq: cluster2rankings[sequence2cluster[seq]] * max_activity_1
)

df2["activity_2"] = df["sequence"].apply(
    lambda seq: cluster2rankings[sequence2cluster[seq]] * max_activity_2
)

df2["expression"] = df["sequence"].apply(
    lambda seq: cluster2rankings[sequence2cluster[seq]] * max_expression
)

df2.to_csv("data/wt_attributes_percentile_scaled.csv", index=False)
df2.to_csv("final_data/wt_attributes_percentile_scaled.csv", index=False)