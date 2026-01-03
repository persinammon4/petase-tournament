import pandas as pd
from seq2ph import net_charge, stability_score, petase_activity

df = pd.read_csv("data/pet-2025-wildtype-cds.csv")
df2 = pd.read_csv("data/solubility_scores.csv")

df["net_charge_pH7"] = df["Wt AA Sequence"].apply(lambda seq: net_charge(seq, 7.0))
df["net_charge_pH5.5"] = df["Wt AA Sequence"].apply(lambda seq: net_charge(seq, 5.5))
df["net_charge_pH9"] = df["Wt AA Sequence"].apply(lambda seq: net_charge(seq, 9.0))
df["stability_score_pH7"] = df["Wt AA Sequence"].apply(lambda seq: stability_score(seq, 7.0))
df["stability_score_pH5.5"] = df["Wt AA Sequence"].apply(lambda seq: stability_score(seq, 5.5))
df["stability_score_pH9"] = df["Wt AA Sequence"].apply(lambda seq: stability_score(seq, 9.0))
df["solubility"] = df2["solubility_scores"].values  
df["predicted_petase_activity_pH7"] = df["Wt AA Sequence"].apply(lambda seq: petase_activity(seq, 7.0))
df["predicted_petase_activity_pH5.5"] = df["Wt AA Sequence"].apply(lambda seq: petase_activity(seq, 5.5))
df["predicted_petase_activity_pH9"] = df["Wt AA Sequence"].apply(lambda seq: petase_activity(seq, 9.0))

df.to_csv("data/wt_attributes.csv", index=False)
df.to_csv("final_data/wt_attributes.csv", index=False)

# with open("wt_sequences.txt", "w") as f:
#     for i, n in enumerate(df["Wt AA Sequence"].values):
#         # Write the FASTA header
#         f.write(f">wt_{i}\n")
#         # Write the sequence, stripping any quotes
#         f.write(n.strip('"') + "\n")
