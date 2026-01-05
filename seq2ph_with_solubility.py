import math
import pandas as pd

# pKa values for ionizable groups
pKa = {
    "N_term": 9.6,
    "C_term": 2.4,
    "D": 3.9,
    "E": 4.1,
    "H": 6.0,
    "C": 8.3,
    "Y": 10.1,
    "K": 10.5,
    "R": 12.5
}

# Calculate net charge at a given pH
def net_charge(sequence, pH):
    charge = 0
    charge += 1 / (1 + 10**(pH - pKa["N_term"]))  # N-term
    charge -= 1 / (1 + 10**(pKa["C_term"] - pH))  # C-term

    for aa in sequence:
        if aa in ["D", "E", "C", "Y"]:
            charge -= 1 / (1 + 10**(pKa[aa] - pH))
        elif aa in ["K", "R", "H"]:
            charge += 1 / (1 + 10**(pH - pKa[aa]))
    return charge

# Stability score based on net charge
def stability_score(sequence, pH, k=0.02):
    Z = net_charge(sequence, pH)
    score = max(0, 1 - k * Z**2)
    return score

def solubility_score(sequence, pH, a=0.8):
    Z = abs(net_charge(sequence, pH))
    return 1 - math.exp(-a * Z)

def solubility_score_saved(sequence):
    score = 1.0
    df = pd.read_csv("data/solubility_scores.csv")

    fastid2seq = {}

    with open("wt_sequences.txt", "r") as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            fasta_id = lines[i].strip().lstrip(">")
            seq = lines[i+1].strip()
            fastid2seq[fasta_id] = seq
    
    seq_list = df["fasta_id"].values
    score_list = df["solubility_scores"].values

    for i, fasta_id in enumerate(seq_list):
        if fastid2seq[fasta_id] == sequence:
            score = score_list[i]
            break
    return score

def catalytic_ph_factor(pH, optimal_pH=8.0, sigma=1.2):
    return math.exp(-((pH - optimal_pH)**2) / (2 * sigma**2))

# PETase activity model (μmol TPA/min·mg enzyme)
def petase_activity(sequence,
                     pH,
                     max_activity=1.0,   # μmol TPA / min·mg <= should this be looked up from another dataset, right now totals never exceed 1 [unit]
                     k_stability=0.02,
                     a_solubility=0.8,
                     optimal_pH=8.0,
                     sigma=1.2):

    stability = stability_score(sequence, pH, k_stability) # refine based on literature, depends on net charge
    solubility = solubility_score_saved(sequence)
    catalysis = catalytic_ph_factor(pH, optimal_pH, sigma) # <= is this necessary?

    print(f"solubility score calculated: {solubility_score(sequence, pH)}, solubility score saved {solubility_score_saved(sequence)}")
    print(f"stability: {stability}, catalytic ph factor: {catalysis}")

    activity = max_activity * stability * solubility * catalysis
    return activity

# Example usage
sequence = "AADNPYQRGPDPTNASIEAATGPFAVGTQPIVGASGFGGGQIYYPTDTSQTYGAVVIVPGFISVWAQLNWLGPRLASQGFVVIGIETSVITDLPDPRGDQALAALDWATTRSPVASRIDRTRLAAAGWSMGGGGLRRAALQRPSLKAIVGMAPWNGERNWSAVTVPTLFFGGSSDAVASPNDHAKPFYNSITRAEKDYIELRNADHFFPTSANTTMAKYFISWLKRWVDNDTRYTQFLCPGPSTGLFAPVSASMNTCPF"
ph = 7.0
print("Net charge:", net_charge(sequence, ph))
print("Stability score:", stability_score(sequence, ph))
print("Predicted PETase activity (μmol TPA/min·mg E):", petase_activity(sequence, ph))
