import math

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

# PETase activity model (μmol TPA/min·mg enzyme)
def petase_activity(sequence, pH, max_activity=1.0, optimal_pH=8.0, sigma=1.5):
    stability = stability_score(sequence, pH)
    # Gaussian factor for pH effect
    ph_factor = math.exp(-((pH - optimal_pH)**2) / (2 * sigma**2))
    activity = max_activity * stability * ph_factor
    return activity

# Example usage
# sequence = "ACDEK"
# ph = 7.0
# print("Net charge:", net_charge(sequence, ph))
# print("Stability score:", stability_score(sequence, ph))
# print("Predicted PETase activity (μmol TPA/min·mg E):", petase_activity(sequence, ph))
