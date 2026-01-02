import math

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

def net_charge(sequence, pH):
    charge = 0

    # N-terminus
    charge += 1 / (1 + 10**(pH - pKa["N_term"]))

    # C-terminus
    charge -= 1 / (1 + 10**(pKa["C_term"] - pH))

    for aa in sequence:
        if aa in ["D", "E", "C", "Y"]:
            charge -= 1 / (1 + 10**(pKa[aa] - pH))
        elif aa in ["K", "R", "H"]:
            charge += 1 / (1 + 10**(pH - pKa[aa]))

    return charge
