# Charge proxy
RES_CHARGE = {
    "ASP": -1, "GLU": -1,
    "LYS": +1, "ARG": +1,
    "HIS": +1,
}

POLAR_RES = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "HIS"}
CHARGED_RES = {"ASP", "GLU", "LYS", "ARG", "HIS"}
AROMATIC_RES = {"PHE", "TYR", "TRP"}

# Kyte–Doolittle hydrophobicity
KYTE_DOOLITTLE = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
    "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
    "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}


def compute_physicochemical_features(res_name):
    charge = int(RES_CHARGE.get(res_name, 0))
    is_polar = int(res_name in POLAR_RES or res_name in CHARGED_RES)
    is_charged = int(res_name in CHARGED_RES)
    is_aromatic = int(res_name in AROMATIC_RES)

    hydrophobicity = float(KYTE_DOOLITTLE.get(res_name, 0.0))
    is_hydrophobic = int(hydrophobicity > 0)

    return {
        "res_charge": charge,
        "is_polar": is_polar,
        "is_charged": is_charged,
        "is_aromatic": is_aromatic,
        "hydrophobicity_kd": hydrophobicity,
        "is_hydrophobic": is_hydrophobic,
    }
