from Bio.PDB import PDBParser
import numpy as np
import pandas as pd

pdb_path = "2egn.pdb"

# --- Physicochemical maps ---
# Charge proxy at ~physiological pH (simple)
RES_CHARGE = {
    "ASP": -1, "GLU": -1,
    "LYS": +1, "ARG": +1,
    "HIS": +1,  # sometimes 0; simple proxy
}

# Polarity sets (simple)
POLAR_RES = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "HIS"}
CHARGED_RES = {"ASP", "GLU", "LYS", "ARG", "HIS"}

# Aromatic residues
AROMATIC_RES = {"PHE", "TYR", "TRP"}

# Kyte–Doolittle hydrophobicity scale (common, interpretable)
KYTE_DOOLITTLE = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
    "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
    "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

parser = PDBParser(QUIET=True)
structure = parser.get_structure("x", pdb_path)
model = next(structure.get_models())

protein_chain = model["A"]
peptide_chain = model["B"]

# Collect all peptide atom coordinates
peptide_atoms = []
for residue in peptide_chain:
    for atom in residue:
        peptide_atoms.append(atom.get_coord())
peptide_atoms = np.array(peptide_atoms)

# Compute peptide center of mass
peptide_com = peptide_atoms.mean(axis=0)

rows = []

for residue in protein_chain:
    if residue.id[0] != " ":
        continue  # skip hetero residues (waters, ligands, etc.)

    res_id = residue.id[1]
    res_name = residue.resname

    atom_coords = np.array([atom.get_coord() for atom in residue])

    # Min distance to peptide
    dists = np.linalg.norm(
        atom_coords[:, None, :] - peptide_atoms[None, :, :],
        axis=2
    )
    min_dist = dists.min()

    # Distance to peptide center of mass
    res_com = atom_coords.mean(axis=0)
    dist_to_com = np.linalg.norm(res_com - peptide_com)

    # Local density: number of peptide atoms within 8 Å
    local_density = int(np.sum(dists < 8.0))

    # --- Physicochemical features ---
    charge = int(RES_CHARGE.get(res_name, 0))
    is_polar = int(res_name in POLAR_RES or res_name in CHARGED_RES)
    is_charged = int(res_name in CHARGED_RES)

    is_aromatic = int(res_name in AROMATIC_RES)

    # Hydrophobicity score; if unknown residue, fall back to 0.0
    hydrophobicity_kd = float(KYTE_DOOLITTLE.get(res_name, 0.0))
    is_hydrophobic = int(hydrophobicity_kd > 0)

    rows.append({
        "chain": "A",
        "res_id": res_id,
        "res_name": res_name,
        "min_dist_peptide": float(min_dist),
        "dist_to_peptide_com": float(dist_to_com),
        "local_density_8A": local_density,
        "res_charge": charge,
        "is_polar": is_polar,
        "is_charged": is_charged,
        "is_aromatic": is_aromatic,
        "hydrophobicity_kd": hydrophobicity_kd,
        "is_hydrophobic": is_hydrophobic,
    })

df = pd.DataFrame(rows)
df.to_csv("features.csv", index=False)

print(df.head())
print("\nFeature table saved as features_aromaticity_and_hydrophobicity.csv")
