from Bio.PDB import PDBParser
import numpy as np
import pandas as pd

pdb_path = "2egn.pdb"

# --- Residue property maps (simple, ML-friendly) ---
# Charge at ~physiological pH (rough but standard for feature engineering)
RES_CHARGE = {
    "ASP": -1, "GLU": -1,
    "LYS": +1, "ARG": +1,
    "HIS": +1,  # sometimes treated as 0; keep +1 for a simple proxy
}

# Polarity (binary)
POLAR_RES = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "HIS"}
CHARGED_RES = {"ASP", "GLU", "LYS", "ARG", "HIS"}
# everything else we'll treat as non-polar by default (ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, GLY)

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
    local_density = np.sum(dists < 8.0)

    # --- NEW: charge + polarity features ---
    charge = RES_CHARGE.get(res_name, 0)

    # polarity as simple binary + optional category
    is_polar = int(res_name in POLAR_RES or res_name in CHARGED_RES)
    is_charged = int(res_name in CHARGED_RES)

    rows.append({
        "chain": "A",
        "res_id": res_id,
        "res_name": res_name,
        "min_dist_peptide": float(min_dist),
        "dist_to_peptide_com": float(dist_to_com),
        "local_density_8A": int(local_density),
        "res_charge": int(charge),
        "is_polar": int(is_polar),
        "is_charged": int(is_charged),
    })

df = pd.DataFrame(rows)
df.to_csv("features.csv", index=False)

print(df.head())
print("\nFeature table saved as features_charge_and_polarity.csv")
