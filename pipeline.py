from Bio.PDB import PDBParser
import pandas as pd

from features.geometry import (
    collect_peptide_atoms,
    compute_peptide_com,
    compute_geometry_features,
)
from features.physicochemical import compute_physicochemical_features


pdb_path = "2egn.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("x", pdb_path)
model = next(structure.get_models())

protein_chain = model["A"]
peptide_chain = model["B"]

# Prepare peptide data
peptide_atoms = collect_peptide_atoms(peptide_chain)
peptide_com = compute_peptide_com(peptide_atoms)

rows = []

for residue in protein_chain:
    if residue.id[0] != " ":
        continue

    res_id = residue.id[1]
    res_name = residue.resname

    geom = compute_geometry_features(residue, peptide_atoms, peptide_com)
    phys = compute_physicochemical_features(res_name)

    row = {
        "chain": "A",
        "res_id": res_id,
        "res_name": res_name,
    }
    row.update(geom)
    row.update(phys)

    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv("features_.csv", index=False)

print(df.head())
print("\nFeature table saved as features.csv")
