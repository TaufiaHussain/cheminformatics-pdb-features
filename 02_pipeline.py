from Bio.PDB import PDBParser
import pandas as pd

from features.geometry import (
    collect_peptide_atoms,
    compute_peptide_com,
    compute_geometry_features,
)
from features.physicochemical import compute_physicochemical_features


PDB_PATH = "2egn.pdb"

# Heuristic thresholds (can be config later)
MIN_PEPTIDE_LEN = 2
MAX_PEPTIDE_LEN = 30
MIN_PROTEIN_LEN = 50


def is_standard_residue(res):
    return res.id[0] == " "


def count_standard_residues(chain):
    return sum(1 for r in chain if is_standard_residue(r))


def classify_chains(model):
    """
    Returns:
      protein_chains: list of chain objects
      peptide_chains: list of chain objects
    """
    protein_chains = []
    peptide_chains = []

    for chain in model:
        n = count_standard_residues(chain)
        if MIN_PEPTIDE_LEN <= n <= MAX_PEPTIDE_LEN:
            peptide_chains.append(chain)
        elif n >= MIN_PROTEIN_LEN:
            protein_chains.append(chain)

    return protein_chains, peptide_chains


def main():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", PDB_PATH)
    model = next(structure.get_models())

    protein_chains, peptide_chains = classify_chains(model)

    if not peptide_chains:
        raise ValueError("No peptide-like chains found. Choose a protein–peptide complex PDB.")
    if not protein_chains:
        raise ValueError("No protein-like chains found.")

    rows = []

    # For each peptide chain, compute features for each residue in each protein chain
    for pep_chain in peptide_chains:
        peptide_atoms = collect_peptide_atoms(pep_chain)
        peptide_com = compute_peptide_com(peptide_atoms)

        for prot_chain in protein_chains:
            for residue in prot_chain:
                if not is_standard_residue(residue):
                    continue

                res_id = residue.id[1]
                icode = residue.id[2].strip() or ""   # insertion code (optional)
                res_name = residue.resname

                geom = compute_geometry_features(residue, peptide_atoms, peptide_com)
                phys = compute_physicochemical_features(res_name)

                row = {
                    "protein_chain": prot_chain.id,
                    "peptide_chain": pep_chain.id,
                    "res_id": res_id,
                    "ins_code": icode,
                    "res_name": res_name,
                }
                row.update(geom)
                row.update(phys)
                rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv("features_all_chains.csv", index=False)

    print(df.head())
    print(f"\nSaved: features_all_chains.csv  (rows={len(df)})")
    print(f"Protein chains: {[c.id for c in protein_chains]}")
    print(f"Peptide chains: {[c.id for c in peptide_chains]}")


if __name__ == "__main__":
    main()
