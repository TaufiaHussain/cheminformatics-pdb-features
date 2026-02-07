from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure("x", "2egn.pdb")

model = next(structure.get_models())
chains = list(model.get_chains())

print("Chains found:", [c.id for c in chains])
for c in chains:
    residues = [r for r in c.get_residues() if r.id[0] == " "]
    print(f"Chain {c.id}: {len(residues)} residues")
