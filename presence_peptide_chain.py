from Bio.PDB import PDBParser

pdb_path = "2egn.pdb"  # <-- set this
parser = PDBParser(QUIET=True)
structure = parser.get_structure("x", pdb_path)
model = next(structure.get_models())

print("Chain summary:")
for chain in model:
    standard = 0
    hetero = 0
    for res in chain:
        if res.id[0] == " ":
            standard += 1
        else:
            hetero += 1
    print(f"Chain {chain.id}: standard={standard}, hetero={hetero}")
