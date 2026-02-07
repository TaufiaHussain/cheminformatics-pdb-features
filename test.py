from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure("2egn", "2egn.pdb")

for model in structure:
    for chain in model:
        for residue in chain:
            print(chain.id, residue.resname, residue.id)
