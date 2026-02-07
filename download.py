from Bio.PDB import PDBList

pdb_list = PDBList()
pdb_id = "4hhb"  # Example: Hemoglobin
pdb_list.retrieve_pdb_file(pdb_id, pdir="PDB_files", file_format="pdb")
