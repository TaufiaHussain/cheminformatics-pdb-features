import numpy as np


def collect_peptide_atoms(peptide_chain):
    atoms = []
    for residue in peptide_chain:
        for atom in residue:
            atoms.append(atom.get_coord())
    return np.array(atoms)


def compute_peptide_com(peptide_atoms):
    return peptide_atoms.mean(axis=0)


def compute_geometry_features(residue, peptide_atoms, peptide_com):
    atom_coords = np.array([atom.get_coord() for atom in residue])

    # Pairwise distances
    dists = np.linalg.norm(
        atom_coords[:, None, :] - peptide_atoms[None, :, :],
        axis=2
    )

    min_dist = dists.min()

    # Residue center of mass
    res_com = atom_coords.mean(axis=0)
    dist_to_com = np.linalg.norm(res_com - peptide_com)

    # Local density: peptide atoms within 8 Å
    local_density = int(np.sum(dists < 8.0))

    return {
        "min_dist_peptide": float(min_dist),
        "dist_to_peptide_com": float(dist_to_com),
        "local_density_8A": local_density,
    }
