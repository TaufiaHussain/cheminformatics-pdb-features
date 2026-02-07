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

    # Pairwise distances: (n_res_atoms, n_pep_atoms)
    dists = np.linalg.norm(
        atom_coords[:, None, :] - peptide_atoms[None, :, :],
        axis=2
    )

    min_dist = float(dists.min())

    # Residue "center-of-mass" proxy (mean of atom coords)
    res_com = atom_coords.mean(axis=0)
    dist_to_com = float(np.linalg.norm(res_com - peptide_com))

    # Local density: peptide atoms within 8 Å (as you had)
    local_density_8A = int(np.sum(dists < 8.0))

    # NEW: interface flag
    interface_5A = int(min_dist < 5.0)

    # NEW: per-residue contact counts within 4 Å and 6 Å
    # We count the number of residue atoms that have ANY peptide atom within the cutoff.
    # This is stable and interpretable.
    contact_atoms_4A = int(np.sum(np.any(dists < 4.0, axis=1)))
    contact_atoms_6A = int(np.sum(np.any(dists < 6.0, axis=1)))

    # NEW: mean B-factor (flexibility proxy)
    # Biopython stores B-factor per atom; average over residue atoms
    b_factors = [atom.get_bfactor() for atom in residue]
    mean_bfactor = float(np.mean(b_factors)) if b_factors else float("nan")

    return {
        "min_dist_peptide": min_dist,
        "dist_to_peptide_com": dist_to_com,
        "local_density_8A": local_density_8A,
        "interface_5A": interface_5A,
        "contact_atoms_4A": contact_atoms_4A,
        "contact_atoms_6A": contact_atoms_6A,
        "mean_bfactor": mean_bfactor,
    }
