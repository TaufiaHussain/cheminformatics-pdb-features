import pandas as pd

FEATURES_PATH = "features_all_chains.csv"
df = pd.read_csv(FEATURES_PATH)

# Required columns check (helps avoid KeyError confusion)
required_cols = [
    "min_dist_peptide", "dist_to_peptide_com", "local_density_8A",
    "res_charge", "is_polar", "is_charged", "is_aromatic", "is_hydrophobic",
    "interface_5A", "contact_atoms_4A", "contact_atoms_6A", "mean_bfactor"
]
missing = [c for c in required_cols if c not in df.columns]
assert not missing, f"Missing expected columns: {missing}"

# Basic sanity checks
assert len(df) > 0, "Feature table is empty."

# Distances must be >= 0
for col in ["min_dist_peptide", "dist_to_peptide_com"]:
    assert (df[col] >= 0).all(), f"Found negative values in {col}"

# Density must be >= 0
assert (df["local_density_8A"] >= 0).all(), "local_density_8A has negative values"

# Charge should be in {-1,0,1} (with our current mapping)
assert df["res_charge"].isin([-1, 0, 1]).all(), "Unexpected charge values"

# Flags should be 0/1
flag_cols = ["is_polar", "is_charged", "is_aromatic", "is_hydrophobic", "interface_5A"]
for c in flag_cols:
    assert df[c].isin([0, 1]).all(), f"{c} contains values other than 0/1"

# Contact counts must be >= 0
for col in ["contact_atoms_4A", "contact_atoms_6A"]:
    assert (df[col] >= 0).all(), f"{col} has negative values"

# B-factor check (allow NaN, otherwise must be >= 0)
bf = df["mean_bfactor"]
assert ((bf.isna()) | (bf >= 0)).all(), "mean_bfactor has negative values (unexpected)"

print("Validation passed: basic sanity checks look good.")
