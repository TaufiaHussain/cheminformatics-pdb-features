# Protein--Peptide Feature Extraction Pipeline

## Overview

This project implements a **Python-only, modular feature extraction
pipeline** for protein--peptide complexes in PDB format.\
It computes residue-level **geometric** and **physicochemical** features
and outputs a structured table suitable for machine learning workflows.

The pipeline is designed to be: - Modular - Reproducible - Easy to
extend with new features - Compatible with downstream ML pipelines


## Project Structure

    project/
    │
    ├── features/
    │   ├── __init__.py
    │   ├── geometry.py
    │   └── physicochemical.py
    │
    ├── pipeline.py
    ├── validate.py
    └── example.pdb


## Inputs

### Required

-   **PDB file** containing a protein--peptide complex.
-   The pipeline automatically detects:
    -   Protein chains (long chains)
    -   Peptide chains (short chains)

Example input:

    2egn.pdb


## Outputs

### Main output

**features_all_chains.csv**

A residue-level feature table containing:

| Column \| Description \|

\|--------\|-------------\| protein_chain \| Chain ID of the protein \|
\| peptide_chain \| Chain ID of the peptide \| \| res_id \| Residue
number \| \| ins_code \| Insertion code \| \| res_name \| Residue name
\| \| min_dist_peptide \| Minimum atom--atom distance to peptide \| \|
dist_to_peptide_com \| Distance to peptide center of mass \| \|
local_density_8A \| Number of peptide atoms within 8 Å \| \|
interface_5A \| 1 if residue is within 5 Å of peptide \| \|
contact_atoms_4A \| Number of residue atoms within 4 Å \| \|
contact_atoms_6A \| Number of residue atoms within 6 Å \| \|
mean_bfactor \| Mean B-factor of residue atoms \| \| res_charge \|
Residue charge proxy (-1, 0, +1) \| \| is_polar \| Polar residue flag \|
\| is_charged \| Charged residue flag \| \| is_aromatic \| Aromatic
residue flag \| \| hydrophobicity_kd \| Kyte--Doolittle hydrophobicity
score \| \| is_hydrophobic \| Hydrophobic residue flag \|


## How to Run

### Step 1: Install dependencies

    pip install biopython numpy pandas

### Step 2: Place your PDB file in the project folder

Example:

    2egn.pdb

### Step 3: Run the pipeline

    python pipeline.py

Output:

    features_all_chains.csv


## Validation

To run basic sanity checks on the output:

    python validate.py

This checks: - No negative distances - Valid binary flags - Correct
charge values - Non-negative contact counts - Reasonable B-factors


## Feature Categories

### Geometric features

-   Distance to peptide
-   Interface flag
-   Contact counts
-   Local density
-   B-factor (flexibility proxy)

### Physicochemical features

-   Charge
-   Polarity
-   Aromaticity
-   Hydrophobicity


## Typical Workflow

1.  Load PDB structure
2.  Identify protein and peptide chains
3.  Compute geometry features per residue
4.  Compute physicochemical features
5.  Merge into a single feature table
6.  Validate outputs


## Intended Use

This pipeline is suitable for:

-   Residue-level ML prediction tasks
-   Protein--peptide interface prediction
-   Structural bioinformatics workflows
-   Feature engineering for downstream models


## Notes

-   The current implementation uses simple, interpretable
    physicochemical encodings.
-   Feature definitions can be adjusted to match reference pipelines.
-   Designed to be extended with additional structural or sequence-based
    features.


## Author

Taufia Hussain\
Computational Biologist & Python Developer
