"""
compute_descriptors.py
----------------------
Computes physicochemical descriptors for a panel of PFAS compounds using RDKit,
with a focus on properties relevant to passive membrane permeation.

Scientific rationale
--------------------
Passive permeation of a small molecule across a lipid bilayer is governed largely
by its lipophilicity (partitioning into the hydrophobic core) and its polar surface
area / hydrogen-bonding capacity (desolvation penalty at the headgroup). For PFAS,
the perfluorinated tail drives lipophilicity while the ionisable headgroup
(carboxylate / sulfonate) creates a strong desolvation penalty. We therefore expect
permeation propensity to scale with fluorinated chain length, modulated by
headgroup chemistry.

This script does NOT compute a true free-energy profile (PMF) — that requires
molecular dynamics with enhanced sampling (e.g. the AWH protocol). It produces a
fast, interpretable descriptor-based proxy that can be used to *prioritise*
candidates before expensive MD, and to frame hypotheses that MD would test.

Author: Sema Nur Bozdağ
"""

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen

DATA = os.path.join(os.path.dirname(__file__), "..", "data", "pfas_compounds.csv")
OUT = os.path.join(os.path.dirname(__file__), "..", "data", "pfas_descriptors.csv")


def compute_for_smiles(smiles: str) -> dict:
    """Return a dict of permeation-relevant descriptors, or raise on invalid SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")

    return {
        "mol_formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 2),
        "clogp": round(Crippen.MolLogP(mol), 3),          # lipophilicity proxy
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),  # polar surface area
        "hbd": rdMolDescriptors.CalcNumHBD(mol),           # H-bond donors
        "hba": rdMolDescriptors.CalcNumHBA(mol),           # H-bond acceptors
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "fraction_F": round(
            sum(a.GetSymbol() == "F" for a in mol.GetAtoms()) / mol.GetNumAtoms(), 3
        ),
    }


def main():
    df = pd.read_csv(DATA)
    records = []
    for _, row in df.iterrows():
        desc = compute_for_smiles(row["smiles"])
        records.append({**row.to_dict(), **desc})
    out = pd.DataFrame(records)

    # Simple, transparent permeation-propensity proxy.
    # Higher cLogP favours partitioning into the bilayer core; higher TPSA penalises
    # it (desolvation cost). This is a heuristic ranking score, NOT a permeability
    # coefficient. It exists to be falsifiable against MD-derived PMF later.
    out["permeation_proxy"] = (out["clogp"] - 0.03 * out["tpsa"]).round(3)

    out.to_csv(OUT, index=False)
    print(f"Wrote {len(out)} rows -> {OUT}")
    print(out[["abbreviation", "n_perfluoro_carbons", "headgroup",
               "clogp", "tpsa", "permeation_proxy"]].to_string(index=False))


if __name__ == "__main__":
    main()
