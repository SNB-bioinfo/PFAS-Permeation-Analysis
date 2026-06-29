"""
compute_descriptors.py
----------------------
Computes physicochemical descriptors for a panel of PFAS compounds using RDKit,
with a focus on properties relevant to passive membrane permeation.

Scientific rationale
--------------------
Passive permeation of a small molecule across a lipid bilayer is governed largely
by its lipophilicity (partitioning into the hydrophobic core) and its polar surface
area / hydrogen-bonding capacity (desolvation penalty at the headgroup). Topological
polar surface area (TPSA) and lipophilicity (logP) are classical correlates of
passive permeability in medicinal chemistry (Lipinski et al., 1997; Veber et al.,
2002). For PFAS, the perfluorinated tail drives lipophilicity while the ionisable
headgroup (carboxylate / sulfonate) creates a strong desolvation penalty. We
therefore expect permeation propensity to scale with fluorinated chain length,
modulated by headgroup chemistry.

This script does NOT compute a true free-energy profile (PMF) — that requires
molecular dynamics with enhanced sampling (e.g. the AWH protocol). It produces a
fast, interpretable descriptor-based proxy that can be used to *prioritise*
candidates before expensive MD, and to frame hypotheses that MD would test. The
descriptors generated here are intended for exploratory analysis and hypothesis
generation, not quantitative prediction.

Data source
-----------
Input SMILES (data/pfas_compounds.csv) are canonical structures for well-known
PFAS, taken from PubChem. No experimental permeability or toxicokinetic data are
used; every numeric descriptor below is computed de novo by RDKit.

Author: Sema Nur Bozdağ
"""

from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data" / "pfas_compounds.csv"
OUT = BASE / "data" / "pfas_descriptors.csv"


def compute_for_smiles(smiles: str) -> dict:
    """Return a dict of permeation-relevant descriptors, or raise on invalid SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")

    return {
        "mol_formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 2),
        "clogp": round(Crippen.MolLogP(mol), 3),          # lipophilicity (Crippen)
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),  # polar surface area
        "hbd": rdMolDescriptors.CalcNumHBD(mol),           # H-bond donors
        "hba": rdMolDescriptors.CalcNumHBA(mol),           # H-bond acceptors
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "fraction_F": round(
            sum(a.GetSymbol() == "F" for a in mol.GetAtoms()) / mol.GetNumAtoms(), 3
        ),
    }


def main() -> None:
    df = pd.read_csv(DATA)
    records = [
        {**row.to_dict(), **compute_for_smiles(row["smiles"])}
        for _, row in df.iterrows()
    ]
    out = pd.DataFrame(records)

    # Simple, transparent permeation-propensity proxy.
    # Higher cLogP favours partitioning into the bilayer core; higher TPSA penalises
    # it (desolvation cost), consistent with TPSA being a long-standing correlate of
    # passive permeability in medicinal chemistry (Veber et al., 2002). This is a
    # heuristic ranking score, NOT a permeability coefficient, and exists to be
    # falsifiable against MD-derived PMF later.
    #
    # The 0.03 weight on TPSA is a rough scaling, NOT a fitted coefficient: it is
    # chosen only so the polar term offsets cLogP on a comparable scale across this
    # panel, for relative ranking within the set. Future work may replace this
    # heuristic weighting with coefficients derived from experimental permeability
    # or MD-derived free-energy barriers.
    #
    # CAVEAT: RDKit's Crippen cLogP is an additive atom-contribution model
    # parameterised mainly on non-fluorinated drug-like molecules. The values here
    # should not be interpreted as experimentally accurate lipophilicity values for
    # PFAS, only as relative trends within this panel.
    out["permeation_proxy"] = (out["clogp"] - 0.03 * out["tpsa"]).round(3)

    out.to_csv(OUT, index=False)
    print(f"Wrote {len(out)} rows -> {OUT}")
    print(out[["abbreviation", "n_perfluoro_carbons", "headgroup",
               "clogp", "tpsa", "permeation_proxy"]].to_string(index=False))


if __name__ == "__main__":
    main()
