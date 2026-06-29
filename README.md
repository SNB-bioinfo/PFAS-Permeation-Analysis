# PFAS Passive Membrane Permeation — A Descriptor-Based Exploratory Analysis

**Research question:** *Can classical molecular descriptors explain differences in passive membrane permeation propensity among per- and polyfluoroalkyl substances (PFAS)?*

```
Workflow

  PFAS SMILES (PubChem)
          |
          v
  RDKit descriptors  (cLogP, TPSA, MW, HBD/HBA, fluorine fraction)
          |
          v
  Exploratory analysis  (chain-length & headgroup trends, correlations)
          |
          v
  Hypothesis generation  (which compounds / factors to test)
          |
          v
  Future MD / AWH validation  (PMF of bilayer crossing)
```

## Why this repository

This repository was developed as a preparatory study before applying for a PhD focused on the *in silico* biophysical characterization of PFAS toxicokinetics and membrane permeation. Its goal is **not** to provide a validated predictive model, but to demonstrate the ability to explore chemical descriptors, formulate mechanistic hypotheses, and build a small, reproducible computational analysis around a concrete biophysical question. It is best read as a working note: a descriptor-level screen of the kind one might run to prioritise compounds *before* committing to expensive molecular dynamics (MD) simulations.

## Background and rationale

Passive permeation of a small molecule across a lipid bilayer is governed largely by two competing factors: lipophilicity, which favours partitioning into the hydrophobic core, and polar surface area / hydrogen-bonding capacity, which imposes a desolvation penalty at the headgroup. For PFAS, the perfluorinated tail drives lipophilicity while the ionisable headgroup (carboxylate or sulfonate) creates a strong desolvation penalty. One therefore expects permeation propensity to scale with fluorinated chain length, modulated by headgroup chemistry. This study tests whether simple RDKit-derived descriptors reproduce that expected trend across a homologous PFAS panel.

## Dataset

A panel of 10 PFAS compounds spanning two structural families and a range of chain lengths:

- **Perfluorocarboxylic acids (PFCAs):** PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFNA, PFDA (3–9 perfluorinated carbons)
- **Perfluorosulfonic acids (PFSAs):** PFBS, PFHxS, PFOS

SMILES strings were validated against known molecular formulas (10/10 match). Descriptors are computed with RDKit (`Crippen.MolLogP`, `CalcTPSA`, MW, H-bond donors/acceptors, rotatable bonds, fluorine fraction). See `data/pfas_descriptors.csv`.

### Data sources

The input data (`data/pfas_compounds.csv`) consists only of compound identities and **canonical SMILES structures for well-known PFAS, obtained from the PubChem database**. These are chemical structures, not measurements. **No experimental permeability, partition-coefficient, or toxicokinetic data are used anywhere in this repository** — every numeric value (cLogP, TPSA, MW, the permeation proxy, etc.) is computed *de novo* by RDKit from the SMILES. This distinction matters: the analysis demonstrates a descriptor-computation and hypothesis-framing workflow, not the fitting or validation of a model against measured data. Consequently, no supervised model fitting or predictive-performance evaluation (e.g. ROC/AUC) is performed or claimed.

## Methods

For each compound, descriptors relevant to passive permeation were computed from its SMILES (`src/compute_descriptors.py`). A deliberately simple, interpretable permeation proxy was defined as:

```
permeation_proxy = cLogP − 0.03 × TPSA
```

This is a heuristic, not a free-energy estimate. It encodes the two opposing physical contributions above — lipophilic partitioning (rewarded via cLogP) against the polar desolvation penalty (penalised via TPSA). The 0.03 weight is a rough scaling chosen so that the TPSA term meaningfully offsets cLogP across this panel, **not** a coefficient fitted to experimental permeability; the proxy is used only for relative ranking within this set. A true permeation free energy would require a potential of mean force (PMF) from MD with enhanced sampling (e.g. AWH), which is outside the scope of this descriptor-level screen.

## Results

**1. Chain length is the dominant lipophilicity driver.** Across the perfluorocarboxylic acid homologous series, cLogP increases linearly with the number of perfluorinated carbons at **+0.635 cLogP units per CF₂** (R² = 1.00 to four decimals, n = 7); across all 10 compounds the same slope holds (R² = 0.99). This is the expected monotonic trend — longer fluorinated tails are more lipophilic and, by this proxy, more permeation-prone (see Figure 1).

> *Caveat (see Limitations): the near-perfect fit largely reflects how RDKit's Crippen cLogP is constructed (an additive atom-contribution model), so this confirms internal consistency rather than independently validating the biophysics.*

**2. Headgroup chemistry shifts the trend.** At matched chain length, carboxylates are consistently slightly more lipophilic than sulfonates (e.g. at 6 perfluorinated carbons, PFHpA cLogP 3.81 vs PFHxS 3.57), and sulfonates carry a higher TPSA (54.4 vs 37.3 Å²) reflecting the additional sulfonate oxygens. The proxy captures this as a modest downward shift for the sulfonate family (see Figure 3).

**3. Descriptor structure.** The permeation proxy is dominated by lipophilicity and chain-length terms (correlation with cLogP r = 0.98, with fluorine fraction r = 0.96, with chain length r = 0.96) and only weakly anti-correlated with TPSA (r = −0.27) over this panel, because TPSA takes only two discrete values here (one per headgroup family). Figures are in `figures/`.

## Limitations

This is an exploratory descriptor screen with several deliberate simplifications that a follow-up study would need to address:

- **The proxy is heuristic, not free-energy based.** It ranks compounds; it does not estimate a permeability coefficient or a PMF barrier.
- **cLogP for fluorinated compounds is unreliable.** RDKit's Crippen cLogP is an atom-additive model parameterised mainly on non-fluorinated drug-like molecules; perfluorinated chains are poorly represented, so absolute cLogP values for PFAS should be treated with caution, and the strong chain-length linearity partly reflects the additive construction of the descriptor itself rather than independent evidence.
- **Ionisation is not modelled.** PFCAs and PFSAs are largely anionic at physiological pH; neutral-form descriptors omit the dominant electrostatic desolvation penalty that governs their real permeation and protein binding.
- **Small, homologous panel.** Ten compounds across two families cannot separate chain length from lipophilicity, which are nearly collinear by construction.
- **No experimental permeability for calibration or validation.**

## Lessons learned

- RDKit descriptors are useful for rapid hypothesis generation but insufficient on their own to characterise PFAS permeation, particularly given the limitations of additive lipophilicity models for fluorinated systems.
- Chain length and lipophilicity are nearly collinear in homologous PFAS series, so descriptor analysis alone cannot disentangle their separate contributions — precisely the kind of question MD/PMF calculations are needed to resolve.
- Descriptor-level screening is most valuable as a *prioritisation and hypothesis-framing* step upstream of molecular dynamics, rather than as a predictive endpoint in itself.

## Hypotheses for MD / enhanced-sampling validation

The descriptor trends above suggest concrete, testable hypotheses for a subsequent MD study:

1. PMF barriers for bilayer crossing decrease with perfluorinated chain length within each headgroup family.
2. At matched chain length, sulfonates show a higher crossing barrier than carboxylates, driven by greater headgroup desolvation cost.
3. Including the anionic (deprotonated) state will substantially raise computed barriers relative to the neutral-form descriptor ranking, potentially reordering the most/least permeant compounds.

## MD trajectory analysis demo

`src/md_trajectory_demo.py` is a short, self-contained demonstration that I can load and analyse MD trajectories with MDAnalysis (computing RMSD and radius of gyration). **It uses the bundled adenylate kinase example trajectory, not a PFAS–lipid-bilayer simulation**, and is included only to show familiarity with the MD analysis ecosystem — not as a result relevant to PFAS permeation. A genuine permeation study would replace this with a PFAS-in-membrane trajectory and a PMF/AWH workflow.

## Repository structure

```
data/        PFAS panel (SMILES) and computed descriptors
src/         descriptor computation, trend plots, MD analysis demo
figures/     generated figures
analysis.ipynb   notebook reproducing the descriptor analysis end to end
```

## Reproducibility

```bash
pip install -r requirements.txt
python src/compute_descriptors.py   # regenerate data/pfas_descriptors.csv
python src/plot_trends.py           # regenerate figures 01–03
python src/md_trajectory_demo.py    # regenerate figure 04 (MD demo)
```

Or open `analysis.ipynb` to run the descriptor analysis interactively.

## References

Methodological and background references for the descriptor choices and the PFAS panel:

- Lipinski, C. A., Lombardo, F., Dominy, B. W., & Feeney, P. J. (1997). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Advanced Drug Delivery Reviews*, 23(1–3), 3–25.
- Veber, D. F., Johnson, S. R., Cheng, H.-Y., Smith, B. R., Ward, K. W., & Kopple, K. D. (2002). Molecular properties that influence the oral bioavailability of drug candidates. *Journal of Medicinal Chemistry*, 45(12), 2615–2623.
- Wildman, S. A., & Crippen, G. M. (1999). Prediction of physicochemical parameters by atomic contributions. *Journal of Chemical Information and Computer Sciences*, 39(5), 868–873. *(basis of RDKit's cLogP)*
- Kim, S., Chen, J., Cheng, T., et al. (2023). PubChem 2023 update. *Nucleic Acids Research*, 51(D1), D1373–D1380. *(source of the SMILES structures)*
- Venable, R. M., Krämer, A., & Pastor, R. W. (2019). Molecular dynamics simulations of membrane permeability. *Chemical Reviews*, 119(9), 5954–5997. *(MD/PMF methodology for the proposed validation)*
- RDKit: Open-source cheminformatics. https://www.rdkit.org
- Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., & Beckstein, O. (2011). MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *Journal of Computational Chemistry*, 32(10), 2319–2327.

Recent work from the host group (Inserm U1248), to situate this preparatory study within the lab's membrane-transport modelling programme:

- Crespi, V., Tóth, Á., Janaszkiewicz, A., Falguières, T., & Di Meo, F. (2025). Membrane-dependent dynamics and dual translocation mechanisms of ABCB4: insights from molecular dynamics simulations. *Computational and Structural Biotechnology Journal*, 27, 1215–1232.
- Tóth, Á., Janaszkiewicz, A., Crespi, V., & Di Meo, F. (2023). On the interplay between lipids and asymmetric dynamics of an NBS degenerate ABC transporter. *Communications Biology*, 6, 149.

## License

Released under the MIT License — see `LICENSE`.
