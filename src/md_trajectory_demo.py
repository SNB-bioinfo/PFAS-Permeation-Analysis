"""
md_trajectory_demo.py
---------------------
A minimal, honest demonstration that I can read and analyse molecular dynamics
output using the MDAnalysis ecosystem. This does NOT run an MD simulation — it
analyses an existing trajectory bundled with MDAnalysis.

Why this system is relevant to the PhD project
----------------------------------------------
The example trajectory is adenylate kinase (AdK), a protein that undergoes a
large open <-> closed domain motion. This is a useful proxy for the conformational
problem at the heart of membrane-transporter modelling: transporters such as
OAT1/OAT3 or P-gp cycle between inward-facing and outward-facing states, and
characterising those transitions (and how a ligand like a PFAS shifts the
equilibrium) is exactly what unbiased and enhanced-sampling MD is used for.

What this script does
---------------------
1. Loads a topology + trajectory.
2. Computes backbone RMSD over time relative to the first frame
   (a standard first check of conformational drift / sampling).
3. Computes the radius of gyration as a coarse measure of compaction
   (open vs closed state).
4. Saves a figure summarising both.

Author: Sema Nur Bozdağ
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.tests.datafiles import PSF, DCD

HERE = os.path.dirname(__file__)
FIGDIR = os.path.join(HERE, "..", "figures")
os.makedirs(FIGDIR, exist_ok=True)


def main():
    u = mda.Universe(PSF, DCD)
    print(f"Loaded system: {len(u.atoms)} atoms, {len(u.trajectory)} frames")

    # --- backbone RMSD relative to first frame ---
    R = rms.RMSD(u, u, select="backbone", ref_frame=0)
    R.run()
    # results.rmsd columns: [frame, time, RMSD]
    rmsd = R.results.rmsd[:, 2]

    # --- radius of gyration per frame ---
    backbone = u.select_atoms("backbone")
    rgyr = []
    for _ in u.trajectory:
        rgyr.append(backbone.radius_of_gyration())
    rgyr = np.array(rgyr)

    frames = np.arange(len(rmsd))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(frames, rmsd, color="#2a6f97")
    ax1.set_xlabel("Frame")
    ax1.set_ylabel("Backbone RMSD (Å)")
    ax1.set_title("Conformational drift from frame 0")
    ax1.grid(alpha=0.3)

    ax2.plot(frames, rgyr, color="#bc4749")
    ax2.set_xlabel("Frame")
    ax2.set_ylabel("Radius of gyration (Å)")
    ax2.set_title("Compaction (open ↔ closed motion)")
    ax2.grid(alpha=0.3)

    fig.suptitle("AdK trajectory analysis — proxy for transporter "
                 "inward/outward-facing transitions")
    fig.tight_layout()
    out = os.path.join(FIGDIR, "04_md_trajectory_analysis.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)

    print(f"RMSD range: {rmsd.min():.2f}–{rmsd.max():.2f} Å")
    print(f"Rgyr range: {rgyr.min():.2f}–{rgyr.max():.2f} Å")
    print(f"Saved figure -> {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
