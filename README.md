# Virtual Ligand Screening Pipeline

## Description
This Python-based program simulates molecular docking, ligand-protein interactions, and drug binding analysis. The program allows virtual screening of ligands and includes several advanced calculations such as Monte Carlo docking, genetic algorithms, free energy calculations, and molecular dynamics.

Virtual Ligand Screening Pipeline

Objective: Create a Python-based program for simulating molecular docking, ligand-protein interactions, and drug binding analysis. The calculations should be performed manually, with OpenMM utilized solely for viewing the results (e.g., 3D visualization and molecular dynamics simulations). This system will include a comprehensive set of molecular docking, force-field calculations, free energy calculations, and relevant bioinformatics functions for the virtual screening of ligands.

File Structure and Calculations:
The pipeline should be modular, with each file dedicated to a specific calculation or utility. Below is a list of the necessary files, their calculations, and their purposes.

1. Calculations Directory (calculations/)
monte_carlo_docking.py: Implements Monte Carlo sampling to explore the ligand conformational space and perform docking based on energy minimization.
genetic_algorithm.py: Applies genetic algorithms for optimizing docking poses and ligand conformations.
poisson_boltzmann.py: Implements the Poisson-Boltzmann equation for calculating electrostatic solvation effects.
generalized_born.py: Applies the Generalized Born model for estimating solvation energy.
normal_mode_analysis.py: Implements Normal Mode Analysis for calculating entropy contributions to the system.
free_energy_perturbation.py: Computes free energy changes based on perturbation theory.
docking_score.py: Evaluates the docking pose using scoring functions such as van der Waals, electrostatic, hydrogen bonding, and torsion energy.
newtons_law.py: Simulates molecular motion based on Newton's laws (F=ma).
lennard_jones.py: Models van der Waals forces using the Lennard-Jones potential.
coulombs_law.py: Models electrostatic interactions using Coulomb's law.
free_energy.py: Calculates binding free energy (ΔG = ΔH - TΔS) for predicting drug binding strength.
molecular_dynamics.py: Models realistic molecular movement using force-field interactions.
sasa.py: Calculates the Solvent Accessible Surface Area (SASA) for solvation effects.
rmsd.py: Computes Root Mean Square Deviation (RMSD) to evaluate the conformational stability of ligand-protein complexes.
radius_of_gyration.py: Computes the Radius of Gyration (Rg) to assess the compactness of molecular structures.
hydrogen_bond_energy.py: Estimates hydrogen bond energies based on distance and interaction strength.
mm_pbsa.py: Implements MM-PBSA for calculating binding affinity using Molecular Mechanics combined with the Poisson-Boltzmann equation.
qm_mm.py: Implements Quantum Mechanics/Molecular Mechanics (QM/MM) for more accurate docking calculations.
hydration_free_energy.py: Calculates the contribution of water molecules in docking (hydration free energy).
2. UI Directory (ui/)
dashboard.py: The main user interface for project selection, controlling simulations, and displaying results.
molecular_input.py: Provides an interface for manual molecular input, including atom and bond definition.
ligand_screening.py: Manages ligand library upload and batch screening results.
molecular_visualization.py: Displays 3D molecular structures, protein-ligand interactions, and results.
analysis_export.py: Manages data export, including binding scores and interaction diagrams.
3. Visualization Directory (visualization/)
py3dmol_viewer.py: Implements the Py3Dmol library for real-time molecular visualization.
rdkit_tools.py: Provides RDKit utilities for molecule rendering and structure manipulation.
ngl_viewer.py: Interfaces with the NGL Viewer for web-based molecular visualization.
pymol_tools.py: Integrates PyMOL for detailed structural analysis and rendering.
heatmap.py: Generates heatmaps for energy distributions and molecular interactions.
interaction_map.py: Visualizes interaction maps (e.g., hydrogen bonds, electrostatic forces).
4. Docking Directory (docking/)
docking_engine.py: Main docking engine for running ligand-protein docking simulations.
scoring_functions.py: Implements various scoring functions for evaluating docking poses.
search_algorithms.py: Contains algorithms for efficient ligand search, including Monte Carlo, genetic algorithms, and others.
5. Utilities Directory (utils/)
file_io.py: Utility for reading and writing molecular files (e.g., .sdf, .mol2, .csv, .pdb).
constants.py: Contains constants and unit conversions for calculations.
helpers.py: General helper functions used across multiple modules.
6. Main Application (app.py)
Entry point to launch the program, set up the UI, and integrate all modules.
7. Supporting Files
requirements.txt: Dependencies required for the project, such as RDKit, Py3Dmol, Matplotlib, OpenMM, etc.
README.md: Project documentation and setup instructions.
.gitignore: To ignore unnecessary files and directories.
Key Requirements:

The system must be flexible to accommodate new scoring functions and docking algorithms.
Provide clear user interfaces for inputting molecular structures, viewing results, and exporting data.
Enable the integration of OpenMM only for visualizing results from the calculations.
Maintain modularity so that each calculation (e.g., Monte Carlo docking, genetic algorithms, etc.) is implemented in its own file.
Optimize for efficient ligand screening and binding affinity prediction.


virtual_ligand_screening/
│── calculations/
│   │── __init__.py
│   │── monte_carlo_docking.py
│   │── genetic_algorithm.py
│   │── poisson_boltzmann.py
│   │── generalized_born.py
│   │── normal_mode_analysis.py
│   │── free_energy_perturbation.py
│   │── docking_score.py
│   │── newtons_law.py
│   │── lennard_jones.py
│   │── coulombs_law.py
│   │── free_energy.py
│   │── molecular_dynamics.py
│   │── sasa.py
│   │── rmsd.py
│   │── radius_of_gyration.py
│   │── hydrogen_bond_energy.py
│   │── mm_pbsa.py
│   │── qm_mm.py
│   │── hydration_free_energy.py
│
│── ui/
│   │── __init__.py
│   │── dashboard.py
│   │── molecular_input.py
│   │── ligand_screening.py
│   │── molecular_visualization.py
│   │── analysis_export.py
│
│── visualization/
│   │── __init__.py
│   │── py3dmol_viewer.py
│   │── rdkit_tools.py
│   │── ngl_viewer.py
│   │── pymol_tools.py
│   │── heatmap.py
│   │── interaction_map.py
│
│── docking/
│   │── __init__.py
│   │── docking_engine.py
│   │── scoring_functions.py
│   │── search_algorithms.py
│
│── utils/
│   │── __init__.py
│   │── file_io.py
│   │── constants.py
│   │── helpers.py
│
│── app.py  # Main entry point
│── requirements.txt  # Dependencies (RDKit, Py3Dmol, Matplotlib, etc.)
│── README.md  # Project documentation
│── .gitignore  # Ignore unnecessary files

