import subprocess

def run_gromacs_simulation(protein_name):
    # Example GROMACS commands using the protein_name
    subprocess.run(["gmx", "pdb2gmx", "-f", f"data/{protein_name}.pdb", "-o", f"data/{protein_name}_processed.gro", "-water", "spc"])
    subprocess.run(["gmx", "editconf", "-f", f"data/{protein_name}_processed.gro", "-o", f"data/{protein_name}_box.gro", "-c", "-d", "1.0", "-bt", "dodecahedron"])
    subprocess.run(["gmx", "solvate", "-cp", f"data/{protein_name}_box.gro", "-cs", "spc216.gro", "-o", f"data/{protein_name}_solv.gro"])
    subprocess.run(["gmx", "grompp", "-f", "gromacs/minim.mdp", "-c", f"data/{protein_name}_solv.gro", "-p", "data/topol.top", "-o", f"data/{protein_name}_minim.tpr"])
    subprocess.run(["gmx", "mdrun", "-v", "-deffnm", f"data/{protein_name}_minim"])