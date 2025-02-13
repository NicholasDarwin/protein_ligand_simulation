import subprocess

def run_gromacs_simulation(protein_name, ligand_name):
    # Process protein and ligand
    subprocess.run(["gmx", "pdb2gmx", "-f", f"data/{protein_name}.pdb", "-o", f"data/{protein_name}_processed.gro", "-water", "spc"])
    subprocess.run(["gmx", "pdb2gmx", "-f", f"data/{ligand_name}.pdb", "-o", f"data/{ligand_name}_processed.gro", "-water", "spc"])
    
    # Set up simulation boxes
    subprocess.run(["gmx", "editconf", "-f", f"data/{protein_name}_processed.gro", "-o", f"data/{protein_name}_box.gro", "-c", "-d", "1.0", "-bt", "dodecahedron"])
    subprocess.run(["gmx", "editconf", "-f", f"data/{ligand_name}_processed.gro", "-o", f"data/{ligand_name}_box.gro", "-c", "-d", "1.0", "-bt", "dodecahedron"])
    
    # Merge protein and ligand
    with open("data/merged.gro", "w") as outfile:
        for filename in [f"data/{protein_name}_box.gro", f"data/{ligand_name}_box.gro"]:
            with open(filename) as infile:
                outfile.write(infile.read())
    
    # Continue with solvation and minimization
    subprocess.run(["gmx", "solvate", "-cp", "data/merged.gro", "-cs", "spc216.gro", "-o", "data/system_solv.gro"])
    subprocess.run(["gmx", "grompp", "-f", "gromacs/minim.mdp", "-c", "data/system_solv.gro", "-p", "data/topol.top", "-o", "data/system_minim.tpr"])
    subprocess.run(["gmx", "mdrun", "-v", "-deffnm", "data/system_minim"])