import tkinter as tk
from tkinter import messagebox
from app.pdb_fetcher import fetch_pdb
from app.simulation import run_gromacs_simulation

def on_search_button_click():
    protein_name = protein_name_entry.get()
    ligand_name = ligand_name_entry.get()
    if protein_name and ligand_name:
        # Fetch PDB for protein
        pdb_file = fetch_pdb(protein_name)
        if pdb_file:
            messagebox.showinfo("Success", f"PDB file for {protein_name} downloaded successfully.")
        else:
            messagebox.showerror("Error", "Failed to fetch the PDB file.")
    else:
        messagebox.showerror("Error", "Please enter both protein and ligand names.")

def on_run_simulation_click():
    protein_name = protein_name_entry.get()
    if protein_name:
        # Trigger GROMACS simulation with the protein name
        run_gromacs_simulation(protein_name)
        messagebox.showinfo("Success", "Simulation started!")
    else:
        messagebox.showerror("Error", "Please enter a protein name first.")

def launch_gui():
    root = tk.Tk()
    root.title("Protein-Ligand MD Simulation")

    # Protein Name Entry
    protein_name_label = tk.Label(root, text="Protein Name:")
    protein_name_label.grid(row=0, column=0)
    global protein_name_entry
    protein_name_entry = tk.Entry(root)
    protein_name_entry.grid(row=0, column=1)

    # Ligand Name Entry
    ligand_name_label = tk.Label(root, text="Ligand Name:")
    ligand_name_label.grid(row=1, column=0)
    global ligand_name_entry
    ligand_name_entry = tk.Entry(root)
    ligand_name_entry.grid(row=1, column=1)

    # Search Button
    search_button = tk.Button(root, text="Fetch PDB", command=on_search_button_click)
    search_button.grid(row=2, columnspan=2)

    # Run Simulation Button
    run_simulation_button = tk.Button(root, text="Run MD Simulation", command=on_run_simulation_click)
    run_simulation_button.grid(row=3, columnspan=2)

    root.mainloop()
