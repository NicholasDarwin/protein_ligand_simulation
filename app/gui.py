import tkinter as tk
from tkinter import messagebox, scrolledtext
from app.pdb_fetcher import fetch_pdb
from app.simulation import run_gromacs_simulation

def launch_gui():
    root = tk.Tk()
    root.title("Protein-Ligand MD Simulation")
    
    # Create and pack widgets
    tk.Label(root, text="Protein Name:").pack()
    global protein_name_entry
    protein_name_entry = tk.Entry(root)
    protein_name_entry.pack()
    
    tk.Label(root, text="Ligand Name:").pack()
    global ligand_name_entry
    ligand_name_entry = tk.Entry(root)
    ligand_name_entry.pack()
    
    # Add text area for output
    global output_text
    output_text = scrolledtext.ScrolledText(root, width=80, height=20)
    output_text.pack(padx=10, pady=10)
    
    tk.Button(root, text="Fetch PDB", command=on_search_button_click).pack()
    tk.Button(root, text="Run MD Simulation", command=on_run_simulation_click).pack()
    
    root.mainloop()

def on_run_simulation_click():
    protein_name = protein_name_entry.get()
    ligand_name = ligand_name_entry.get()
    if protein_name and ligand_name:
        # Trigger GROMACS simulation with both protein and ligand names
        run_gromacs_simulation(protein_name, ligand_name)
        messagebox.showinfo("Success", "Simulation started!")
    else:
        messagebox.showerror("Error", "Please enter both protein and ligand names.")