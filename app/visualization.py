# In visualization.py
import pymol

def visualize_structure(pdb_file):
    pymol.finish_launching()
    pymol.cmd.load(pdb_file)
    pymol.cmd.show("cartoon")
    pymol.cmd.color("blue")
    pymol.cmd.show("sticks", "hetatm")