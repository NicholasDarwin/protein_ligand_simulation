# visualization/py3dmol_viewer.py
import py3Dmol
import streamlit.components.v1 as components

def visualize_structure(mol):
    """Visualize molecule using py3Dmol"""
    view = py3Dmol.view(width=800, height=600)
    
    # Convert molecule to PDB format for visualization
    pdb = Chem.MolToPDBBlock(mol)
    view.addModel(pdb, "pdb")
    
    # Style the visualization
    view.setStyle({'stick':{}})
    view.zoomTo()
    
    # Display in Streamlit
    components.html(view._repr_html_(), height=600)