# visualization/py3dmol_viewer.py
import py3Dmol

def visualize_structure(pdb_file_path):
    with open(pdb_file_path, 'r') as file:
        pdb_data = file.read()
    
    # Create a viewer
    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(pdb_data, "pdb")
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.show()
