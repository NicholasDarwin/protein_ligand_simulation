from rdkit import Chem
from rdkit.Chem import Draw

def render_molecule(molecule_data):
    mol = Chem.MolFromMolBlock(molecule_data)
    img = Draw.MolToImage(mol)
    return img
