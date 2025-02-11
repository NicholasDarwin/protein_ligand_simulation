import openmm
from openmm import app
from openmm.unit import *
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_ligand_parameters(mol2_path):
    """Generate parameters for the ligand using RDKit"""
    # Read mol2 file
    mol = Chem.MolFromMol2File(mol2_path, removeHs=False)
    if mol is None:
        raise ValueError("Failed to read mol2 file")
    
    # Add hydrogens if needed
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates if not present
    if not mol.GetConformer().Is3D():
        AllChem.EmbedMolecule(mol)
    
    # Optimize the geometry
    AllChem.MMFFOptimizeMolecule(mol)
    
    return mol

def create_system_with_custom_ligand(pdb, ligand_mol, forcefield):
    """Create OpenMM system with custom ligand parameters"""
    # Create a new topology that includes the ligand
    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    # Add ligand parameters to the force field
    # This is a simplified version - in practice, you'd need more sophisticated parameter generation
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        rigidWater=True,
        removeCMMotion=True
    )
    
    return system, modeller

def run_simulation(input_file, force_field="amber99sbildn.xml", temperature=300*kelvin, time_step=2*femtoseconds, steps=50000):
    try:
        # Handle uploaded file
        if hasattr(input_file, 'name'):
            with tempfile.NamedTemporaryFile(suffix='.mol2', delete=False) as tmp_file:
                tmp_file.write(input_file.getvalue())
                mol2_path = tmp_file.name
        else:
            mol2_path = input_file

        # Generate ligand parameters
        ligand_mol = generate_ligand_parameters(mol2_path)
        
        # Convert to PDB format for OpenMM
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as pdb_file:
            writer = Chem.PDBWriter(pdb_file.name)
            writer.write(ligand_mol)
            writer.close()
            pdb = app.PDBFile(pdb_file.name)

        # Load force field
        try:
            forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')
        except Exception as e:
            raise ValueError(f"Failed to load force field: {str(e)}")

        # Create system with custom ligand
        system, modeller = create_system_with_custom_ligand(pdb, ligand_mol, forcefield)

        # Set up integrator
        integrator = openmm.LangevinIntegrator(temperature, 1/picosecond, time_step)
        
        # Create simulation
        simulation = app.Simulation(modeller.topology, system, integrator)
        
        # Set initial positions
        simulation.context.setPositions(modeller.positions)
        
        # Minimize energy
        print("Minimizing energy...")
        simulation.minimizeEnergy()

        # Set up reporters
        simulation.reporters.append(app.DCDReporter('trajectory.dcd', 1000))
        simulation.reporters.append(app.StateDataReporter(
            'output.txt', 
            1000, 
            step=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=steps,
            separator='\t'
        ))

        # Run simulation
        print(f"Running simulation for {steps} steps...")
        simulation.step(steps)
        print("Simulation completed successfully")
        
        return simulation

    except Exception as e:
        raise Exception(f"Simulation failed: {str(e)}")
    finally:
        # Clean up temporary files
        if 'mol2_path' in locals():
            try:
                os.unlink(mol2_path)
            except:
                pass
        if 'pdb_file' in locals():
            try:
                os.unlink(pdb_file.name)
            except:
                pass