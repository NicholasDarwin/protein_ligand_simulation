def poisson_boltzmann(ligand, receptor):
    # Constants
    epsilon_solvent = 78.5  # Water dielectric constant
    epsilon_protein = 4.0   # Protein dielectric constant
    kappa = 0.1            # Debye-Hückel parameter
    
    # Calculate potential using linearized PB equation
    phi = solve_pb_equation(ligand, receptor, epsilon_solvent, 
                          epsilon_protein, kappa)
    
    # Calculate solvation energy
    energy = 0.5 * sum(q * p for q, p in zip(get_charges(ligand, receptor), phi))
    return energy