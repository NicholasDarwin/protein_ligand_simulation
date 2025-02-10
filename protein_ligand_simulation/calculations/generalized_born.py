def generalized_born(ligand, receptor):
    # Calculate Born radii
    ligand_radii = calculate_born_radii(ligand)
    receptor_radii = calculate_born_radii(receptor)
    
    # Calculate GB energy
    epsilon_solvent = 78.5  # Water dielectric constant
    epsilon_vacuum = 1.0
    
    energy = sum(calculate_gb_term(q1, q2, r, r1, r2, epsilon_solvent, epsilon_vacuum)
                for q1, q2, r, r1, r2 in get_atom_pairs_with_radii(ligand, receptor,
                                                                  ligand_radii, receptor_radii))
    return energy