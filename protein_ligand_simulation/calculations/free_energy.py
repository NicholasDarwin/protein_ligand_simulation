def free_energy(enthalpy, entropy, temperature):
    # Calculate free energy change
    return enthalpy - temperature * entropy
