def coulombs_law(q1, q2, r):
    # Coulomb's law for electrostatic interaction
    k = 8.99e9  # Coulomb constant (N·m²·C⁻²)
    return k * (q1 * q2) / r**2
