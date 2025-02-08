import unittest
from src.virus_simulation import lennard_jones_potential, simulate_virus_dynamics

class TestVirusSimulation(unittest.TestCase):
    
    def test_lennard_jones_potential(self):
        # Test the Lennard-Jones potential calculation
        result = lennard_jones_potential(1.0)
        self.assertAlmostEqual(result, 0.0, places=1)  # Expected value for distance = 1
    
    def test_simulate_virus_dynamics(self):
        # Test the virus simulation logic
        positions = simulate_virus_dynamics()
        self.assertEqual(positions.shape[1], 3)  # Ensure 3D positions
        self.assertEqual(positions.shape[0], 10)  # Should have 10 particles

if __name__ == '__main__':
    unittest.main()
