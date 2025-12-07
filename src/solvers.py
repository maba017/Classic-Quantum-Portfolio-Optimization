import dimod
from dwave.samplers import SimulatedAnnealingSampler
from dimod import ExactSolver
import time

class PortfolioSolver:
    """
    Clase para resolver el problema de optimización de carteras usando
    samplers clásicos de D-Wave.
    """

    def __init__(self, h, J):
        self.h = h
        self.J = J
        self.bqm = dimod.BinaryQuadraticModel.from_ising(h, J)

    def solve_simulated_annealing(self, num_reads=1000, num_sweeps=1000):
        """
        Resuelve usando Simulated Annealing.
        Devuelve: (mejor_solucion, mejor_energia, sampleset_completo)
        """
        print(f"Ejecutando Simulated Annealing Clásico ({num_reads} muestras)...")
        start_time = time.time()
        
        sampler = SimulatedAnnealingSampler()
        sampleset = sampler.sample(
            self.bqm, 
            num_reads=num_reads, 
            num_sweeps=num_sweeps
        )
        
        elapsed = time.time() - start_time
        print(f"Terminado en {elapsed:.4f} segundos.")
        
        best_solution = sampleset.first.sample
        best_energy = sampleset.first.energy
        
        # AHORA DEVOLVEMOS TAMBIÉN EL 'sampleset' COMPLETO
        return best_solution, best_energy, sampleset

    def solve_exact(self):
        """
        Resuelve usando Fuerza Bruta (ExactSolver).
        """
        print("🕵️ Ejecutando Exact Solver (Validación)...")
        if len(self.bqm.variables) > 20:
            print("⚠️ ADVERTENCIA: Demasiadas variables. Se puede bloquear.")
            
        sampler = ExactSolver()
        sampleset = sampler.sample(self.bqm)
        
        best_solution = sampleset.first.sample
        best_energy = sampleset.first.energy
        
        return best_solution, best_energy, sampleset