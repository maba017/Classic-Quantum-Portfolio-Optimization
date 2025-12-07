import numpy as np
import pandas as pd
import time
from itertools import chain
from functools import reduce
import quimb as qu
from quimb import ikron, pauli

class IsingEvolver:
    """
    Motor de Dinámica Cuántica (Optimizado).
    Calcula el Gap Relativo: (E1 - E0) / (Emax - E0).
    """

    def __init__(self, h, J, num_qubits, g=2.0, k=4):
        self.h = h
        self.J = J
        self.n = num_qubits
        self.g = g
        self.k = k
        self.dim = 2**self.n

        print(f"Motor Quimb (Gap Relativo): {self.n} qubits")

        # Operadores
        self.Sx = pauli("X", sparse=True).astype(np.float32)
        self.Sz = pauli("Z", sparse=True).astype(np.float32)
        
        self.H_z = self._get_H_z()
        self.H_x = self._get_H_x()

    def _get_H_z(self):
        dims = [2] * self.n
        partsZ = (val * ikron(self.Sz, dims, inds=[i]) for i, val in self.h.items())
        partsZZ = (val * ikron(self.Sz, dims, inds=list(pair)) for pair, val in self.J.items())
        parts = chain(partsZ, partsZZ)
        try:
            H = reduce(lambda x, y: x + y, parts)
        except TypeError:
            H = qu.sparse.eye(self.dim) * 0.0
        return H.astype(np.float32)

    def _get_H_x(self):
        dims = [2] * self.n
        if abs(self.g) > 0:
            partsX = (-self.g * ikron(self.Sx, dims, inds=[i]) for i in range(self.n))
            H = reduce(lambda x, y: x + y, partsX)
        else:
            H = qu.sparse.eye(self.dim) * 0.0
        return H.astype(np.float32)

    def _get_schedule(self, steps):
        """Schedule logarítmico invertido (zoom al final)."""
        _start = 1e-8
        sched = _start + 1 - np.geomspace(_start, 1.0, num=steps)
        sched = sched[::-1]
        
        z_sched = 1.0 * sched
        z_sched[0] = 0.0
        x_sched = 1.0 * (1 - sched)
        x_sched[0] = 1.0
        
        return sched, z_sched, x_sched

    def analyze_adiabatic_path(self, steps=50):
        t_values, z_sched, x_sched = self._get_schedule(steps)
        results = []
        
        print(f"Simulando trayectoria (Gap Relativo, {steps} pasos)...")
        start_t = time.time()

        # Warm start solo para el estado fundamental (el crítico)
        v0 = None

        for i in range(len(t_values)):
            H_s = z_sched[i] * self.H_z + x_sched[i] * self.H_x
            
            try:
                # 1. Buscar autovalores bajos (E0, E1...) con Warm Start
                eig_vals, eig_vecs = qu.eigensystem_partial(
                    H_s, k=self.k, isherm=True, which="SA", 
                    tol=1e-3, v0=v0, return_vecs=True
                )
                
                # 2. Buscar autovalor máximo (Emax)
                # Usamos 'LA' (Largest Algebraic). No hace falta warm start, converge rápido.
                eig_max, _ = qu.eigensystem_partial(
                    H_s, k=1, isherm=True, which="LA", tol=1e-3
                )
                
                # Actualizar warm start
                v0 = eig_vecs[:, 0]

                # --- CÁLCULOS FÍSICOS ---
                E0 = eig_vals[0]
                E1 = eig_vals[1]
                Emax = eig_max[0]
                
                # Fórmula del Gap Relativo
                gap_abs = E1 - E0
                width = Emax - E0
                gap_rel = gap_abs / width if width != 0 else 0.0

                # Decodificar solución (recordando invertir bits si es necesario)
                best_idx = np.argmax(np.abs(v0)**2)
                raw_bits = [int(b) for b in format(best_idx, f'0{self.n}b')]
                # Invertimos bits (0->1, 1->0) para alinear con D-Wave si es necesario
                # Si en tu nb original no invertías, quita el "1-b"
                corrected_bits = [1 - b for b in raw_bits]

                row = {
                    't': t_values[i],
                    'gap_rel': gap_rel,     # (E1-E0)/(Emax-E0)
                    'gap_abs': gap_abs,     # E1-E0
                    'E0': E0,
                    'Emax': Emax,
                    'best_bits': corrected_bits
                }
                
                for j in range(self.k):
                    row[f'E{j}'] = eig_vals[j]
                
                results.append(row)

            except Exception as e:
                print(f"Error paso {i}: {e}")
                if results: results.append(results[-1])

        print(f"Completado en {time.time() - start_t:.2f} s")
        return pd.DataFrame(results)