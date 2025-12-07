import numpy as np
import pandas as pd
import time
from itertools import chain
from functools import reduce
import quimb as qu
from quimb import ikron, pauli
import scipy.sparse.linalg as sla

class IsingEvolver:
    """
    Motor de Dinámica Cuántica (Schedule Lineal).
    Evoluciona H(t) linealmente de 0 a 1.
    """

    def __init__(self, h, J, num_qubits, g=2.0, k=4):
        self.h = h
        self.J = J
        self.n = num_qubits
        self.g = g
        self.k = k
        self.dim = 2**self.n
        
        print(f"⚛️ Motor Quimb Lineal: {self.n} qubits (Dim: {self.dim})")
        
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
        """
        Schedule estrictamente LINEAL de 0 a 1.
        s(t) = t/T
        """
        # np.linspace genera puntos equidistantes incluyendo 0 y 1
        sched = np.linspace(0.0, 1.0, num=steps)
        
        z_sched = sched          # Sube de 0 a 1
        x_sched = 1.0 - sched    # Baja de 1 a 0
        
        return sched, z_sched, x_sched

    def analyze_adiabatic_path(self, steps=50):
        t_values, z_sched, x_sched = self._get_schedule(steps)
        results = []
        
        print(f"🌊 Simulando trayectoria lineal ({steps} pasos)...")
        start_t = time.time()

        # Warm start vector (Inicializamos con superposición uniforme que es el GS de Hx)
        v0 = np.ones((self.dim, 1)) / np.sqrt(self.dim)
        v0 = v0.astype(np.float32)

        for i in range(len(t_values)):
            H_s = z_sched[i] * self.H_z + x_sched[i] * self.H_x
            
            try:
                # Diagonalización parcial (k primeros)
                vals, vecs = sla.eigsh(
                    H_s, k=self.k, which='SA', v0=v0, tol=1e-3, return_eigenvectors=True
                )
                
                # Ordenar
                idx = np.argsort(vals)
                vals = vals[idx]
                vecs = vecs[:, idx]
                
                # Autovalor máximo (para gap relativo)
                vals_max, _ = sla.eigsh(
                    H_s, k=1, which='LA', tol=1e-3, return_eigenvectors=False
                )
                
                # Actualizar warm start
                v0 = vecs[:, 0]

                # Datos
                E0 = vals[0]
                E1 = vals[1]
                Emax = vals_max[0]
                
                gap_abs = E1 - E0
                width = Emax - E0
                gap_rel = gap_abs / width if width > 1e-9 else 0.0

                # Decodificar
                best_idx = np.argmax(np.abs(v0)**2)
                raw_bits = [int(b) for b in format(best_idx, f'0{self.n}b')]
                
                # Inversión de bits si es necesaria (depende de tu modelo QUBO)
                corrected_bits = [1 - b for b in raw_bits] 
                corrected_bits = raw_bits

                row = {
                    't': t_values[i],
                    'gap_rel': gap_rel,     
                    'gap_abs': gap_abs,     
                    'E0': E0,
                    'Emax': Emax,
                    'best_bits': corrected_bits
                }
                
                for j in range(self.k):
                    row[f'E{j}'] = vals[j]
                
                results.append(row)

            except Exception as e:
                print(f"⚠️ Error paso {i}: {e}")
                if results: results.append(results[-1])

        print(f"✅ Completado en {time.time() - start_t:.2f} s")
        return pd.DataFrame(results)