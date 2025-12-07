import numpy as np
import math
import dimod

class QuboFormulation:
    """
    Clase para convertir el problema de optimización de carteras (Markowitz)
    en un modelo QUBO/Ising compatible con D-Wave.
    
    Basado en la formulación de:
    Burges Bruna, G. (2025). Computación Cuántica Adiabática para Finanzas.
    """

    def __init__(self, mu, sigma, budget, stock_prices, risk_aversion=1.0, budget_penalty=1.0):
        """
        Args:
            mu (list): Retornos esperados de cada activo.
            sigma (list of lists): Matriz de covarianzas.
            budget (float): Presupuesto total.
            stock_prices (list): Precios actuales.
            risk_aversion (float): Alpha.
            budget_penalty (float): Lambda.
        """
        self.mu = mu
        self.sigma = sigma
        self.budget = budget
        self.prices = stock_prices
        self.alpha = risk_aversion
        self.lam = budget_penalty
        
        self.num_assets = len(mu)
        # Ec. (2.6): Cálculo de bits necesarios por empresa
        self.max_shares = [int(budget / p) for p in stock_prices]
        self.bits_per_asset = [math.floor(math.log2(n)) + 1 if n > 0 else 1 for n in self.max_shares]
        
        # Mapeo global de qubits
        self.qubit_map = self._create_qubit_map()
        self.total_qubits = sum(self.bits_per_asset)

    def _create_qubit_map(self):
        mapping = []
        global_idx = 0
        for i in range(self.num_assets):
            asset_qubits = []
            for k in range(self.bits_per_asset[i]):
                asset_qubits.append({
                    'global_index': global_idx,
                    'power': 2**k,
                    'asset_index': i
                })
                global_idx += 1
            mapping.append(asset_qubits)
        return mapping

    def get_ising_coefficients(self):
        """Devuelve h (lineal) y J (cuadrático) para el Hamiltoniano de Ising."""
        Q = self._build_qubo_matrix()
        bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
        h, J, offset = bqm.to_ising()
        return h, J, offset

    def _build_qubo_matrix(self):
        """Construye la matriz Q completa."""
        Q = {}
        for i in range(self.num_assets):
            for bit_i in self.qubit_map[i]:
                idx_i = bit_i['global_index']
                pow_i = bit_i['power']
                price_i = self.prices[i]
                
                for j in range(self.num_assets):
                    for bit_j in self.qubit_map[j]:
                        idx_j = bit_j['global_index']
                        pow_j = bit_j['power']
                        price_j = self.prices[j]
                        
                        # Términos base
                        coeff = 0.0
                        # Riesgo: alpha * sigma_ij * n_i * n_j
                        coeff += self.alpha * self.sigma[i][j] * pow_i * pow_j
                        # Presupuesto: lambda * p_i * p_j * 2^k * 2^l
                        coeff += self.lam * price_i * price_j * pow_i * pow_j
                        
                        if idx_i == idx_j:
                            # Diagonal: Restar beneficio (lineal) y término lineal de presupuesto
                            coeff -= self.mu[i] * pow_i
                            coeff -= 2 * self.lam * self.budget * price_i * pow_i
                            Q[(idx_i, idx_i)] = Q.get((idx_i, idx_i), 0.0) + coeff
                        elif idx_i < idx_j:
                            # Fuera diagonal: Multiplicar por 2 (simetría)
                            Q[(idx_i, idx_j)] = Q.get((idx_i, idx_j), 0.0) + 2 * coeff
        return Q
    def decode_sample(self, sample):
        """Convierte una muestra de espines (+1/-1) a cantidad de acciones (enteros)."""
        shares = [0] * self.num_assets
        for i in range(self.num_assets):
            current_shares = 0
            for bit_info in self.qubit_map[i]:
                # Convertir spin a bit: +1 -> 1, -1 -> 0
                idx = bit_info['global_index']
                # Manejar si el sample viene como dict o lista
                if isinstance(sample, dict):
                    spin_val = sample[idx]
                else:
                    # Si es lista/array, asumimos orden por índice
                    spin_val = sample[idx]
                
                bit_val = 1 if spin_val > 0 else 0
                current_shares += bit_val * bit_info['power']
            shares[i] = current_shares
        return shares