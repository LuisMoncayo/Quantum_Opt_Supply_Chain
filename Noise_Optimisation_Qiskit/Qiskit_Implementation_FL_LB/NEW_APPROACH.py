#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Noise-free QAOA for unbalanced-penalisation Knapsack (Qiskit 1.x)
Author: Luis Moncayo
"""
import numpy as np
from qiskit_optimization import QuadraticProgram
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import QAOAAnsatz
from qiskit_aer.primitives import Estimator, Sampler
from qiskit_algorithms.optimizers import COBYLA

# ------------------------------------------------------------
# 1. Problem data
# ------------------------------------------------------------
values  = [8, 47, 10, 5, 16]
weights = [3, 11, 14, 19, 5]
max_w   = 26
lam1, lam2 = 0.5,0.01
n = len(values)

# ------------------------------------------------------------
# 2. Build QUBO (unbalanced penalisation)
# ------------------------------------------------------------
p = QuadraticProgram()
for i in range(n):
    p.binary_var(name=f"x{i}")

# base objective
p.minimize(linear={f"x{i}": -values[i] for i in range(n)})

# initialise quadratic
for i in range(n):
    for j in range(n):
        p.objective.quadratic[(f"x{i}", f"x{j}")] = 0.0

# λ₁ and λ₂ terms
for i in range(n):
    p.objective.linear[f"x{i}"] += lam1 * weights[i]
p.objective.constant += -lam1 * max_w

for i in range(n):
    for j in range(n):
        p.objective.quadratic[(f"x{i}", f"x{j}")] += lam2 * weights[i] * weights[j]
for i in range(n):
    p.objective.linear[f"x{i}"] += -2 * lam2 * max_w * weights[i]
p.objective.constant += lam2 * max_w**2

# ------------------------------------------------------------
# 3. Convert to Ising
# ------------------------------------------------------------
H, offset = p.to_ising()
if len(H) == 0:
    H = SparsePauliOp.from_list([("I"*n, 0.0)])
print(f"Hamiltonian built with {len(H)} Pauli terms")

# ------------------------------------------------------------
# 4. Build ansatz and Estimator (noise-free)
# ------------------------------------------------------------
reps = 20
ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
estimator = Estimator()

def expected_energy(params):
    """Compute ⟨ψ(θ)|H|ψ(θ)⟩ using Estimator primitive."""
    job = estimator.run([ansatz], [H], [params])
    return job.result().values[0]

# ------------------------------------------------------------
# 5. Optimise γ, β with COBYLA
# ------------------------------------------------------------
opt = COBYLA(maxiter=3000)
theta0 = np.random.uniform(0, 2*np.pi, ansatz.num_parameters)
res = opt.minimize(fun=expected_energy, x0=theta0)
theta_star = res.x
print("\nBest expected energy:", res.fun)
print("Optimal parameters (θ*):", np.round(theta_star, 4))

# ------------------------------------------------------------
# 6. Sample bitstrings (no noise)
# ------------------------------------------------------------
qc = ansatz.assign_parameters(theta_star)
qc.measure_all()
sampler = Sampler()
job = sampler.run([qc], shots=5000)
counts = job.result().quasi_dists[0].binary_probabilities()

# ------------------------------------------------------------
# 7. Display results
# ------------------------------------------------------------
print("\nTop bitstrings:")
for k, v in sorted(counts.items(), key=lambda x: -x[1]):
    print(f"{k}: {v:.4f}")

best = max(counts, key=counts.get)
x_vec = np.array(list(map(int, best[::-1])))
value = sum(v for v, xi in zip(values, x_vec) if xi)
weight = sum(w for w, xi in zip(weights, x_vec) if xi)
print("\nBest bitstring:", best)
print("Decoded x:", x_vec.tolist())
print(f"Total value  = {value}")
print(f"Total weight = {weight}")
