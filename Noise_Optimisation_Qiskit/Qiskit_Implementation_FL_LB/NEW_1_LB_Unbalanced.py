#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Line Balancing QAOA — Unbalanced Penalisation Sweep
Author: Luis Moncayo
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.converters import QuadraticProgramToQubo
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import QAOAAnsatz
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Estimator, Sampler
from qiskit_algorithms.optimizers import COBYLA

# ============================================================
# 1. Load instance
# ============================================================
path = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Data_Files/Data_Line_Balancing/data_excel.xlsx"
df = pd.read_excel(path, sheet_name="LB_7", header=None)

nu_variables, nu_constraints = 21, 18
A = df.iloc[6:6+nu_constraints, 0:nu_variables].to_numpy(dtype=float)
sense = df.iloc[6:6+nu_constraints, nu_variables].astype(str).values
b = df.iloc[6:6+nu_constraints, nu_variables+1].astype(float).values

print(f"--- Dimensions ---\nA: {A.shape}, sense: {sense.shape}, b: {b.shape}")

# ============================================================
# 2. Row normalisation (important for stable penalties)
# ============================================================
row_norms = np.maximum(np.sum(np.abs(A), axis=1), 1e-9)
A_norm = A / row_norms[:, None]
b_norm = b / row_norms

# ============================================================
# 3. Penalty sweep grid
# ============================================================
lambda1_list = [1.0]
lambda2_list = [8,10,12,14,16,18]
reps = 3

results_summary = []

# ============================================================
# 4. Common simulator & primitives
# ============================================================
backend = AerSimulator(method="matrix_product_state")
estimator = Estimator(run_options={"backend": backend})
sampler = Sampler(run_options={"backend": backend})

# ============================================================
# 5. Sweep loop
# ============================================================
for λ1 in lambda1_list:
    for λ2 in lambda2_list:

        print(f"\n=== λ₁={λ1}, λ₂={λ2} ===")

        # ---------- build QUBO manually ----------
        n = nu_variables
        Q = np.zeros((n, n))
        c = np.zeros(n)
        offset = 0.0

        for i in range(nu_constraints):
            a = A_norm[i]
            bi = b_norm[i]

            if sense[i].strip() == "<=":
                Q += λ2 * np.outer(a, a)
                c += -λ1 * a - 2 * λ2 * bi * a
                offset += λ1 * bi + λ2 * bi**2
            elif sense[i].strip() == "=":
                Q += λ2 * np.outer(a, a)
                c += -2 * λ2 * bi * a
                offset += λ2 * bi**2
            else:
                raise ValueError(f"Unknown sense {sense[i]}")

        # symmetrise and scale off-diagonals (QAOA expects upper-triangular)
        Q = np.triu(Q) + np.triu(Q, 1).T
        Q[np.triu_indices_from(Q, 1)] *= 2

        # ---------- build Hamiltonian ----------
        paulis, coeffs = [], []
        for i in range(n):
            if abs(c[i]) > 1e-12:
                paulis.append("Z" * (n - i - 1) + "Z" + "I" * i)
                coeffs.append(0.5 * c[i])
            for j in range(i + 1, n):
                if abs(Q[i, j]) > 1e-12:
                    z_string = ["I"] * n
                    z_string[i] = z_string[j] = "Z"
                    paulis.append("".join(z_string[::-1]))
                    coeffs.append(0.25 * Q[i, j])

        H = SparsePauliOp.from_list(list(zip(paulis, coeffs)))
        print(f"✅ QUBO via unbalanced penalisation: {len(H)} Pauli terms")

        # ---------- QAOA ansatz ----------
        ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
        print("Circuit depth:", ansatz.decompose().depth())

        # ---------- short multi-start COBYLA ----------
        best_fun, best_theta = np.inf, None
        for attempt in range(5):
            theta0 = np.random.uniform(0, 2*np.pi, ansatz.num_parameters)
            opt = COBYLA(maxiter=300)
            try:
                res = opt.minimize(lambda p: estimator.run([ansatz],[H],[p]).result().values[0], x0=theta0)
                if res.fun < best_fun:
                    best_fun, best_theta = res.fun, res.x
            except Exception as e:
                print("  Estimator failed:", e)
        print(f"✅ Best energy: {best_fun:.6f}")

        # ---------- sample ----------
        qc = ansatz.assign_parameters(best_theta)
        qc.measure_all()
        result = sampler.run([qc], shots=4096).result()
        counts = result.quasi_dists[0].binary_probabilities()

        # ---------- feasibility analysis ----------
        rows = []
        for bitstring, prob in counts.items():
            x_vec = np.array(list(map(int, bitstring[::-1])))
            lhs = A @ x_vec
            feas = True
            viol_sum = 0.0
            for i, s in enumerate(sense):
                r = lhs[i] - b[i]
                if s.strip() == "<=":
                    feas &= (r <= 1e-6)
                    viol_sum += max(0, r)
                elif s.strip() == "=":
                    feas &= abs(r) <= 1e-6
                    viol_sum += abs(r)
            rows.append({"bitstring": "'" + bitstring, "prob": prob,
                         "violation": viol_sum, "feasible": feas})
        df_r = pd.DataFrame(rows)
        feas_mass = df_r.loc[df_r.feasible, "prob"].sum()
        n_feas = df_r.feasible.sum()

        print(f"Feasible mass={feas_mass:.2e}, n_feasible={n_feas}")

        # ---------- save & plot ----------
        tag = f"L1_{λ1}_L2_{λ2}_{datetime.now().strftime('%H%M%S')}"
        df_r.to_csv(f"LB5_QAOA_Unbal_{tag}.csv", index=False)

        plt.figure(figsize=(6,4))
        plt.scatter(df_r.violation, df_r.prob, s=10)
        plt.yscale("log")
        plt.xlabel("Total constraint violation (sum)")
        plt.ylabel("Probability (log scale)")
        plt.title(f"λ₁={λ1}, λ₂={λ2}")
        plt.tight_layout()
        plt.savefig(f"LB5_QAOA_Unbal_{tag}.png", dpi=200)
        plt.close()

        results_summary.append((λ1, λ2, feas_mass, n_feas, best_fun))

# ============================================================
# 6. Summary table
# ============================================================
df_sum = pd.DataFrame(results_summary, columns=["λ1","λ2","feas_mass","n_feas","best_energy"])
df_sum.to_csv("LB5_QAOA_PenaltySweep_Summary.csv", index=False)
print("\n=== Sweep summary ===")
print(df_sum)
