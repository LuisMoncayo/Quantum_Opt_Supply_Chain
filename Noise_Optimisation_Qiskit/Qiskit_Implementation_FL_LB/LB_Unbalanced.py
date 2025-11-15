#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Line Balancing QAOA — Unbalanced Penalisation (Statistical Analysis)
Performs N_RUNS independent initialisations for fixed λ₁–λ₂.
Author: Luis Moncayo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import QAOAAnsatz
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Estimator, Sampler
from qiskit_algorithms.optimizers import COBYLA
import time  # Add this at the top of your script

# ============================================================
# 1. Parameters
# ============================================================
N_RUNS = 30
λ1 = 1.0
λ2 = 12          # ← choose your preferred λ₂
reps = 20
shots = 5000

# ============================================================
# 2. Load instance
# ============================================================
path = (
    "/Users/luismoncayo/Library/CloudStorage/Dropbox/"
    "Python/Qiskit_Implementation_FL_LB/Data_Files/"
    "Data_Line_Balancing/data_excel.xlsx"
)
df = pd.read_excel(path, sheet_name="LB_5", header=None)

nu_variables, nu_constraints = 15, 12
A = df.iloc[6:6 + nu_constraints, 0:nu_variables].to_numpy(float)
sense = df.iloc[6:6 + nu_constraints, nu_variables].astype(str).values
b = df.iloc[6:6 + nu_constraints, nu_variables + 1].to_numpy(float)
print(f"--- Dimensions ---\nA: {A.shape}, sense: {sense.shape}, b: {b.shape}")

# ============================================================
# 3. Row normalisation
# ============================================================
row_norms = np.maximum(np.sum(np.abs(A), axis=1), 1e-9)
A_norm = A / row_norms[:, None]
b_norm = b / row_norms

# ============================================================
# 4. Build QUBO
# ============================================================
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

Q = np.triu(Q) + np.triu(Q, 1).T
Q[np.triu_indices_from(Q, 1)] *= 2

# ============================================================
# 5. Build Hamiltonian
# ============================================================
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
print(f"✅ QUBO built: {len(H)} Pauli terms")

# ============================================================
# 6. Simulator setup
# ============================================================
backend = AerSimulator(method="matrix_product_state")
estimator = Estimator(run_options={"backend": backend})
sampler = Sampler(run_options={"backend": backend})

# ============================================================
# 7. Statistical loop
# ============================================================
results = []

ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
print("Circuit depth:", ansatz.decompose().depth())

for run in range(1, N_RUNS + 1):
    print(f"\n================ RUN {run:02d} / {N_RUNS} ================")
    
    start_time = time.time()
    start_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    theta0 = np.random.uniform(0, 2*np.pi, ansatz.num_parameters)
    opt = COBYLA(maxiter=300)

    try:
        res = opt.minimize(lambda p: estimator.run([ansatz], [H], [p]).result().values[0], x0=theta0)
        best_fun = res.fun
        best_theta = res.x
    except Exception as e:
        print("  Estimator failed:", e)
        continue

    # ---------- Sampling ----------
    qc = ansatz.assign_parameters(best_theta)
    qc.measure_all()
    result = sampler.run([qc], shots=shots).result()
    counts = result.quasi_dists[0].binary_probabilities()
    
    
    end_time = time.time()
    end_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    elapsed = end_time - start_time


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
        rows.append({"bitstring": "'" + bitstring, "prob": prob, "violation": viol_sum, "feasible": feas})
    df_r = pd.DataFrame(rows)
    
    # ---------- Save per-shot data ----------
    bit_feas_map = df_r.set_index("bitstring")["feasible"].to_dict()

    shot_records = []
    for bitstring, prob in counts.items():
        n_occurrences = int(round(prob * shots))  # number of times observed
        x_vec = np.array(list(map(int, bitstring[::-1])))
        sum_x = int(x_vec.sum())

        shot_records.append({
            "run": run,
            "bitstring": "'" + bitstring,
            "probability": prob,
            "count": n_occurrences,
            "sum_x": sum_x,
            "feasible": bit_feas_map.get("'" + bitstring, False)
        })

    df_shots = pd.DataFrame(shot_records)
    shot_csv = f"LB5_QAOA_Run{run:02d}_Shots_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    df_shots.to_csv(shot_csv, index=False)
    print(f"🧾 Saved per-shot file: {shot_csv}")



    feas_mass = df_r.loc[df_r.feasible, "prob"].sum()
    n_feas = df_r.feasible.sum()
    print(f"✅ Energy={best_fun:.6f} | FeasMass={feas_mass:.3e} | nFeas={n_feas}")
    print(f"🕒 Started: {start_dt} | Ended: {end_dt} | Elapsed: {elapsed:.2f} sec")

    results.append((run, best_fun, feas_mass, n_feas))

    # ============================================================
    # 8. Summary statistics
    # ============================================================
    df_runs = pd.DataFrame(results, columns=["run", "best_energy", "feas_mass", "n_feas"])
    df_runs["energy_norm"] = df_runs["best_energy"] / λ2
    summary = df_runs.describe().T
    

