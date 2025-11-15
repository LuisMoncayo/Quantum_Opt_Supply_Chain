#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Facility Location / Line Balancing QAOA (Ideal Simulator, Multi-start)
Exports all counts and feasibility summaries
Author: Luis Moncayo
"""

import pandas as pd
import numpy as np
from scipy.optimize import minimize, linprog
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import QAOAAnsatz
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Estimator
from datetime import datetime
import math, sys

pd.set_option("display.max_columns", None)
np.set_printoptions(linewidth=500, suppress=True)

# ============================================================
# 1. Load instance
# ============================================================
path_to_instance = "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Data_Files/Data_Line_Balancing/data_excel.xlsx"
df = pd.read_excel(path_to_instance, sheet_name="LB_5", header=None)

Cycle_Time = int(df.iloc[0, 1])
nu_Operations = int(df.iloc[1, 1])
nu_Precedences = int(df.iloc[2, 1])
nu_Cells = math.ceil(df.iloc[4, 1:nu_Operations + 1].sum() / Cycle_Time)
nu_variables = nu_Operations * nu_Cells
nu_constraints = nu_Cells + nu_Operations + nu_Precedences

A = df.iloc[6:6 + nu_constraints, 0:nu_variables].to_numpy()
c = np.zeros(nu_variables)
sense = df.iloc[6:6 + nu_constraints, nu_variables].astype(str).values
b = df.iloc[6:6 + nu_constraints, nu_variables + 1].astype(int).values

print("\n--- Dimension verification ---")
print(f"A: {A.shape}, c: {c.shape}, sense: {sense.shape}, b: {b.shape}")

# ============================================================
# 2. Build QUBO (unbalanced penalisation)
# ============================================================
def build_qubo_unbalanced(c, A, sense, b, lam1, lam2):
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    n = len(c)
    Q = np.zeros((n, n))
    const = 0.0

    for i in range(n):
        Q[i, i] += c[i]

    for i in range(A.shape[0]):
        a = A[i]
        bi = float(b[i])
        lin_sign = -1.0
        Q[np.arange(n), np.arange(n)] += lin_sign * (-lam1 * a)
        const += lam1 * lin_sign * bi
        const += lam2 * bi**2
        Q[np.arange(n), np.arange(n)] += lam2 * (a**2) - 2 * lam2 * bi * a
        for j in range(n):
            for k in range(j + 1, n):
                Q[j, k] += 2 * lam2 * a[j] * a[k]

    Q = np.triu(Q)
    Q = (Q + Q.T) / 2     # enforce symmetry
    return Q, const


def qubo_to_pauli(Q, const):
    n = Q.shape[0]
    I = const
    Zi = np.zeros(n)
    ZZ = {}

    for i in range(n):
        I += Q[i, i] / 2
        Zi[i] -= Q[i, i] / 2

    for i in range(n):
        for j in range(i + 1, n):
            val = Q[i, j]
            if val != 0.0:
                I += val / 4
                Zi[i] -= val / 4
                Zi[j] -= val / 4
                ZZ[(i, j)] = ZZ.get((i, j), 0.0) + val / 4

    paulis = [("I" * n, I)]
    for i, coeff in enumerate(Zi):
        if coeff != 0.0:
            s = ["I"] * n
            s[i] = "Z"
            paulis.append(("".join(s), coeff))
    for (i, j), coeff in ZZ.items():
        if coeff != 0.0:
            s = ["I"] * n
            s[i] = "Z"; s[j] = "Z"
            paulis.append(("".join(s), coeff))
    return SparsePauliOp.from_list(paulis)


lambda1, lambda2 = 2.0, 0.1
Q, const = build_qubo_unbalanced(c, A, sense, b, lambda1, lambda2)
H = qubo_to_pauli(Q, const)

# ============================================================
# 3. Ideal (noiseless) simulator
# ============================================================
print("\n✅ Running ideal (noise-free) QAOA simulation")
ideal_backend = AerSimulator(method="statevector")

# ============================================================
# 4. Feasibility check (classical)
# ============================================================
print("\n--- Feasibility check (classical) ---")
A_ub, b_ub, A_eq, b_eq = [], [], [], []
for i, s in enumerate(sense):
    if s.strip() == "<=":
        A_ub.append(A[i]); b_ub.append(b[i])
    elif s.strip() == ">=":
        A_ub.append(-A[i]); b_ub.append(-b[i])
    elif s.strip() == "=":
        A_eq.append(A[i]); b_eq.append(b[i])
A_ub = np.array(A_ub) if A_ub else None
b_ub = np.array(b_ub) if b_ub else None
A_eq = np.array(A_eq) if A_eq else None
b_eq = np.array(b_eq) if b_eq else None
res = linprog(c=np.zeros(A.shape[1]), A_ub=A_ub, b_ub=b_ub,
              A_eq=A_eq, b_eq=b_eq, bounds=(0, 1), method="highs")
if res.success:
    print("✅ Classical check: feasible region exists.")
else:
    print("❌ No feasible region. Verify A, b, or sense.")
    sys.exit()

# ============================================================
# 5. QAOA execution (multi-start)
# ============================================================
N_RUNS = 5
reps = 10
n_shots = 500

print(f"\nConfiguration:")
print(f"  N_RUNS  = {N_RUNS}")
print(f"  reps    = {reps}")
print(f"  n_shots = {n_shots}")
print(f"  lambda1 = {lambda1}")
print(f"  lambda2 = {lambda2}")

ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
estimator = Estimator()

def expected_energy(params):
    job = estimator.run([ansatz], [H], [params])
    return job.result().values[0]

def check_feasibility(counts, A, b, sense, reverse=False):
    rows = []
    for key, freq in counts.items():
        bits = key[::-1] if reverse else key
        x_vec = np.array(list(map(int, bits)))
        lhs = A @ x_vec
        feasible = True
        for i, s in enumerate(sense):
            s = s.strip()
            if s == "<=":
                ok = lhs[i] <= b[i]
            elif s == "=":
                ok = np.isclose(lhs[i], b[i])
            elif s == ">=":
                ok = lhs[i] >= b[i]
            else:
                raise ValueError(f"Invalid sense: {s}")
            feasible &= ok
        rows.append({
            "string": key,
            "freq": freq,
            "sum_string": np.sum(x_vec),
            "feasible": feasible,
            "reverse": reverse
        })
    return pd.DataFrame(rows)

# ============================================================
# 6. Runs
# ============================================================
for run_id in range(1, N_RUNS + 1):
    print(f"\n================ RUN {run_id:02d} / {N_RUNS} ================")
    start_time = datetime.now()

    p = ansatz.reps
    best_fun = np.inf
    best_theta = None
    print("Running multi-start optimisation with 10 attempts...")

    for attempt in range(10):
        gam = np.random.uniform(0, 2 * np.pi, size=p)
        bet = np.random.uniform(0, np.pi / 2, size=p)
        theta0 = np.ravel(np.column_stack([gam, bet]))
        res = minimize(expected_energy, theta0, method="COBYLA",
                       options={"maxiter": 200, "rhobeg": 0.1, "disp": False})
        print(f"  Attempt {attempt+1:02d}: energy = {res.fun:.4f}")
        if res.fun < best_fun:
            best_fun = res.fun
            best_theta = res.x

    theta_star = best_theta
    print(f"Best expected energy after multi-start: {best_fun:.4f}")
    print("Optimal parameters (theta_star):", np.round(theta_star, 4))

    # ============================================================
    # 7. Sampling and feasibility (with all counts)
    # ============================================================
    qc_opt = ansatz.assign_parameters(theta_star)
    qc_opt.measure_all()
    qc_t = transpile(qc_opt, backend=ideal_backend, optimization_level=1)
    job = ideal_backend.run(qc_t, shots=n_shots)
    counts = job.result().get_counts()

    # ---- Save raw counts ----
    counts_df = pd.DataFrame(list(counts.items()), columns=["string", "freq"])
    counts_df["string"] = "'" + counts_df["string"].astype(str)
    counts_output_csv = f"counts_run{run_id:02d}_ideal.csv"
    counts_df.to_csv(counts_output_csv, index=False)
    print(f"Raw counts saved: {counts_output_csv} ({len(counts_df)} bitstrings)")

    # ---- Feasibility analysis ----
    df_rev = check_feasibility(counts, A, b, sense, reverse=True)
    df_nrev = check_feasibility(counts, A, b, sense, reverse=False)
    df_all = pd.concat([df_rev, df_nrev], ignore_index=True)

    n_feas_rev = df_rev["feasible"].sum()
    n_feas_nrev = df_nrev["feasible"].sum()
    sum_freq_rev = df_rev.loc[df_rev["feasible"], "freq"].sum()
    sum_freq_nrev = df_nrev.loc[df_nrev["feasible"], "freq"].sum()

    print(f"Feasible (reverse): {n_feas_rev} | Feasible (non-rev): {n_feas_nrev}")
    print(f"Sum(freq feasible) reverse={sum_freq_rev}, non-rev={sum_freq_nrev}")

    # ============================================================
    # 8. Export feasibility results
    # ============================================================
    output_csv = f"feasibility_summary_run{run_id:02d}_ideal.csv"
    df_all["string"] = "'" + df_all["string"].astype(str)
    df_all.to_csv(output_csv, index=False)
    print(f"Feasibility results saved: {output_csv}")

    elapsed = (datetime.now() - start_time).total_seconds()
    print(f"Elapsed: {elapsed:.2f}s")
