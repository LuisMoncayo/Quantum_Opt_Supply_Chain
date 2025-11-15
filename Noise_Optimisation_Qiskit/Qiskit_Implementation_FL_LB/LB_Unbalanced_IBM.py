#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 12:52:10 2025

@author: luismoncayo
"""
import geopandas as gpd
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import QAOAAnsatz
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit_aer.primitives import Estimator
from qiskit_ibm_runtime import QiskitRuntimeService
from data.load_data import LoadData
from qubo_matrix.qubo import QuboTransformation
import time
from datetime import datetime

pd.set_option('display.max_columns', None)
np.set_printoptions(linewidth=500, suppress=True)


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
c = np.zeros(nu_variables, dtype=int)
A = df.iloc[6:6 + nu_constraints, 0:nu_variables].to_numpy(int)
sense = df.iloc[6:6 + nu_constraints, nu_variables].astype(str).values
b = df.iloc[6:6 + nu_constraints, nu_variables + 1].to_numpy(int)
print(f"--- Dimensions ---\nA: {A.shape}, sense: {sense.shape}, b: {b.shape}")

# ============================================================
# 1. Build QUBO matrix (unbalanced penalisation)
# ============================================================
def build_qubo_unbalanced(c, A, sense, b, lam1, lam2):
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    n = len(c)
    Q = np.zeros((n, n))
    const = 0.0

    # Objective diagonal
    for i in range(n):
        Q[i, i] += c[i]

    for i in range(A.shape[0]):
        a = A[i]
        bi = float(b[i])
        s = sense[i].strip()

        if s == "<=":
            lin_sign = -1.0
        elif s == "=":
            lin_sign = -1.0
        elif s == ">=":
            lin_sign = -1.0
        else:
            raise ValueError("sense must be one of '<=', '=', '>='")

        # Linear part: -λ₁(a·x - b)
        Q[np.arange(n), np.arange(n)] += lin_sign * (-lam1 * a)
        const += lam1 * lin_sign * bi

        # Quadratic part: +λ₂(a·x - b)²
        const += lam2 * bi**2
        Q[np.arange(n), np.arange(n)] += lam2 * (a**2) - 2 * lam2 * bi * a
        for j in range(n):
            for k in range(j + 1, n):
                Q[j, k] += 2 * lam2 * a[j] * a[k]

    # make symmetric first
    Q = np.triu(Q)
    Q = Q + Q.T - np.diag(np.diag(Q))
    return Q, const

# ============================================================
# 2. Convert QUBO → SparsePauliOp
# ============================================================
def qubo_to_pauli(Q, const):
    n = Q.shape[0]
    I = const
    Zi = np.zeros(n)
    ZZ = {}

    # Linear terms
    for i in range(n):
        I += Q[i, i] / 2
        Zi[i] -= Q[i, i] / 2

    # Quadratic terms
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


# ============================================================
# 3. Model data
# ============================================================

lambda1, lambda2 = 1,12 

Q, const = build_qubo_unbalanced(c, A, sense, b, lambda1, lambda2)

# --- enforce upper-triangular 2Qij convention ---
Q = np.triu(Q)
for i in range(Q.shape[0]):
    for j in range(i + 1, Q.shape[1]):
        Q[i, j] *= 2.0

print("Upper-triangular QUBO matrix:\n", Q)
print("Constant term:", const)

H = qubo_to_pauli(Q, const)
print("\nSparsePauliOp built successfully.")
print(H)

# ============================================================
# 3. Connect to IBM Quantum and build noise model
# ============================================================
try:
    service = QiskitRuntimeService()  # requires saved IBM token
    backends = service.backends(simulator=False)
    backend_real = next((b for b in backends if "brisbane" in b.name), backends[0])
    print(f"\n✅ Using backend: {backend_real.name} ({backend_real.configuration().num_qubits} qubits)")
    noise_model = NoiseModel.from_backend(backend_real)
    noisy_backend = AerSimulator(
        noise_model=noise_model,
        basis_gates=noise_model.basis_gates,
        method="matrix_product_state"
        )

    print("Noise model successfully loaded from IBM calibration data.")
except Exception as e:
    print("\n⚠️ Could not connect to IBM Quantum or load noise model. Falling back to ideal AerSimulator.")
    print("Error:", e)
    noise_model = None
    noisy_backend = AerSimulator(method="statevector")
    backend_real = None

# ============================================================
# 4–7. Repeat 30 runs with real-device noise
# ============================================================
N_RUNS = 1
reps = 1
n_shots = 1#256
ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
estimator = Estimator()

print(f"\nConfiguration:")
print(f"  N_RUNS  = {N_RUNS}")
print(f"  reps    = {reps}")
print(f"  n_shots = {n_shots}")
print(f"  lambda1 = {lambda1}")
print(f"  lambda2 = {lambda2}\n")

def expected_energy(params):
    job = estimator.run([ansatz], [H], [params])
    return job.result().values[0]

for run_id in range(1, N_RUNS + 1):
    print(f"\n================ RUN {run_id:02d} / {N_RUNS} ================")
    start_time = datetime.now()
    print(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    # ---- QAOA optimisation (ideal) ----
    p = ansatz.reps
    gam = np.random.uniform(0, 2 * np.pi, size=p)
    bet = np.random.uniform(0, np.pi / 2, size=p)
    theta0 = np.ravel(np.column_stack([gam, bet]))
    res = minimize(expected_energy, theta0, method="COBYLA",
                   options={"maxiter": 8, "rhobeg": 0.5, "disp": False})
    theta_star = res.x
    print(f"Best expected energy (ideal): {res.fun:.4f}")
    print("Optimal parameters (theta_star):", np.round(theta_star, 4))

    # ---- Sampling with real-device noise ----
    qc_opt = ansatz.assign_parameters(theta_star)
    qc_opt.measure_all()
    qc_t = transpile(qc_opt, noisy_backend, optimization_level=1)
    print("Transpiled circuit depth:", qc_t.depth())
    job = noisy_backend.run(qc_t, shots=n_shots)
    counts = job.result().get_counts()
    sorted_counts = dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))
    end_time = datetime.now()
    elapsed = end_time - start_time
    print(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed time: {elapsed.total_seconds():.2f} seconds")
    
    # ---- Decode and feasibility check ----
    n_vars = A.shape[1]
    rows = []
    for bitstring, freq in sorted_counts.items():
        x = np.array(list(bitstring[::-1]), dtype=int)  # reverse order
        lhs = A @ x
        feasible = np.all(
            (sense == ">=") & (lhs >= b) |
            (sense == "<=") & (lhs <= b) |
            (sense == "=") & (lhs == b)
        )
        rows.append({
            "bitstring": bitstring,
            "count": freq,
            "feasible": bool(feasible),
            "violations": np.sum(~(
                ((sense == ">=") & (lhs >= b)) |
                ((sense == "<=") & (lhs <= b)) |
                ((sense == "=") & (lhs == b))
            ))
        })
    
    df_res = pd.DataFrame(rows)
    # ---- Add column with sum of bits ----
    df_res["sum_bits"] = df_res["bitstring"].apply(lambda s: sum(int(ch) for ch in s))
    n_feas = df_res["feasible"].sum()
    print(f"\nNumber of feasible solutions: {n_feas} / {len(df_res)}")
    
    # ---- Save to CSV ----
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"LB5_QAOA_Run{run_id:02d}_Results_{timestamp}.csv"
    df_res.to_csv(filename, index=False)

    
    
    
    
    
    
    
    
