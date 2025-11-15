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
from qiskit.visualization import circuit_drawer
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit_aer.primitives import Estimator
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit import transpile

from data.load_data import LoadData
from qubo_matrix.qubo import QuboTransformation

import time
from datetime import datetime

pd.set_option('display.max_columns', None)
np.set_printoptions(linewidth=500, suppress=True)



# ============================================================
# FUNCTION: FULLY DECOMPOSE QAOA INTO GATE-BY-GATE CIRCUIT
# ============================================================
def expand_qaoa_circuit(ansatz, params):
    """Return a fully decomposed QAOA circuit (no opaque QAOA blocks)."""
    # Decompose internal structure of QAOAAnsatz
    qc = ansatz.decompose()

    # Apply parameters
    qc = qc.assign_parameters(params)

    # Deep unrolling: break U-gates → RZ/RX, break ZZ → CX-RZ-CX
    for _ in range(5):
        qc = qc.decompose()

    return qc



# ============================================================
# LOAD GEOJSON + DATA PREPARATION
# ============================================================
geojson_data = gpd.read_file("/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Version_LB_FL/Data_Files/uk-postcode-polygons-master/geojson/EX.geojson")
codes = ['EX5','EX9','EX10','EX11','EX12','EX13']

gdf = geojson_data[geojson_data['name'].apply(lambda x: x in codes)]

problem_data = LoadData(gdf)
constraints, variables, set_constraints = problem_data.compute_set_constriants()
obj_func = np.concatenate((np.ones(gdf.shape[0], dtype=int), np.zeros((variables-gdf.shape[0]), dtype=int)))
set_obj_func = obj_func[np.newaxis, :]

n_variables = len(codes)
print(f"\nThe problem has {constraints} constraints and {n_variables} variables\n")



# ============================================================
# 1. Build QUBO matrix (unbalanced)
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
    Q = Q + Q.T - np.diag(np.diag(Q))
    return Q, const



# ============================================================
# 2. Convert QUBO → Pauli operator
# ============================================================
def qubo_to_pauli(Q, const):
    n = Q.shape[0]
    I = const
    Zi = np.zeros(n)
    ZZ = {}

    for i in range(n):
        I += Q[i, i]/2
        Zi[i] -= Q[i, i]/2

    for i in range(n):
        for j in range(i + 1, n):
            val = Q[i, j]
            if val != 0:
                I += val/4
                Zi[i] -= val/4
                Zi[j] -= val/4
                ZZ[(i,j)] = ZZ.get((i,j), 0) + val/4

    paulis = [("I"*n, I)]

    for i, coeff in enumerate(Zi):
        if coeff != 0:
            s = ["I"]*n
            s[i] = "Z"
            paulis.append(("".join(s), coeff))

    for (i,j), coeff in ZZ.items():
        if coeff != 0:
            s = ["I"]*n
            s[i] = "Z"; s[j] = "Z"
            paulis.append(("".join(s), coeff))

    return SparsePauliOp.from_list(paulis)



# ============================================================
# QUBO DATA
# ============================================================
c = set_obj_func[0, :n_variables]
A = set_constraints[:, :n_variables]
sense = [">="] * n_variables
b = [1] * n_variables

lambda1, lambda2 = 0.5, 0.01

Q, const = build_qubo_unbalanced(c, A, sense, b, lambda1, lambda2)
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
# 3. IBM BACKEND + NOISE MODEL
# ============================================================
try:
    service = QiskitRuntimeService()
    backends = service.backends(simulator=False)
    backend_real = next((b for b in backends if "brisbane" in b.name), backends[0])
    print(f"\nUsing backend: {backend_real.name}")

    noise_model = NoiseModel.from_backend(backend_real)
    noisy_backend = AerSimulator(
        noise_model=noise_model,
        basis_gates=noise_model.basis_gates,
        method="density_matrix"
    )
except:
    print("\n⚠ No IBM backend, using ideal simulator.")
    noisy_backend = AerSimulator(method="statevector")



# ============================================================
# 4. QAOA CONFIG
# ============================================================
N_RUNS = 1
reps = 1
n_shots = 1

ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
ansatz.draw("mpl", fold=120).savefig("QAOA_initial_circuit.png", dpi=300)



# ============================================================
# 5. Optimisation + FULL CIRCUIT DECOMPOSITION
# ============================================================
estimator = Estimator()

def expected_energy(params):
    job = estimator.run([ansatz], [H], [params])
    return job.result().values[0]


for run_id in range(1, N_RUNS + 1):
    print(f"\n========== RUN {run_id} ==========")

    p = ansatz.reps
    gam = np.random.uniform(0, 2*np.pi, size=p)
    bet = np.random.uniform(0, np.pi/2, size=p)
    theta0 = np.ravel(np.column_stack([gam, bet]))

    res = minimize(expected_energy, theta0, method="COBYLA",
                   options={"maxiter": 8})
    theta_star = res.x
    print("Optimal parameters:", theta_star)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # FULL QAOA EXPANSION
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    qc_param = expand_qaoa_circuit(ansatz, theta_star)
    #qc_param.measure_all()
    qc_param.remove_final_measurements()

    # Save expanded diagram (with U gate boxes)
    qc_param.draw("mpl", fold=180).savefig(
        "QAOA_expanded_layer.png",
        dpi=300,
        bbox_inches="tight"
    )

    # === PRINT FULLY DECOMPOSED CIRCUIT ===
    print("\n\n===== FULLY DECOMPOSED CIRCUIT =====\n")
    print(qc_param)
    print("\n====================================\n")

    print("\nSaved fully decomposed circuit as QAOA_expanded_layer.png")
