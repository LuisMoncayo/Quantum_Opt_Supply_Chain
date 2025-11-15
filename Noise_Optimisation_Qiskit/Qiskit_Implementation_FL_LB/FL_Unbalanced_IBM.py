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
from qiskit.visualization import circuit_drawer
from data.load_data import LoadData
from qubo_matrix.qubo import QuboTransformation
import time
from datetime import datetime

pd.set_option('display.max_columns', None)
np.set_printoptions(linewidth=500, suppress=True)


# import geopandas as gpd
# import pandas as pd
# pd.set_option('display.max_columns', None)  # Show all columns
# import numpy as np
# np.set_printoptions(linewidth=500, suppress=True)
# from data.load_data import LoadData
# from qubo_matrix.qubo import QuboTransformation
# #from optimisation.compute_solution import ComputeSolution
# #from problems.problems import ProblemsToSolve
# import os

# from scipy.optimize import minimize
# from qiskit.quantum_info import SparsePauliOp
# from qiskit_algorithms.optimizers import SPSA
# from qiskit.circuit.library import QAOAAnsatz
# from qiskit_aer.primitives import Estimator
# from qiskit_aer import AerSimulator
# from qiskit import transpile


geojson_data = gpd.read_file("/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Version_LB_FL/Data_Files/uk-postcode-polygons-master/geojson/EX.geojson")
#geojson_data = gpd.read_file("/home/luismoncayo/Dropbox/Python/PennyLane_QUBO/Data_Files/uk-postcode-polygons-master/geojson/EX.geojson")
codes = ['EX5','EX9','EX10','EX11','EX12','EX13']
#codes = ['EX1','EX2','EX3','EX4','EX5','EX6','EX7']
#codes = ['EX9','EX10','EX11','EX12','EX13','EX14','EX15','EX24']

gdf = geojson_data[geojson_data['name'].apply(lambda x: x in codes)]

problem_data = LoadData(gdf)
#problem_data.print_map()
constraints, variables, set_constraints = problem_data.compute_set_constriants()
obj_func = np.concatenate((np.ones(gdf.shape[0], dtype=int), np.zeros((variables-gdf.shape[0]), dtype=int)))
set_obj_func = obj_func[np.newaxis, :]

n_variables = len(codes)

print(f"\nThe problem has {constraints} constraints and {n_variables} variables\n")
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
c = set_obj_func[0, :n_variables]
A = square_matrix = set_constraints[:, :n_variables]
sense = [">="] * n_variables
b = [1] * n_variables

lambda1, lambda2 = 0.5, 0.01 

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
        method="density_matrix"
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
n_shots =1#256
ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
# --- Visualise and save initial QAOA circuit ---
ansatz.draw("mpl", fold=100).savefig("QAOA_initial_circuit.png", dpi=300, bbox_inches="tight")
print("Initial QAOA circuit saved as QAOA_initial_circuit.png")


estimator = Estimator()

print(f"\nConfiguration:")
print(f"  N_RUNS  = {N_RUNS}")
print(f"  reps    = {reps}")
print(f"  n_shots = {n_shots}")
print(f"  lambda2 = {lambda2}")
print(f"  lambda1 = {lambda1}\n")

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
    # --- Visualise and save optimised QAOA circuit ---
    # qc_opt.decompose().draw("mpl", fold=100).savefig("QAOA_optimised_circuit.png", dpi=300, bbox_inches="tight")
    # fig = circuit_drawer(qc_opt, output='mpl')
    # fig.savefig("qaoa_qiskit_expanded.png", dpi=300, bbox_inches='tight')
    # print(qc_opt)
    # print(qc_opt.decompose().draw('text'))
    # print("Optimised QAOA circuit saved as qaoa_qiskit_expanded.png")

    qc_t = transpile(qc_opt, noisy_backend, optimization_level=1)
    print("Transpiled circuit depth:", qc_t.depth())
    job = noisy_backend.run(qc_t, shots=n_shots)
    counts = job.result().get_counts()
    sorted_counts = dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))
    end_time = datetime.now()
    elapsed = end_time - start_time
    print(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed time: {elapsed.total_seconds():.2f} seconds")
    
    # ---- Decode and evaluate ----
    n_vars = gdf.shape[0]
    rows = []
    for s, c in sorted_counts.items():
        arr = np.array(list(s[::-1]), dtype=int)
        x = arr[:n_vars]
        z = arr
        energy_M = int(x.sum())
        #energy_Q = float(z @ QT @ z + offset)
        feasible = bool(np.all(set_constraints[:, :n_vars] @ x >= 1))
        rows.append(("'" + s, c, energy_M, "feasible" if feasible else "infeasible"))

    df = pd.DataFrame(rows, columns=["string", "counts", "energy_M", "feasible"])
    df_sorted = df.sort_values(by=["feasible", "energy_M"], ascending=[False, True]).reset_index(drop=True)

    # ---- Save file with run index ----
    file_name = f"facility_location_results_noisy_run{run_id:02d}.csv"
    df_sorted.to_csv(file_name, index=False)
    print(f"File saved as {file_name}")

    # ---- Summary ----
    count_opt = df[(df["energy_M"] == 2) & (df["feasible"] == "feasible")]
    sum_count_opt = count_opt["counts"].sum()
    feasibles = df[df["feasible"] == "feasible"]
    sum_count_feas = feasibles["counts"].sum()
    print(f"Feasible solutions with energy_M = 2 (total shots): {sum_count_opt}")
    print(f"Total feasible measurements (total shots): {sum_count_feas}")



























