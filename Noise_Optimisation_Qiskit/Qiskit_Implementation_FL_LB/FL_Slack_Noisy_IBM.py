#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Facility Location QAOA with real IBM backend noise model
Multiple QAOA runs (independent initialisations)
Author: Luis Moncayo
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

# ============================================================
# 1. Load geographic data
# ============================================================
geojson_data = gpd.read_file(
    "/Users/luismoncayo/Library/CloudStorage/Dropbox/Python/Qiskit_Implementation_FL_LB/Data_Files/uk-postcode-polygons-master/geojson/EX.geojson"
)
#codes = ['EX5','EX9','EX10','EX11','EX12','EX13']
#odes = ['EX1','EX2','EX3','EX4','EX5','EX6','EX7']
codes = ['EX9','EX10','EX11','EX12','EX13','EX14','EX15','EX24']
gdf = geojson_data[geojson_data['name'].apply(lambda x: x in codes)]

problem_data = LoadData(gdf)
constraints, variables, set_constraints = problem_data.compute_set_constriants()

obj_func = np.concatenate(
    (np.ones(gdf.shape[0], dtype=int),
     np.zeros((variables - gdf.shape[0]), dtype=int))
)
set_obj_func = obj_func[np.newaxis, :]
matrix_computation = QuboTransformation(set_obj_func, set_constraints)
lagrange = 1
offset, QT = matrix_computation.compute_Q(lagrange)

# ============================================================
# 2. QUBO → SparsePauliOp
# ============================================================
def qubo_to_pauli(Q, const):
    n = Q.shape[0]
    I = const
    Zi = np.zeros(n)
    ZZ = {}
    for i in range(n):
        I += Q[i, i] / 2.0
        Zi[i] -= Q[i, i] / 2.0
    for i in range(n):
        for j in range(i + 1, n):
            val = Q[i, j]
            if val != 0.0:
                I += val / 4.0
                Zi[i] -= val / 4.0
                Zi[j] -= val / 4.0
                ZZ[(i, j)] = ZZ.get((i, j), 0.0) + val / 4.0
    paulis = [("I" * n, I)]
    for i, c in enumerate(Zi):
        if c:
            s = ["I"] * n
            s[i] = "Z"
            paulis.append(("".join(s), c))
    for (i, j), c in ZZ.items():
        if c:
            s = ["I"] * n
            s[i] = "Z"
            s[j] = "Z"
            paulis.append(("".join(s), c))
    return SparsePauliOp.from_list(paulis)

H = qubo_to_pauli(QT, offset)
print("\nSparsePauliOp built successfully.")

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
        noise_model=noise_model.simplify(),   # removes tiny Kraus terms
        basis_gates=noise_model.basis_gates,
        method="density_matrix",
        max_parallel_threads=8,               # use all cores
        max_parallel_experiments=2,
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
N_RUNS = 30
reps = 3
n_shots = 250#256
ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
estimator = Estimator()

print(f"\nConfiguration:")
print(f"  N_RUNS  = {N_RUNS}")
print(f"  reps    = {reps}")
print(f"  n_shots = {n_shots}")
print(f"  lagrange = {lagrange}\n")

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
    
    # ---- Decode and evaluate ----
    n_vars = gdf.shape[0]
    rows = []
    for s, c in sorted_counts.items():
        arr = np.array(list(s[::-1]), dtype=int)
        x = arr[:n_vars]
        z = arr
        energy_M = int(x.sum())
        energy_Q = float(z @ QT @ z + offset)
        feasible = bool(np.all(set_constraints[:, :n_vars] @ x >= 1))
        rows.append(("'" + s, c, energy_M, energy_Q, "feasible" if feasible else "infeasible"))

    df = pd.DataFrame(rows, columns=["string", "counts", "energy_M", "energy_Q", "feasible"])
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