#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Line Balancing QAOA — Unbalanced Penalisation (full Pauli-string form)
Compatible with Qiskit 2.2.2 + Aer 0.17.2
Author: Luis Moncayo
"""

import pandas as pd
import numpy as np
from datetime import datetime
pd.set_option("display.max_columns", None)
np.set_printoptions(linewidth=500, suppress=True)

from qiskit.circuit.library import QAOAAnsatz
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Estimator, Sampler
from qiskit_algorithms.optimizers import COBYLA

# ============================================================
# 1. Load instance
# ============================================================
path = (
    "/Users/luismoncayo/Library/CloudStorage/Dropbox/"
    "Python/Qiskit_Implementation_FL_LB/Data_Files/"
    "Data_Line_Balancing/data_excel.xlsx"
)
df = pd.read_excel(path, sheet_name="LB_5", header=None)

nu_variables = 15
nu_constraints = 12
A = df.iloc[6:6 + nu_constraints, 0:nu_variables].to_numpy(float)
sense = df.iloc[6:6 + nu_constraints, nu_variables].astype(str).values
b = df.iloc[6:6 + nu_constraints, nu_variables + 1].to_numpy(float)

print(f"--- Dimensions ---\nA: {A.shape}, sense: {sense.shape}, b: {b.shape}")

# ============================================================
# 2. QUBO via Unbalanced Penalisation
# ============================================================
λ1, λ2 = 2,0.5
Q = np.zeros((nu_variables, nu_variables))
c = np.zeros(nu_variables)
offset = 0.0

for i in range(nu_constraints):
    a = A[i, :]
    bi = b[i]
    if sense[i].strip() == ">=":
        a, bi = -a, -bi
    Q += λ2 * np.outer(a, a)
    c += -λ1 * a - 2 * λ2 * bi * a
    offset += λ1 * bi + λ2 * bi**2

Q = np.triu(Q)
# for i in range(Q.shape[0]):
#     for j in range(i + 1, Q.shape[1]):
#         Q[i, j] *= 2.0
print(f"\n✅ QUBO constructed via unbalanced penalisation")
print(f"λ₁={λ1}, λ₂={λ2}, Q shape={Q.shape}")

# ============================================================
# 3. Convert to Ising Hamiltonian — full Pauli strings
# ============================================================
def make_Z_string(i, n):
    """Return 'Z' on qubit i, 'I' elsewhere."""
    s = ["I"] * n
    s[i] = "Z"
    return "".join(s)

def make_ZZ_string(i, j, n):
    s = ["I"] * n
    s[i] = "Z"
    s[j] = "Z"
    return "".join(s)

paulis, coeffs = [], []

for i in range(nu_variables):
    coeffs.append(0.5 * Q[i, i])
    paulis.append(make_Z_string(i, nu_variables))

for i in range(nu_variables):
    for j in range(i + 1, nu_variables):
        if abs(Q[i, j]) > 1e-9:
            coeffs.append(0.25 * Q[i, j])
            paulis.append(make_ZZ_string(i, j, nu_variables))

H = SparsePauliOp.from_list(list(zip(paulis, coeffs)))
print(f"✅ Hamiltonian created with {len(H)} Pauli terms")

# ============================================================
# 4. Simulator and primitives
# ============================================================
backend = AerSimulator(method="matrix_product_state", device="CPU")
estimator = Estimator(run_options={"backend": backend})
sampler = Sampler(run_options={"backend": backend})

# ============================================================
# 5. QAOA ansatz + COBYLA optimisation
# ============================================================
reps = 20
ansatz = QAOAAnsatz(cost_operator=H, reps=reps)
print("QAOA circuit depth:", ansatz.decompose().depth())

def expected_energy(params):
    try:
        job = estimator.run([ansatz], [H], [params])
        return float(job.result().values[0])
    except Exception as e:
        print("Simulation failed:", e)
        return 1e6

best_fun, best_theta = np.inf, None
print("\n--- Multi-start COBYLA optimisation ---")
for attempt in range(3):  # reduce runs to save time
    theta0 = np.random.uniform(0, 2 * np.pi, ansatz.num_parameters)
    opt = COBYLA(maxiter=500)
    res = opt.minimize(fun=expected_energy, x0=theta0)
    print(f" Attempt {attempt+1:02d}: energy = {res.fun:.6f}")
    if res.fun < best_fun:
        best_fun, best_theta = res.fun, res.x

print(f"\n✅ Best energy: {best_fun:.6f}")
print("Optimal θ* =", np.round(best_theta, 4))

# ============================================================
# 6. Sampling
# ============================================================
qc = ansatz.assign_parameters(best_theta)
qc.measure_all()
result = sampler.run([qc], shots=5000).result()
qd = result.quasi_dists[0]
counts = qd.binary_probabilities()

# ============================================================
# 7. Feasibility + CSV output
# ============================================================
def qubo_energy(x, Q, c, offset=0.0):
    Qsym = Q + Q.T - np.diag(np.diag(Q))
    return float(x @ Qsym @ x + c @ x + offset)

rows = []
for bitstring, prob in counts.items():
    x = np.array(list(map(int, bitstring[::-1])))
    lhs = A @ x
    feas = True
    for i, s in enumerate(sense):
        s = s.strip()
        if s == "<=":
            ok = lhs[i] <= b[i] + 1e-6
        elif s == "=":
            ok = np.isclose(lhs[i], b[i], atol=1e-6)
        elif s == ">=":
            ok = lhs[i] >= b[i] - 1e-6
        else:
            ok = False
        feas &= ok
    E = qubo_energy(x, Q, c, offset)
    # prepend apostrophe so Excel treats as text
    bit_txt = "'" + bitstring
    rows.append({
        "bitstring (x1→x15)": bit_txt,
        "probability": prob,
        "sum_x": int(x.sum()),
        "energy": E,
        "feasible": feas
    })

df = pd.DataFrame(rows).sort_values("probability", ascending=False)
print(f"\n✅ Feasible bitstrings: {df['feasible'].sum()} / {len(df)}")

out = f"LB5_QAOA_Unbalanced_{datetime.now().strftime('%Y%m%d_%H%M')}.csv"
df.to_csv(out, index=False)
print(f"✅ Results saved to {out} (bitstrings forced as text for Excel)")

# ============================================================
# 8. Diagnostic plot — probability vs. constraint violation
# ============================================================
import matplotlib.pyplot as plt

violations = []
for bitstring in df["bitstring (x1→x15)"]:
    x = np.array(list(map(int, bitstring.strip("'")[::-1])))  # remove ' and reverse
    lhs = A @ x
    v = 0.0
    for i, s in enumerate(sense):
        s = s.strip()
        if s == "<=":
            v += max(0, lhs[i] - b[i])
        elif s == "=":
            v += abs(lhs[i] - b[i])
        elif s == ">=":
            v += max(0, b[i] - lhs[i])
    violations.append(v)

df["violation_sum"] = violations

plt.figure(figsize=(7,5))
plt.scatter(df["violation_sum"], df["probability"], alpha=0.6, s=15)
plt.yscale("log")
plt.xlabel("Total constraint violation (sum of magnitudes)")
plt.ylabel("Probability (log scale)")
plt.title("Constraint violation vs. probability of sampled states")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()