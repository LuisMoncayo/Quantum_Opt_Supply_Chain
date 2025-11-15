#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 10:01:29 2025

@author: luismoncayo
"""
"""
Line Balancing QUBO — full sample export
Author: Luis Moncayo
"""

import pandas as pd
import numpy as np
from datetime import datetime
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.converters import QuadraticProgramToQubo

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
A = df.iloc[6:6 + nu_constraints, 0:nu_variables].to_numpy()
sense = df.iloc[6:6 + nu_constraints, nu_variables].astype(str).values
b = df.iloc[6:6 + nu_constraints, nu_variables + 1].astype(float).values

print(f"--- Dimensions ---\nA: {A.shape}, sense: {sense.shape}, b: {b.shape}")

# ============================================================
# 2. QuadraticProgram → QUBO
# ============================================================
qp = QuadraticProgram()
for j in range(nu_variables):
    qp.binary_var(name=f"x{j+1}")

qp.minimize(linear={f"x{j+1}": 0.0 for j in range(nu_variables)})

for i in range(nu_constraints):
    lin_expr = {f"x{j+1}": float(A[i, j]) for j in range(nu_variables) if A[i, j] != 0}
    qp.linear_constraint(linear=lin_expr, sense=sense[i], rhs=float(b[i]), name=f"c{i}")

converter = QuadraticProgramToQubo()
qubo = converter.convert(qp)

Q = qubo.objective.quadratic.to_array(symmetric=True)
c = qubo.objective.linear.to_array()
offset = qubo.objective.constant
n_qubo = Q.shape[0]

print(f"\n✅ QUBO built successfully.\n  Q shape: {Q.shape}\n  Linear term: {c.shape}")

# ============================================================
# 3. Multi-start classical sampling
# ============================================================
N_ATTEMPTS = 5000
rows = []
best_energy = np.inf
best_x = None

for attempt in range(N_ATTEMPTS):
    x = np.random.randint(0, 2, n_qubo)
    energy = float(x @ Q @ x + c @ x + offset)
    if energy < best_energy:
        best_energy, best_x = energy, x

    # Map back to original 15 vars
    x_real = x[:nu_variables]
    lhs = np.array(A @ x_real, dtype=float)
    feasible = True
    for i, s in enumerate(sense):
        s = s.strip()
        if s == "<=":
            ok = lhs[i] <= b[i] + 1e-9
        elif s == "=":
            ok = np.isclose(lhs[i], b[i])
        elif s == ">=":
            ok = lhs[i] >= b[i] - 1e-9
        else:
            ok = False
        feasible &= ok

    rows.append({
        "attempt": attempt + 1,
        "bitstring_full": "".join(map(str, x.tolist())),
        "bitstring_real": "".join(map(str, x_real.tolist())),
        "energy": energy,
        "sum_x": int(x_real.sum()),
        "feasible": feasible
    })

    if attempt % 1000 == 0:
        print(f"  ... {attempt}/{N_ATTEMPTS}, best energy = {best_energy:.3f}")

# ============================================================
# 4. Results summary
# ============================================================
df_results = pd.DataFrame(rows)
n_feas = df_results["feasible"].sum()
best_row = df_results.loc[df_results["energy"].idxmin()]

print("\n--- Summary ---")
print(f"Feasible samples: {n_feas} / {len(df_results)}")
print(f"Best energy: {best_row['energy']:.6f}")
print(f"Best bitstring (real): {best_row['bitstring_real']}")

# ============================================================
# 5. Export all results
# ============================================================
fname = f"LB5_QUBO_all_samples_{datetime.now().strftime('%Y%m%d_%H%M')}.csv"
df_results.to_csv(fname, index=False)
print(f"\n✅ All {len(df_results)} samples saved to {fname}")




