#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 11:40:26 2025

@author: luismoncayo
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 20:00:12 2025

@author: luismoncayo
"""
import numpy as np
np.set_printoptions(linewidth=500)
import pandas as pd
from data.data_line_balancing import UploadDataLineBalacing
from qubo_matrix.qubo import QuboTransformation
from optimisation.compute_solution import ComputeSolution
from problems.problems import ProblemsToSolve
import os

#/home/luismoncayo/Downloads/Balancing_QUBO/Data_Files/Data_Line_Balancing/toy_examples
path_to_instance = '/home/luismoncayo/Dropbox/Python/Balancing_QUBO_U/Data_Files/Data_Line_Balancing/toy_examples/instance_toy_n=11.txt'
load_data = UploadDataLineBalacing(path_to_instance)
tasks, tasks_times, precedence = load_data.upload_data()

nu_stations = 3
cycle_time = 7
nu_constraints = nu_stations+len(tasks)+len(precedence) 
nu_variables = nu_stations*len(tasks)

var_name = [f'x{i}{j}' for j in range(1, nu_stations+1) for i in range(1, len(tasks)+1)]
constraints_matrix = np.zeros((nu_constraints, nu_variables))

#-----------------------------------------------------------------------------
row = 0
for j in range(1,nu_stations+1):
    k=0
    for i in range(1,len(tasks)+1):
        col = var_name.index(f'x{i}{j}')
        constraints_matrix[row,col] = tasks_times[k]
        k = k+1
    row = row+1

for i in tasks:
    for j in range(1,nu_stations+1):
        col = var_name.index(f'x{i}{j}')
        constraints_matrix[row,col]=1
    row = row+1
    
for p in precedence:
    #print(f'i = {p[0]}, j = {p[1]}')
    for s in range(1, nu_stations+1):
        col_i = var_name.index(f'x{p[0]}{s}')
        constraints_matrix[row,col_i]=1
        #print(f'i = x{p[0]}{s}')
        col_j = var_name.index(f'x{p[1]}{s}')
        constraints_matrix[row,col_j]=-1
        #print(f'j = x{p[1]}{s}')
    row = row+1

#-----------------------------------------------------------------------------
def slack_binary(upper_bound):
    u = upper_bound#upper bound
    lst = []# = 2^{k+1}-1 -------
    k = 0
    while True:
        val = 2**(k+1) - 1
        lst.append(val)
        if val >= u:
            break
        k += 1
    return [2**i for i in range(len(lst))]
# -----------------------------------------------------------------------------

for s in range(nu_stations):
    slack_var = slack_binary(cycle_time-min(tasks_times))
    zero_matrix = np.zeros((nu_constraints,len(slack_var)))
    zero_matrix[s, :] = slack_var
    constraints_matrix = np.concatenate((constraints_matrix, zero_matrix), axis=1)
    for n in range(len(slack_var)):
        var_name.append(f's{s+1}{n+1}')

from_row = nu_stations+len(tasks)
to_row = len(precedence)
for p in range(from_row,from_row+to_row):
    slack_var = slack_binary(nu_stations-1)
    zero_matrix = np.zeros((nu_constraints,len(slack_var)))
    zero_matrix[p, :] = slack_var
    constraints_matrix = np.concatenate((constraints_matrix, zero_matrix), axis=1)
    for n in range(len(slack_var)):
        var_name.append(f's{p+1}{n+1}')

# See the constraint matrix
df = pd.DataFrame(constraints_matrix, columns=var_name)
df.to_excel('constrains_matrix.xlsx', index=False)

set_obj_func = np.zeros((1,constraints_matrix.shape[1]))
rhs = np.concatenate([np.full(nu_stations, -cycle_time), np.full(len(tasks), -1), np.zeros(len(precedence))])
set_constraints = np.hstack((constraints_matrix, rhs.reshape(-1,1)))

matrix_computation = QuboTransformation(set_obj_func,set_constraints)
lagrange = 10;
offset, QT = matrix_computation.compute_Q(lagrange)

parameters = [[10,10,200],
              [15,15,200]] #[betas,gammas,shots,lambda_1,lambda_2]

par = parameters[0]

optimise = ComputeSolution()
problem = ProblemsToSolve()
betas = np.linspace(0, 1, par[0])[::-1]
gammas = np.linspace(0, 1, par[1])
shots = par[2] # Number of samples used
###############################################################################
### Attempt transforming the problem into a QUBO model
###############################################################################
n_qubits = len(QT)
h, J, zoffset = optimise.from_Q_to_Ising(QT, offset)
samples_QT = optimise.qaoa_circuit(gammas, betas, h, J, n_qubits, shots)

results_QT = pd.DataFrame({'Shots': range(1, len(samples_QT) + 1), 'Binary_String': [''.join(map(str, row)) for row in samples_QT]})
count_results_QT = results_QT['Binary_String'].value_counts().reset_index()
count_results_QT['Energy_Ising'] = count_results_QT['Binary_String'].apply(lambda x: optimise.energy_Ising(x, h, J, zoffset))

sorted_samples_QT = problem.summary_QT(count_results_QT, set_constraints)
sorted_samples_QT.to_excel("sorted_samples_QT.xlsx") 
###############################################################################














