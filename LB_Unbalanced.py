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
from optimisation.compute_solution import ComputeSolution
from problems.problems import ProblemsToSolve
import os

#/home/luismoncayo/Downloads/Balancing_QUBO/Data_Files/Data_Line_Balancing/toy_examples
path_to_instance = '/home/luismoncayo/Dropbox/Python/Balancing_QUBO_U/Data_Files/Data_Line_Balancing/toy_examples/instance_toy_n=11.txt'
load_data = UploadDataLineBalacing(path_to_instance)
tasks, tasks_times, precedence = load_data.upload_data()

nu_stations = 3
cycle_time = 45
nu_constraints = nu_stations+len(tasks)+len(precedence) 
nu_variables = nu_stations*len(tasks)

var_name = [f'x{i}{j}' for j in range(1, nu_stations+1) for i in range(1, len(tasks)+1)]
constraints_matrix = np.zeros((nu_constraints, nu_variables))


###############################################################################
#[betas,gammas,shots,lambda_1,lambda_2]
parameters = [#[15,15,5000,0.0001,0.00001],
              #[15, 15, 3000, 0.001, 0.0001],
              #[20, 20, 3000, 0.0001, 0.00001],
              #[15, 15, 5000, 0.001, 0.0001],
              #[20, 20, 5000, 0.0001, 0.00001],
              #[15, 15, 5000, 0.005, 0.0005],
              #[15, 15, 5000, 0.009, 0.0009],
              #[15, 15, 5000, 0.01, 0.001],
              #[15, 15, 5000, 0.1, 0.01],
              [20, 20, 5000, 1, 0.9]
              ] 

total_experiments = 2
for p in range(len(parameters)):
    par = parameters[p]
    summary = []
    
    for e in range(1,total_experiments+1):
        optimise = ComputeSolution()
        problem = ProblemsToSolve()
        betas = np.linspace(0, 1, par[0])[::-1]
        gammas = np.linspace(0, 1, par[1])
        shots = par[2] # Number of samples used
        
        ##############################################################################################
        ### Attempt using Quantum Approximate Optimization Algorithm (QAOA)
        ##############################################################################################
        nu_var, mdl = problem.balance_problem(nu_stations, cycle_time, tasks, tasks_times, precedence)
        lambda_1 = par[3]
        lambda_2 = par[4]
        nu_qubits = nu_var
        h_new, J_new = optimise.sol_qaoa(mdl, lambda_1, lambda_2)
        cpu_time, samples_qaoa = optimise.qaoa_circuit(gammas, betas, h_new, J_new, nu_qubits, shots)
        ##############################################################################################
        
        var_names = [var.name for var in mdl.iter_variables()]
        solutions = pd.DataFrame(samples_qaoa,columns=var_names)
        
        store_type = []
        for row in range(len(solutions)):
            row_feasibility = []
            stations = range(1,nu_stations+1)
            for j in stations:
                sum_stat = 0
                for i in tasks:
                    sum_stat = sum_stat + solutions.loc[row,f"x_{i}_{j}"]*tasks_times[i-1]
                if sum_stat <= cycle_time:
                    row_feasibility.append(1)
                else:
                    row_feasibility.append(0)
            for i in tasks:
                sum_tasks = 0
                for j in stations:
                    sum_tasks = sum_tasks + solutions.loc[row,f"x_{i}_{j}"]
                if sum_tasks == 1:
                    row_feasibility.append(1)
                else:
                    row_feasibility.append(0)
            for pre in precedence:
                a = pre[0]
                b = pre[1]
                sum_pre = 0
                for j in stations:
                    sum_pre = sum_pre + j*solutions.loc[row,f"x_{a}_{j}"] - j*solutions.loc[row,f"x_{b}_{j}"]
                if sum_pre <= 0:
                    row_feasibility.append(1)
                else:
                    row_feasibility.append(0)
                
            store_type.append(all(x == 1 for x in row_feasibility))
        
        solutions["Type"] = store_type
        folder_name_e = "Outputs_Balancing/Solutions/p_"+str(par)
        os.makedirs(folder_name_e, exist_ok=True)
       
        file_name = folder_name_e+"/qaoa_p_"+str(p)+"_e_"+str(e)+"_bal.xlsx"
        solutions.to_excel(file_name, index=False)
       
        true_sol = solutions[solutions["Type"] == True]
        this_results =[e,cpu_time,true_sol.shape[0]]
        summary.append(this_results)
    
    folder_name_s = "Outputs_Balancing/Summary"
    os.makedirs(folder_name_s, exist_ok=True)
    summary_df = pd.DataFrame(summary, columns=['Experiment','CPU_Time','Nu_Opt_Sol'])
    summary_df.to_excel("Outputs_Balancing/Summary/qaoa_"+str(par)+".xlsx", index=False)
       
        
    

    
