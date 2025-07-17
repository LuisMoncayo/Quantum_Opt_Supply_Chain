#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 12:52:10 2025

@author: luismoncayo
"""
import geopandas as gpd
import pandas as pd
pd.set_option('display.max_columns', None)  # Show all columns
import numpy as np
from data.load_data import LoadData
from qubo_matrix.qubo import QuboTransformation
from optimisation.compute_solution import ComputeSolution
from problems.problems import ProblemsToSolve
import os

geojson_data = gpd.read_file("/home/luismoncayo/Dropbox/Python/PennyLane_QUBO/Data_Files/uk-postcode-polygons-master/geojson/EX.geojson")
#codes = ['EX9','EX10','EX11','EX12','EX13']
#codes = ['EX5','EX9','EX10','EX11','EX12','EX13']
#codes = ['EX1','EX2','EX3','EX4','EX5','EX6','EX7']
codes = ['EX9','EX10','EX11','EX12','EX13','EX14','EX15','EX24']

gdf = geojson_data[geojson_data['name'].apply(lambda x: x in codes)]

problem_data = LoadData(gdf)
#problem_data.print_map()
constraints, variables, set_constraints = problem_data.compute_set_constriants()
# matrix_df = pd.DataFrame(set_constraints)
# matrix_df.to_excel("constraints_location.xlsx", index=False, header=False)

obj_func = np.concatenate((np.ones(gdf.shape[0], dtype=int), np.zeros((variables-gdf.shape[0]), dtype=int)))
set_obj_func = obj_func[np.newaxis, :]


###############################################################################
#[betas,gammas,shots,lambda_1,lambda_2]
parameters = [[10,10,10,1.0,0.5],
              #[30,30,8000,1.0,0.9],
              #[30,30,8000,1.5,1.0],
              #[30,30,8000,1.5,1.4],
              #[30,30,8000,0.5,0.01]
              ] 

total_experiments = 1
for par in parameters:

    summary = []
    for e in range(1,total_experiments+1):
        folder_name_e = "Outputs_Balancing/Solutions/p_"+str(par)
        os.makedirs(folder_name_e, exist_ok=True)
        
        optimise = ComputeSolution()
        problem = ProblemsToSolve()
        betas = np.linspace(0, 1, par[0])[::-1]
        gammas = np.linspace(0, 1, par[1])
        shots = par[2]
        
        ###############################################################################
        ### Attempt using Quantum Approximate Optimization Algorithm (QAOA)
        ###############################################################################
        mdl = problem.location_problem(obj_func, set_constraints)
        lambda_1 = par[3]
        lambda_2 = par[4]
        nu_qubits = set_constraints.shape[0]
        print(f"Qubits: {nu_qubits}" )
        h_new, J_new = optimise.sol_qaoa(mdl, lambda_1, lambda_2)
        cpu_time, samples_qaoa = optimise.qaoa_circuit(gammas, betas, h_new, J_new, nu_qubits, shots)
        ###############################################################################
        
        all_shots = pd.DataFrame({'Solution': [''.join(map(str, row)) for row in samples_qaoa]})
        all_shots['Values_Obj_Func'] = all_shots['Solution'].apply(lambda x: sum(int(digit) for digit in x))
        summary_all = problem.summary_QAOA(all_shots, set_constraints)# <<-----------
        
        file_name = folder_name_e+"/U_p_"+str(par)+"_e_"+str(e)+"_FL.xlsx"
        summary_all.to_excel(file_name, index=False)
    
        feasibles_sol = summary_all[summary_all['Type_Sol'] == True]
        #number_fesiable = feasibles_sol.shape[0]
        summary.append([e,cpu_time,feasibles_sol.shape[0]])
        # summary_df = pd.DataFrame(summary, columns=['Experiment','CPU_Time','Nu_Opt_Sol'])
        # summary_df.to_excel("Outputs_Balancing/Summary/qaoa_"+str(par)+".xlsx", index=False)
        
    folder_name_s = "Outputs_Balancing/Summary"
    os.makedirs(folder_name_s, exist_ok=True)
    summary_df = pd.DataFrame(summary, columns=['Experiment','CPU_Time','Nu_Opt_Sol'])
    summary_df.to_excel("Outputs_Balancing/Summary/U_qaoa_"+str(par)+".xlsx", index=False)



