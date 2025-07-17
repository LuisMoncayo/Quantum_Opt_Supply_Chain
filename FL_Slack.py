"""
Created on Wed Feb 12 12:37:59 2025

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

#/Users/luismoncayo/Dropbox/Python/PennyLane_QUBO/Data_Files/uk-postcode-polygons-master/geojson/EX.geojson
#/home/luismoncayo/Dropbox/Python/PennyLane_QUBO/Data_Files/uk-postcode-polygons-master/geojson
geojson_data = gpd.read_file("/Users/luismoncayo/Dropbox/Python/PennyLane_QUBO/Data_Files/uk-postcode-polygons-master/geojson/EX.geojson")
codes = ['EX9','EX10','EX11','EX12','EX13']
#codes = ['EX5','EX9','EX10','EX11','EX12','EX13']
#codes = ['EX9','EX10','EX11','EX12','EX13','EX14','EX15','EX24']
#codes = ['EX1','EX2','EX3','EX4','EX5','EX6','EX7']
gdf = geojson_data[geojson_data['name'].apply(lambda x: x in codes)]

problem_data = LoadData(gdf)
problem_data.print_map()
constraints, variables, set_constraints = problem_data.compute_set_constriants()
# matrix_df = pd.DataFrame(set_constraints)
# matrix_df.to_excel("constraints_location.xlsx", index=False, header=False)

obj_func = np.concatenate((np.ones(gdf.shape[0], dtype=int), np.zeros((variables-gdf.shape[0]), dtype=int)))
set_obj_func = obj_func[np.newaxis, :]

matrix_computation = QuboTransformation(set_obj_func,set_constraints)


#[betas,gammas,shots, lagrange]
parameters = [[20,20,10,0.5],
              #[20,20,5000,1],
              #[20,20,5000,10],
              #[10,10,3000,1],
              #[20,20,3000,5],
              ]

for par in parameters:

    total_experiments = 1
    summary = []
    for e in range(1,total_experiments+1):
        folder_name_e = "Outputs_Balancing/Solutions/p_"+str(par)
        os.makedirs(folder_name_e, exist_ok=True)
        
        lagrange = par[3];
        offset, QT = matrix_computation.compute_Q(lagrange)
        # df = pd.DataFrame(QT)
        # df.to_excel('QT_location.xlsx', index=False, header=False)
        
        optimise = ComputeSolution()
        problem = ProblemsToSolve()
        betas = np.linspace(0, 1, par[0])[::-1]
        gammas = np.linspace(0, 1, par[1])
        shots = par[2]
        
        ###############################################################################
        ### Attempt transforming the problem into a QUBO model
        ###############################################################################
        n_qubits = len(QT)
        print(n_qubits)
        h, J, zoffset = optimise.from_Q_to_Ising(QT, offset)
        cpu_time, samples_QT = optimise.qaoa_circuit(gammas, betas, h, J, n_qubits, shots)
        
        #all_shots = pd.DataFrame({'Solution': [''.join(map(str, row)) for row in samples_QT[1]]})
        all_shots = pd.DataFrame({'Solution': [''.join(map(str, row)) for row in samples_QT]})
        all_shots['Energy_Ising'] = all_shots['Solution'].apply(lambda x: optimise.energy_Ising(x, h, J, zoffset))
        summary_all = problem.summary_QT(all_shots, set_constraints)
        
        file_name = folder_name_e+"/p_"+str(par)+"_e_"+str(e)+"_FL.xlsx"
        summary_all.to_excel(file_name, index=False)
    
                  
        
        feasibles_sol = summary_all[summary_all['Type_Sol'] == True]
        #number_fesiable = feasibles_sol.shape[0]
        summary.append([e,cpu_time,feasibles_sol.shape[0]])
        # summary_df = pd.DataFrame(summary, columns=['Experiment','CPU_Time','Nu_Opt_Sol'])
        # summary_df.to_excel("Outputs_Balancing/Summary/qaoa_"+str(par)+".xlsx", index=False)
        
    folder_name_s = "Outputs_Balancing/Summary"
    os.makedirs(folder_name_s, exist_ok=True)
    summary_df = pd.DataFrame(summary, columns=['Experiment','CPU_Time','Nu_Opt_Sol'])
    summary_df.to_excel("Outputs_Balancing/Summary/qaoa_"+str(par)+".xlsx", index=False)













































###############################################################################
### Attempt using Quantum Approximate Optimization Algorithm (QAOA)
###############################################################################
# mdl = problem.location_problem(obj_func, set_constraints)
# lambda_1 = 3.5
# lambda_2 = 2.8
# nu_qubits = set_constraints.shape[0]
# h_new, J_new = optimise.sol_qaoa(mdl, lambda_1, lambda_2)
# samples_qaoa = optimise.qaoa_circuit(gammas, betas, h_new, J_new, nu_qubits, shots)

# results_qaoa = pd.DataFrame({'Shots': range(1, len(samples_qaoa) + 1), 'Binary_String': [''.join(map(str, row)) for row in samples_qaoa]})
# count_results_qaoa = results_qaoa['Binary_String'].value_counts().reset_index()
# count_results_qaoa['Values_Obj_Func'] = count_results_qaoa['Binary_String'].apply(lambda x: sum(int(digit) for digit in x))

# sorted_samples = problem.summary_QAOA(count_results_qaoa, set_constraints)
# sorted_samples.to_excel("sorted_samples_qaoa.xlsx")

###############################################################################
### Graphs ----
###############################################################################
# to_plot_QT = results_QT
# to_plot_QT['Energy_Ising'] = to_plot_QT['Binary_String'].apply(lambda x: optimise.energy_Ising(x, h, J, zoffset))
# optimise.energy_plot(to_plot_QT,"QT")

# to_plot_QAOA = results_qaoa
# to_plot_QAOA['Values_Obj_Func'] = to_plot_QAOA['Binary_String'].apply(lambda x: sum(int(digit) for digit in x))
# optimise.energy_plot(to_plot_QAOA,"QAOA")























