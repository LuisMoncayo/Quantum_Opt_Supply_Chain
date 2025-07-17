#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 13:13:15 2025

@author: luismoncayo
"""
from docplex.mp.model import Model
import numpy as np
import pandas as pd

class ProblemsToSolve:
    
    def __init__(self):
        print("The problem is created")
        
    def location_problem(self,obj_func,set_constraints):
        nu_locations = set_constraints.shape[0]
        #nu_locations = 25
        
        mdl = Model()
        x = mdl.binary_var_list(nu_locations, name="x")
        
        objective = mdl.sum(x[i] * obj_func[i] for i in range(nu_locations))
        mdl.minimize(objective)
        
        for j in range(nu_locations):
            mdl.add_constraint(mdl.sum(x[k] * set_constraints[j,k] for k in range(nu_locations)) >= 1)
        
        num_vars = mdl.number_of_variables
        print(mdl.export_as_lp_string())
        return mdl
    
    def balance_problem(self,nu_stations,cycle_time,tasks,tasks_times,precedence):
        # Define sets
        stations = range(1,nu_stations+1)
        # Create model
        mdl = Model(name="Task_Station_Assignment")

        # Create binary decision variables x_ij
        x = mdl.binary_var_dict(((t, s) for t in tasks for s in stations), name="x")
        y = mdl.binary_var_dict(stations, name="y")

        objective = mdl.sum(y[j] for j in stations)
        mdl.minimize(objective)

        for j in stations:#for stations
            mdl.add_constraint(sum(tasks_times[i-1]*x[i,j] for i in tasks)- cycle_time*y[j] <= 0, "sta_"+str(j))

        for i in tasks:#for tasks
            mdl.add_constraint(sum(x[i,j] for j in stations) == 1, "tas_"+str(i))

        for p in precedence:
            a = p[0]
            b = p[1]
            mdl.add_constraint(sum(j*x[a,j]-j*x[b,j] for j in stations) <= 0, "pre"+str(p))
        
        num_vars = mdl.number_of_variables
        #print(mdl.export_as_lp_string())
        return num_vars, mdl
    
    def summary_QT(self, all_samples, set_constraints):

        types_solution = []
        for s in range(all_samples.shape[0]):
            rhs = []
            binary_string = all_samples.iloc[s,0]
            array = [int(digit) for digit in binary_string]
            for c in range(set_constraints.shape[0]):
                result = sum(array[i] * set_constraints[c][i] for i in range(set_constraints.shape[1]-1))
                rhs.append(result)
            all_ones = all(x == 1 for x in rhs)
            types_solution.append(all_ones)
        all_samples['Type_Sol'] = types_solution
        df_sorted = all_samples.sort_values(by="Energy_Ising")
        df_sorted[df_sorted['Type_Sol'] == 'True']
        
        return df_sorted
    
    def summary_QAOA(self, all_samples, set_constraints):

        types_solution = []
        for s in range(all_samples.shape[0]):
            rhs = []
            binary_string = all_samples.iloc[s,0]
            array = [int(digit) for digit in binary_string]
            for c in range(set_constraints.shape[0]):
                result = sum(array[i] * set_constraints[c][i] for i in range(set_constraints.shape[0]))
                rhs.append(result)
            all_ones = all(x >= 1 for x in rhs)
            types_solution.append(all_ones)
        all_samples['Type_Sol'] = types_solution
        df_sorted = all_samples.sort_values(by="Values_Obj_Func")
        df_sorted[df_sorted['Type_Sol'] == 'True']
        
        return df_sorted
    
        
