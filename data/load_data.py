#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:28:05 2025

@author: luismoncayo
"""
import matplotlib.pyplot as plt
import pandas as pd
import re
pd.set_option('display.max_columns', None)  # Show all columns

class LoadData:
    
    def __init__(self,gdf):
        self.gdf = gdf
        self.gdf['Solution'] = [0 for _ in range(self.gdf.shape[0])]
        
    def print_map(self):
        #gdf = gpd.read_file('Data_Files/uk-postcode-polygons-master/geojson/EX.geojson')
    
        fig, ax = plt.subplots(figsize=(18, 16))
        self.gdf.plot(ax=ax, edgecolor='black',column="Solution",cmap="gist_yarg")
        for idx, row in self.gdf.iterrows():
            # Get the centroid of the polygon
            centroid = row['geometry'].centroid
            # Add the label at the centroid coordinates
            ax.text(centroid.x, centroid.y, row['name'], fontsize=16, ha='center', color='black')
            
    def compute_set_constriants(self):
        adjacency_results = []
        # Loop through each pair of polygons to check if they are adjacent
        for idx1, poly1 in self.gdf.iterrows():
            this_poly = [poly1['name']]
            for idx2, poly2 in self.gdf.iterrows():
                if idx1 != idx2:  # Skip self-comparison
                    is_adjacent = poly1['geometry'].touches(poly2['geometry'])
                    if is_adjacent:
                        this_poly.append(poly2['name'])
                        #adjacency_results.append((poly1['name'], poly2['name'], is_adjacent))
            adjacency_results.append(this_poly)
        print(adjacency_results)
        
        # Print adjacency results
        u = []
        constraints = []
        for result in adjacency_results:
            sorted_numbers = sorted(int(ex[2:]) for ex in result)
            u.append(len(sorted_numbers))
            constraints.append(sorted_numbers)
            print(sorted_numbers)
        
        lst = []# = 2^{k+1}-1 -------
        k = 0
        while True:
            val = 2**(k+1) - 1
            lst.append(val)
            if val >= max(u)-1:
                break
            k += 1
        # -----------------------------

        add_var = []
        for v in u:
            for ii, val in enumerate(lst):
                if val >= v-1:
                    add_var.append(ii)
                    break
        
        bin_variables = []
        for o in add_var:
            sequence = [2**i for i in range(o+1)]
            bin_variables.append(sequence)
        
        aux_df_var = pd.DataFrame({'Cons': constraints, 'Add_S': bin_variables})
        
        aux_df_var["EX_Value"]=[int(re.search(r"\d+", item).group(0)) if re.search(r"\d+", item) else None for item in self.gdf["name"]]
        #aux_df_var["EX_Value"] = self.gdf["name"].apply(lambda x: int(re.sub(r"\D", "", x)))
        #aux_df_var["EX_Value"] = [9,10,11,12,13,24]
        #print(self.gdf)
        
        all_indexed = []
        for con in constraints:
            index_conts =[]
            for code in con:
                i, j = aux_df_var.eq(code).stack().idxmax()
                index_conts.append(i)
            all_indexed.append(index_conts)
        
        aux_df_var["Index_Cons"] = all_indexed
        
        #print(aux_df_var)
        total_length = aux_df_var['Add_S'].str.len().sum()
        num_cons = len(aux_df_var)
        
        # Define the matrix size
        import numpy as np
        np.set_printoptions(threshold=np.inf)
        matrix = np.zeros((num_cons, num_cons), dtype=int)
        
        # Fill the matrix
        for row_idx, values in enumerate(aux_df_var["Index_Cons"], start=0):
            for col_idx in values:
                if col_idx < num_cons:
                    matrix[row_idx, col_idx] = 1
        
        aux_df_var["EX_Value"].sum()
        
        matrix_S = np.zeros((num_cons, total_length), dtype=int)
        
        # Fill the matrix
        col = 0
        for row_idx, values in enumerate(aux_df_var["Add_S"]):
            #print(f"row={row_idx} and value = {values}")
            for i in range(len(values)):
                matrix_S[row_idx,col] = -values[i]
                col = col+1
        
        matrix_A = np.concatenate((matrix, matrix_S), axis=1)
        set_constraints = np.hstack((matrix_A, np.full((matrix_A.shape[0], 1), -1)))

        
        rows, columns = set_constraints.shape
        print(f"The problem has {rows} constraints and {columns-1} variables")
        return rows, columns-1, set_constraints
