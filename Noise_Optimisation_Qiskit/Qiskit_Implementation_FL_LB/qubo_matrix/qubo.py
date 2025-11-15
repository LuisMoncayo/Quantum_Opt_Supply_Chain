#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 11:49:02 2025

@author: luismoncayo
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:01:11 2025

@author: luismoncayo
"""
import numpy as np

class QuboTransformation:
    
    def __init__ (self, set_obj_func, set_constraints):
        self.set_obj_func = set_obj_func
        self.set_constraints = set_constraints
        
    def compute_Q(self, lagrange):
        #print(self.set_obj_func)
        #print(self.set_constraints)
        #print(lagrange)
        
        #constraints  = len(self.set_constraints)
        variables = len(self.set_constraints[0])-1

        matrix = np.zeros((variables,variables))
        offset = 0

        for cons in self.set_constraints:
            #print(cons)
            for i in range(len(cons)):
                if cons[i] != 0:
                    if i == variables:
                        offset += cons[i]*cons[i]
                    else:
                        for j in range(i,len(cons)):
                            #value = 2*coefficient*cons[j]
                            if i == j:
                                matrix[i][i] += cons[i]*cons[i]
                            elif j == variables:
                                matrix[i][i] += 2*cons[i]*cons[j]
                            else:
                                matrix[i][j] += 2*cons[i]*cons[j]
                            #print(f"i={i}, j={j}")
                        #print(f"matrix = {matrix}")
                        #print("------")
            #print("++++++++++++++")
        #print(matrix)
        #print(offset)

        matrix = lagrange*matrix
        offset = lagrange*offset

        for e in range(len(self.set_obj_func[0])):
            matrix[e][e] += self.set_obj_func[0][e] 
            
        # for r in range(variables):
        #     for c in range(r+1,variables):
        #         matrix[r][c] = (matrix[r][c]+matrix[c][r])/2
        #         matrix[c][r] = matrix[r][c]
                #print(f"r={r}, c={c}")
        #print(matrix)
        return offset, matrix

        
        
        
        
    
    
