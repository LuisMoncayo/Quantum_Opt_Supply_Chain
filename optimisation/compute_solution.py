#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 11:22:17 2025

@author: luismoncayo
"""
import pennylane as qml
import numpy as np
from collections import defaultdict
import time
from openqaoa.problems import FromDocplex2IsingModel

class ComputeSolution:
    
    def __init__(self):
        print("Optimising using PennyLane")
        
    def from_Q_to_Ising(self, Q, offset):
        """Convert the matrix Q of Eq.3 into Eq.13 elements J and h"""
        n_qubits = len(Q)  # Get the number of qubits (variables) in the QUBO matrix
        # Create default dictionaries to store h and pairwise interactions J
        h = defaultdict(int)
        J = defaultdict(int)

        # Loop over each qubit (variable) in the QUBO matrix
        for i in range(n_qubits):
            # Update the magnetic field for qubit i based on its diagonal element in Q
            h[(i,)] -= Q[i, i] / 2
            # Update the offset based on the diagonal element in Q
            offset += Q[i, i] / 2
            # Loop over other qubits (variables) to calculate pairwise interactions
            for j in range(i + 1, n_qubits):
                # Update the pairwise interaction strength (J) between qubits i and j
                J[(i, j)] += Q[i, j] / 4
                # Update the magnetic fields for qubits i and j based on their interactions in Q
                h[(i,)] -= Q[i, j] / 4
                h[(j,)] -= Q[i, j] / 4
                # Update the offset based on the interaction strength between qubits i and j
                offset += Q[i, j] / 4
        # Return the magnetic fields, pairwise interactions, and the updated offset
        return h, J, offset
        
    def qaoa_circuit(self, gammas, betas, h, J, num_qubits, shots):
        #dev = qml.device("default.qubit", shots=shots)
        #dev = qml.device("lightning.qubit", wires=num_qubits, shots=shots)
        dev = qml.device("lightning.gpu", wires=num_qubits, shots=shots)
        @qml.qnode(dev)
        def circuit():
            wmax = max(
                np.max(np.abs(list(h.values()))), np.max(np.abs(list(J.values())))
            )  # Normalizing the Hamiltonian

            p = len(gammas)

            # Apply initial Hadamard gates
            for i in range(num_qubits):
                qml.Hadamard(wires=i)

            # Apply p layers of the QAOA circuit
            for layer in range(p):
                # ---------- COST HAMILTONIAN ----------
                for ki, v in h.items():  # Single-qubit terms
                    qml.RZ(2 * gammas[layer] * v / wmax, wires=ki[0])
                for kij, vij in J.items():  # Two-qubit terms
                    qml.CNOT(wires=[kij[0], kij[1]])
                    qml.RZ(2 * gammas[layer] * vij / wmax, wires=kij[1])
                    qml.CNOT(wires=[kij[0], kij[1]])

                # ---------- MIXER HAMILTONIAN ----------
                for i in range(num_qubits):
                    qml.RX(-2 * betas[layer], wires=i)
            
            return qml.sample()
        
        # Measure CPU time
        start_time = time.process_time()
        result = circuit()
        end_time = time.process_time()
        
        cpu_time = end_time - start_time
        #print(f"CPU time: {cpu_time:.6f} seconds")
        return cpu_time, result  # Call the qnode function
    
    def samples_dict(self, samples, n_items):
        """Just sorting the outputs in a dictionary"""
        results = defaultdict(int)
        for sample in samples:
            #print(f"Inside samples_dict {sample}")
            results["".join(str(i) for i in sample)[:n_items]] += 1
        return results
    
    def energy_Ising(self, z, h, J, offset):
        """
        Calculate the energy of an Ising model given spin configurations.
    
        Parameters:
        - z: A dictionary representing the spin configurations for each qubit.
        - h: A dictionary representing the magnetic fields for each qubit.
        - J: A dictionary representing the pairwise interactions between qubits.
        - offset: An offset value.
    
        Returns:
        - energy: The total energy of the Ising model.
        """
        if isinstance(z, str):
            z = [(1 if int(i) == 0 else -1) for i in z]
            #print(f"The z = {z}")
        # Initialize the energy with the offset term
        energy = offset
        # Loop over the magnetic fields (h) for each qubit and update the energy
        for k, v in h.items():
            energy += v * z[k[0]]
        # Loop over the pairwise interactions (J) between qubits and update the energy
        for k, v in J.items():
            energy += v * z[k[0]] * z[k[1]]
        # Return the total energy of the Ising model
        return energy
        
    def sol_qaoa(self,mdl,lam_1, lam_2):
        lambda_1, lambda_2 = (
            lam_1,
            lam_2,
        )  # Parameters of the unbalanced penalization function (They are in the main paper)
        ising_hamiltonian = FromDocplex2IsingModel(
            mdl,
            unbalanced_const=True,
            #strength_ineq=[lambda_1,lambda_2,0.1,0.1],  # https://arxiv.org/abs/2211.13914
            strength_ineq=[lambda_1,lambda_2],  # https://arxiv.org/abs/2211.13914
        ).ising_model

        h_new = {
            tuple(i): w for i, w in zip(ising_hamiltonian.terms, ising_hamiltonian.weights) if len(i) == 1
        }
        J_new = {
            tuple(i): w for i, w in zip(ising_hamiltonian.terms, ising_hamiltonian.weights) if len(i) == 2
        }
        
        return h_new, J_new
    