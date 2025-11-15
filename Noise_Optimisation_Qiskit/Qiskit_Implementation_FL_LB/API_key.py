#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 19:40:28 2025

@author: luismoncayo
"""
import qiskit

print(qiskit.__version__)
#bkxGF7RiEdxenWrwixTt4qORofLAYzH8sHrX8GpbfNVh

# from qiskit_ibm_runtime import QiskitRuntimeService

# service = QiskitRuntimeService(channel="ibm_quantum_platform", token = "bkxGF7RiEdxenWrwixTt4qORofLAYzH8sHrX8GpbfNVh")

# backend = service.backend(name="ibm_brisbane")
# backend.num_qubits

from qiskit_ibm_runtime import QiskitRuntimeService
 
QiskitRuntimeService.save_account(
  token="bkxGF7RiEdxenWrwixTt4qORofLAYzH8sHrX8GpbfNVh", # Use the 44-character API_KEY you created and saved from the IBM Quantum Platform Home dashboard
  #name="<account-name>", # Optional
  instance="crn:v1:bluemix:public:quantum-computing:us-east:a/943b8b92cf0240d5b73592aaee7c63ff:92e116d9-8b6f-4273-a4df-041d5174209c::", # Optional
  set_as_default=True, # Optional
  overwrite=True, # Optional
)

from qiskit_ibm_runtime import QiskitRuntimeService
service = QiskitRuntimeService()
print(service.backends())


