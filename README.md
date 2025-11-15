Quantum Optimisation for Supply Chain
QUBO Formulations and QAOA Implementations

This repository provides the full implementation of quantum-inspired optimisation methods for classical supply-chain problems. It includes QUBO formulations and QAOA pipelines for Facility Location (FL) and Load Balancing (LB), implemented using PennyLane and Qiskit. Both noiseless and noise-aware simulations are supported, including experiments calibrated with real IBM Quantum backends.

Contents

QUBO models for FL and LB (slack-variable and unbalanced penalisation versions)

QAOA implementations in PennyLane and Qiskit

Noise-model simulations using AerSimulator and hardware-derived calibration data

Statistical analysis tools for evaluating feasible solutions, energy distributions, and convergence

Data loaders for the LB instances and parameter files

Plots and visualisation scripts for publication-ready figures

Utilities for reproducibility, logging, and multi-run experiments

Features

Hardware-aligned noisy simulations (IBM backends such as ibm_brisbane)

Classical optimisation loops (COBYLA, SPSA)

Full reproducibility through deterministic seeding and saved output files

Structured pipeline for generating, running, and analysing QAOA experiments

Purpose

This repository accompanies the research work on applying QAOA to foundational supply-chain problems, providing a transparent and open implementation suitable for benchmarking, teaching, and reproducible research.
