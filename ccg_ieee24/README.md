# C&CG for Robust Energy and Reserve Dispatch (IEEE 24-Bus)

This folder implements the **Column-and-Constraint Generation (C&CG)** algorithm for a robust energy and reserve dispatch problem. The dispatch model is based on the framework proposed in:

> **Zugno, M., & Conejo, A. J. (2015).** A robust optimization approach to energy and reserve dispatch in electricity markets. *European Journal of Operational Research*, 247(2), 659-671.

## Problem Description
The model addresses the joint clearing of energy and reserve auctions in systems with high renewable penetration. Key features include:
- **Two-Stage Robust Framework**: Hedging against wind power uncertainty.
- **First-Stage**: Day-ahead commitments for energy production and reserve capacity.
- **Second-Stage**: Real-time re-dispatch and reserve deployment under worst-case wind scenarios.
- **Network Constraints**: Full implementation on the **IEEE 24-bus system**.

## Algorithm: Column-and-Constraint Generation
Instead of the traditional dual-based cutting plane method, this project utilizes the **C&CG algorithm** (Zeng & Zhao, 2013). By explicitly introducing recourse variables and constraints associated with identified worst-case scenarios into the master problem, the C&CG method achieves:
- Faster convergence compared to Benders decomposition.
- Stronger lower bounds due to the primal-cut structure.

## Requirements
- MATLAB YALMIP
- Optimization Solver: Gurobi
