# Benders Decomposition for Two-Stage Robust Optimization

This folder contains a reproduction of the **Benders Decomposition (Cutting Plane)** algorithm applied to a two-stage robust optimization problem, following the numerical baseline described in:

> **Zeng, B., & Zhao, L. (2013).** Solving two-stage robust optimization problems using a column-and-constraint generation method. *Operations Research Letters*, 41(5), 457-461.

## Problem Context
The implementation focuses on the **Location-Transportation Problem** under demand uncertainty. The problem is formulated as a min-max-min model where first-stage decisions (facility locations) are made before the realization of uncertain demand, followed by optimal second-stage transportation dispatch.

## Methodology
The algorithm utilizes the traditional Benders-style cutting plane method:
1. **Master Problem**: Solves for first-stage location variables and a lower bound of the recourse cost.
2. **Subproblem**: Identifies the worst-case demand scenario via the dual of the second-stage problem.
3. **Iteration**: Dynamically generates optimality cuts to refine the master problem until the gap converges.

## Implementation Details
- **Test Case**: Small-scale location-transportation network.
- **Solver**: Gurobi / CPLEX.
- **Key Outcome**: Demonstration of the convergence properties and limitations of the standard Benders method in robust frameworks.
