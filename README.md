# CCG-Benders-PowerDispatch

Independent reproductions of advanced decomposition algorithms (C&CG, Benders) for network-constrained power system dispatch and operations. 

This repository serves as a numerical validation toolkit for scalable bi-level model of the sequential decision in the European electricity market.

## 📂 Repository Structure

```text
├── models/             # Core mathematical formulations (MATLAB + Gurobi)
├── algorithms/         # Implementations of decomposition methods (C&CG, Benders)
├── data/               # Network parameters and standard test cases (e.g., IEEE 33/69-bus)
├── examples/           # MATLAB Live Scripts (.mlx) and main scripts demonstrating step-by-step execution
├── README.md           # Project overview, methodology, and setup instructions
└── LICENSE             # MIT License
```
## ⚙️ Environment & Dependencies
To execute the models and decomposition algorithms in this repository, the following environment is required:
* **MATLAB** (R2021b or later recommended)
* **Gurobi Optimizer** (v11.0 or later)
* **YALMIP Toolbox** 
