# C&CG for Robust Energy and Reserve Dispatch (IEEE 24-Bus)

This folder implements the **Column-and-Constraint Generation (C&CG)** algorithm for a robust energy and reserve dispatch problem. The dispatch model is based on the framework proposed in:

> 1. **Zugno, M., & Conejo, A. J. (2015).** A robust optimization approach to energy and reserve dispatch in electricity markets. *European Journal of Operational Research*, 247(2), 659-671.
> 2. IEC 61400-1 Wind turbines Part 1: Design requirements IEC, Switzerland, 2005. [Online]. Available: https://dlbargh.ir/mbayat/46.pdf.
> 3. wind-turbine-models. "Vestas V112-3.3." https://en.wind-turbine-models.com/turbines/693-vestas-v112-3.3 (accessed 2025).
> 4. R. L. Burden and J. D. Faires, "Chapter 3: Interpolation and Polynomial Approximation," in Numerical Analysis, 9th ed. Australia: Cengage Learning, 2011.

## Problem Description
The model addresses the joint clearing of energy and reserve auctions in systems with high renewable penetration. Key features include:
- **Two-Stage Robust Optimization**: Hedging against wind power uncertainty.
- **CVaR**: A Stochastic Programming with a CVaR objective function work as the benchmark.
- **First-Stage**: Day-ahead commitments for energy production and reserve capacity.
- **Second-Stage**: Real-time re-dispatch and reserve deployment under worst-case wind scenarios.
- **Network Constraints**: DC power flow; Full implementation on the **IEEE 24-bus system**.
- **Out-of-Sample Comparison**: Real world wind data is employed to compare dispatch cost of Robust Optimization and Stochastic Programming.

## Algorithm: Column-and-Constraint Generation
Instead of the dual-based cutting plane method (Benders Decomposition), this project utilizes the **C&CG algorithm** (Zeng & Zhao, 2013). By explicitly introducing recourse variables and constraints associated with identified worst-case scenarios into the master problem, the C&CG method achieves:
- Faster convergence compared to Benders decomposition.
- Stronger lower bounds due to the primal-cut structure.

## Requirements
- MATLAB YALMIP
- Optimization Solver: Gurobi.

## Numerical Data
The simulation is grounded in real-world high-fidelity data from an offshore wind farm in the East China Sea.
### 1. Data Source and Location
- **Acquisition**: Data retrieved locally via the cdsapi Python interface.
- **Target Location**: 30°48′34.404″ N, 121°53′21.699″ E (East China Sea Offshore Wind Farm).
- **Dataset**: ERA5 hourly data on single levels (1940–present).
- **Parameters**: 10m $u$ and $v$ wind components are utilized to calculate the hourly wind speed.
### 2. Hub-Center Wind Speed
ERA5 provides wind speeds at a 10m reference height, the wind speed is extrapolated to the hub height ($H = 119$ m) of the Vestas V112-3.3 turbine using the Power Law:

$$
\displaystyle v_{\text{hub}} = v_0 \cdot \left( \frac{H_{\text{hub}}}{H_0} \right)^{\alpha}
$$

- **Wind Shear Exponent ($\alpha$)**: Set to 0.15 (representing neutral stability conditions), which provides a more realistic estimation for offshore environments than the standard IEC 0.2 recommendation [2].
- **Wind Turbine**: Vestas-v112-3.3 wind turbine parameter is used for an assuming wind farm.
- **Wind Speed**: Wind speed from dataset at 10 meters is transfered to the windspeed at the hubcenter based on simple exponential function.
-  **Power Curve**: Newton interpolation is used to fit the wind power curve based on ref data.
### 3. Turbine Power Curve 
The power output model is based on the Vestas V112-3.3 MW turbine [3]. To overcome the sparsity of manufacturer data points, the power curve is reconstructed as follows:
- **Interpolation**: Newton’s Divided-Difference Interpolation is applied to achieve a smooth power-speed mapping [4].
- **Numerical Stability**: To suppress oscillations (Runge's phenomenon) near the cut-in speed ($3$–$4$ m/s), a piecewise linear interpolation is adopted.
- **Operating Range**: Cut-in speed at $3$ m/s, with rated power output achieved at $14$ m/s and above.
