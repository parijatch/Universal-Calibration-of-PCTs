# Universal_Calibration_of_PCTs

This repository contains the code used to generate the simulation results presented in [1].

## Repository structure

- **`calibration heavy-tailed tests/`**
  Contains code for simulating p-values from a multivariate *t* copula with zero location parameter. These simulations illustrate that the **Pareto Combination Test (PCT)** achieves universal (asymptotic) calibration, while the **Cauchy Combination Test (CCT)** becomes conservative under stronger tail dependence.

- **`power simulations/`**
  Contains code for simulating p-values from a multivariate *t* copula with a non-null location parameter. The simulations demonstrate that PCT is uniformly at least as powerful as CCT, and strictly more powerful in most settings considered.

- **`data splitting and FCT/`**
  Contains code implementing the **Fréchet Combination Test (FCT)**. Simulation results in this folder verify that FCT is asymptotically calibrated only under tail independence.

- **`nhanes/`**
  Contains code used to illustrate the PCT using NHANES data.

## References

[1] Chakraborty, P., Guo, F. R., Shedden, K., & Stoev, S. (2025).
    On the universal calibration of Pareto-type linear combination tests.
    arXiv preprint arXiv:2509.12066.
