# Universal-Calibration-of-PCTs

This repository contains the code used for simulations done in [1].

The folder named *calibration - heavy-tailed tests* contains the code used for simulating p-values from a multivariate t copula with location parameter as the origin and illustrating that PCT (Pareto Combination Test) retains universal (asymptotic) calibration whereas CCT (Cauchy Combination Test) is conservative for stronger tail-dependence.

The folder named *power simulations* contains code that simulates p-values from a multivariate t copula with non-null location parameter and illustrates that PCT is always at least as powerful as CCT, with it being more powerful in most cases.

The folder named *data splitting and FCT* contains the code for the Fr\'etchet combination test and verifies using simulations that FCT is asymptotically calibrated only under tain independence.

##References

[1] Chakraborty, P., Guo, F. R., Shedden, K., & Stoev, S. (2025).
    On the universal calibration of Pareto-type linear combination tests.
    arXiv preprint arXiv:2509.12066.
