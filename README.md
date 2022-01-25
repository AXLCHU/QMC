# QMC

## Var & CVaR calculations by QMC

### Objectives: 
- Implementation of low-discrepancy sequences & different sampling methods
- Implementation of stochastic algorithms converging to VaR/CVaR
- Add Imporance Sampling procedure and adaptive confidence level

### Results:
- Better coverage for QRNG vs PRNG
- Variance reduction when using QRNG vs PRNG - depends on the sampling method
- (Huge) Variance reduction ratio when adding Importance Sampling procedure - in particular for high confidence level

### Example:

1) Compute estimate of the optimal IS parameters with N/100 (M=100000) using a stepwise constant confidence level
2) Estimate VaR & CVaR with those paramters with M iterations

### Use of randomized Halton & Sobol sequences vs MT19937
###
### Normal algo vs IS 

![QMC_VaR_2](https://user-images.githubusercontent.com/56386159/150958118-9ef0bee0-123c-4cde-81df-2d491d2a8a46.PNG)

### Calculation for a 1 Put option portfolio with underlying following a GBM

![QMC_VaR_3](https://user-images.githubusercontent.com/56386159/150955969-1aaf15e6-ce89-450e-84da-17f4557202ee.PNG)

### Variance reduction for IS procedure

![QMC_VaR_4](https://user-images.githubusercontent.com/56386159/150955987-9eacbcba-af3c-4c61-9538-2f87f913871f.PNG)
