# QMC

## Var & CVaR calculations by QMC

### Objectives: 
- Implementation of low-discrepancy sequences
- Implementation of different sampling methods
- Implementation of stochastic algorithms converging to VaR/CVaR
- Add Imporance Sampling procedure and adaptive confidence level

### Results:
- Better coverage for QRNG vs PRNG
- Variance reduction when using QRNG vs PRNG - depends on the sampling method
- Variance reduction when adding Importance Sampling procedure

### Example:

# Stochastic algorithm : normal vs Importance Sampling and adaptative confidence level
## - Compute estimate of the optimal IS parameters with N/100 (M=100000) using a stepwise constant confidence level
## - Estimate VaR & CVaR with those paramters with M iterations

![QMC_VaR_2](https://user-images.githubusercontent.com/56386159/150955943-eef14992-2d68-4e1d-9c53-cb2fc3b20d3d.PNG)

#

![QMC_VaR_3](https://user-images.githubusercontent.com/56386159/150955969-1aaf15e6-ce89-450e-84da-17f4557202ee.PNG)

#

![QMC_VaR_4](https://user-images.githubusercontent.com/56386159/150955987-9eacbcba-af3c-4c61-9538-2f87f913871f.PNG)
