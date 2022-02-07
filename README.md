# VaR & CVaR calculations by QMC

### Objectives: 
- Implementation of low-discrepancy sequences & different sampling methods
- Implementation of stochastic algorithms converging to VaR & CVaR
- Add recursive Importance Sampling procedure with adaptive confidence level

#
### Initial algorithm:

- Variance reduction when using QRNG vs PRNG 

![QMC_VaR_5](https://user-images.githubusercontent.com/56386159/150991822-d6465847-2b43-4814-9c77-1a72efd961d7.PNG)

#
### Initial algo vs IS procedure:

- Variance reduction when adding Importance Sampling & adaptive confidence level (using randomized Halton sequence)

![QMC_VaR_2](https://user-images.githubusercontent.com/56386159/150958118-9ef0bee0-123c-4cde-81df-2d491d2a8a46.PNG)

#
### Example 1: 1 Put option portfolio:

- Using randomized Sobol sequence & GBM for underlying asset
- Important Variance Reduction ratios between the 2 procedures: +100% for CVaR & confidence level close to 1

![QMC_VaR_4](https://user-images.githubusercontent.com/56386159/150955987-9eacbcba-af3c-4c61-9538-2f87f913871f.PNG)

#
### Example 2: 10 Calls & 10 Puts on each of 5 uncorrelated GBM assets

- Variance Reduction of 40% for CVaR with nbr of simulaion = 1E5 & confidence level = 0.995

