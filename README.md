# Var & CVaR calculations by QMC

### Objectives: 
- Implementation of low-discrepancy sequences & different sampling methods
- Implementation of stochastic algorithms converging to VaR & CVaR
- Add recursive Importance Sampling procedure with adaptive confidence level

### Initial algorithm:

- Variance reduction when using randomized Halton & Sobol sequences vs PRNG 

![QMC_VaR_5](https://user-images.githubusercontent.com/56386159/150991822-d6465847-2b43-4814-9c77-1a72efd961d7.PNG)


### Initial algo vs IS procedure with randomized Halton sequence:

- Variance reduction when adding Importance Sampling & adaptive confidence level

![QMC_VaR_2](https://user-images.githubusercontent.com/56386159/150958118-9ef0bee0-123c-4cde-81df-2d491d2a8a46.PNG)


### Example for a 1 Put option portfolio:

- Using randomized Sobol sequence & GBM for underlying asset

![QMC_VaR_4](https://user-images.githubusercontent.com/56386159/150955987-9eacbcba-af3c-4c61-9538-2f87f913871f.PNG)


### Variance Reduction ratios between the 2 procedures:

- Huge Variance Reduction ratios: in particular for CVaR & high confidence level


