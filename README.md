# Supplemental Code for "Joint Spatiotemporal Modeling of Zooplankton and Whale Abundance in a Dynamic Marine Environment"
Authors: Bokgyeong Kang, Erin M. Schliep, Alan E. Gelfand, Christopher W. Clark, Christine A. Hudak, Charles A. Mayo, Ryan Schosberg, Tina M. Yack, and Robert S. Schick

Detailed instructions for implementing the proposed joint species distribution models, as well as for generating the tables and figures presented in the paper, are provided.

## Required packages
The code has been tested using R version 4.4.0, “Puppy Cup.” To ensure successful execution, the following R packages must be installed: 
`tigris`, `tidyverse`, `sf`, `sp`, `geodist`, `foreach`


## Simulation study

### Step 1. Simulate data
- `/sim/dataNHPP.R`: Generate a dataset from the model (i) NHPP
- `/sim/dataNHPPSE.R`: Generate a dataset from the model (ii) NHPP+E
- `/sim/dataLGCP.R`: Generate a dataset from the model (iii) NHPP+GP
- `/sim/dataLGCPSE.R`: Generate a dataset from the model (iv) NHPP+GP+E
- The simulated datasets are saved in the directory `/sim/data`

### Step 2. Fit the multivariate Hawkes process models to the simulated datasets 
- `/sim/fitNHPP.R`: Fit the model (i) NHPP to the simulated datasets
- `/sim/fitNHPPSE.R`: Fit the model (ii) NHPP+E to the simulated datasets
- `/sim/fitLGCP.R`: Fit the model (iii) NHPP+GP to the simulated datasets
- `/sim/fitLGCPSE.R`: Fit the model (iv) NHPP+GP+E to the simulated datasets
- The resulting posterior samples for model parameters are saved in the directory `/sim/fit`

### Step 3. Compare models via deviance information criterion (DIC) 

#### Compute loglikelihood
- `/sim/loglikNHPP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (i) NHPP fitted to the simulated datasets
- `/sim/loglikNHPPSE.R`:Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (ii) NHPP+E fitted to the simulated datasets
- `/sim/loglikLGCP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iii) NHPP+GP fitted to the simulated datasets
- `/sim/loglikLGCPSE.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The resulting posterior samples for the loglikelihood  are saved in the directory `/sim/loglik`
- We suggest determining the burn-in period by examining the trace plot of the loglikelihood chain

#### Obtain DIC and summarize results
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/sim/sumDIC.R`: Evaluate DIC for each model and create Table 2 included in the paper

### Step 4. Assess model adequacy via random time change theorem (RTCT)

#### Compute the expected number of calls occurring between two consecutive events
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/sim/rtctNHPP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (i) NHPP fitted to the simulated datasets
- `/sim/rtctNHPPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (ii) NHPP+E fitted to the simulated datasets
- `/sim/rtctLGCP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iii) NHPP+GP fitted to the simulated datasets
- `/sim/rtctLGCPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The results are saved in the directory `/sim/rtct`

#### Summarize results
- `/sim/sumRTCT.R`: Calculate the posterior mean estimates of the order statistics $\{d_{(i)}^{*}\}$, along with their associated uncertainties. Generate Q-Q plots (Figure 4) and calculate the mean squared difference (Table 1) for models (i) to (iv) fitted to each of the simulated datasets

### Step 5. Perform inference using a compensator

#### Evaluate the expected number of calls
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/sim/numNHPP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each MARU for the model (i) NHPP fitted to the simulated datasets
- `/sim/numNHPPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each MARU for the model (ii) NHPP+E fitted to the simulated datasets
- `/sim/numLGCP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each MARU for the model (iii) NHPP+GP fitted to the simulated datasets
- `/sim/numLGCPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each MARU for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The results are saved in the directory `/sim/num`

#### Summarize results
- `/sim/sumNum.R`: Evaluate the expected counts for total calls, contact calls, and countercalls received across MARUs, and obtain the corresponding empirical posterior distributions (Figure 5)
