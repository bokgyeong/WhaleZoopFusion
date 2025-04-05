# Supplemental Code for "Joint Spatiotemporal Modeling of Zooplankton and Whale Abundance in a Dynamic Marine Environment"
Authors: Bokgyeong Kang, Erin M. Schliep, Alan E. Gelfand, Christopher W. Clark, Christine A. Hudak, Charles A. Mayo, Ryan Schosberg, Tina M. Yack, and Robert S. Schick

Detailed instructions for implementing the proposed joint species distribution models, as well as for generating the tables and figures presented in the paper, are provided.

## Required packages
The code has been tested using R version 4.4.0, “Puppy Cup.” To ensure successful execution, the following R packages must be installed: 
`tigris`, `tidyverse`, `sf`, `sp`, `geodist`, `foreach`, `Rcpp`, `RcppArmadillo`, `xtable`, `rnaturalearth`, `egg`, `scales`

## Step 1. Simulate data
- `data_etat_psi.R`: Generate a dataset from the model (i)
- `data_etat_psit.R`: Generate a dataset from the model (ii)
- The simulated datasets are saved in the directory `/data`
- `plot_data.R`: Generate Figure 1 (The location of Cape Cod Bay, MA and locations of zooplankton collection sites), Figure 2 (Distance sampling and passive acoustic monitoring data), and Figure 7 (a) (Zooplankton and whale data). The figures are saved in the directory `/fig`

## Step 2. Fit the proposed models to the simulated datasets
- `fitJoint_etat_psi.R`: Fit the joint model (i) to the simulated datasets
- `fitJoint_etat_psit.R`: Fit the joint model (ii) to the simulated datasets
- `fitWhale_psit.R`: Fit a whale model (right side of Figure 4 in the manuscript) to the simulated datasets
- `fitZoop_etat.R`: Fit a zoop model (left side of Figure 4 in the manuscript) to the simulated datasets
- The resulting posterior samples for model parameters are saved in the directory `/fit`

## Step 3. Summarize the posterior samples
- Compute posterior mean abundance surfaces. Compute posterior samples for abundance, log-likelihood, and the number of observed whales and calls
- `sumJoint_etat_psi.R`: Compute the above quantities based on the results from the joint model (i)  
- `sumJoint_etat_psit.R`: Compute the above quantities based on the results from the joint model (ii)  
- `sumWhale_psit.R`: Compute the above quantities based on the results from the whale model  
- `sumZoop_etat.R`: Compute the above quantities based on the results from the zoop model
- The results are saved in the directory `/sum`

## Step 4. Compare models via loglikelihood and continous ranked probability score (CRPS)
- `table_modelselection.R`: Generate Table 5 (Posterior median estimate and 95% HPD interval for the loglikelihood function as well as
CRPS for whale abundance detected by distance sampling) for the joint models (i) and (ii)

## Step 5. Summarize inference results
- `table_sum.R`: Generate Table 6 (Posterior median estimates and 95% HPD intervals for model parameters) and Table 7 (Posterior median estimates and 95% HPD intervals for average zooplankton abundance and for total whale abundance)
- `plot_sum.R`: Generate Figure 7 (b) (Posterior mean estimates of the average zooplankton abundance surfaces) and Figure 7 (c) (Posterior mean estimates of the whale intensity surfaces). The figures are saved in the directory `/fig`
