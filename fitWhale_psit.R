rm(list = ls())
library(Rcpp); library(RcppArmadillo); library(foreach)
library(tidyverse); library(geodist)

runID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# runID = 1

# directory and file names ----
fold = 'github/'

path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
ifelse(!dir.exists(path.fit), dir.create(path.fit, recursive = T), FALSE)

path.r = paste0(fold, '/src/RFtns.R')
path.cpp = paste0(fold, '/src/RcppFtns.cpp')

data = c('sim_etat_psi', 'sim_etat_psit')
datai = data[runID]

# datai = 'sim_etat_psi'
# datai = 'sim_etat_psit'

fiti = 'Whale_psit'
filename = paste0(path.fit, datai, '_', fiti, '.RData')

# =============================================================================-
# fit the model ----
# =============================================================================-

load(paste0(path.data, datai, '.RData'))

Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
distmat = geodist(coord) / 1000
dim(distmat)
max(distmat)

numWobs = dat.dist %>%
  group_by(day) %>%
  summarise(num = n())


# fixed values ----
shape_c = 9/4
scale_c = 10/7.5


## starting values ----
startvals_beta0 = log(numWobs$num / (cellarea*nrow(coord)))

startval = list(
  beta = c(startvals_beta0, rnorm(ncol(Xwhale[[1]]))),
  kappa_psi = 1,
  psi = foreach(i = 1:length(Xwhale), .combine = cbind) %do% {
    dummy = t(chol(exp(-distmat/phi_eta))) %*% rnorm(nrow(Xwhale[[1]]))
    dummy - mean(dummy)
  },
  c = rep(2.46, nrow(Tw))
)

source(path.r)

set.seed(2)

Rf_fitWhale_psit(
  dat.dist, dat.pam, Xwhale,
  dur_dist, f_dist_obs, f_dist, p_pam,
  coord, cellarea, pii, phi_psi, shape_c, scale_c,
  n.iter = 1000, n.keep = 500, n.save = 500, n.adap = 200,
  new.run = T,
  startval, path.cpp, filename)


# n.iter = 100; n.keep = 50; n.save = 50; n.adap = 200; new.run = T
