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

fiti = 'Zoop_etat'
filename = paste0(path.fit, datai, '_', fiti, '.RData')

# =============================================================================-
# fit the model ----
# =============================================================================-

load(paste0(path.data, datai, '.RData'))

Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
distmat = geodist(coord) / 1000
dim(distmat)
max(distmat)

dat.obl = dat.zoop %>%
  filter(type == 'obl') %>%
  dplyr::select(month, day, type, lon, lat, zoop)
logYobl = log(dat.obl$zoop)


# fixed values ----
sig2_obl = trpar$sig2_obl
tau2 = trpar$tau2


## starting values ----
startval_alpha0 = rnorm(nrow(Tz), mean(logYobl), sqrt(tau2))

startval = list(
  alpha0tilde = mean(startval_alpha0),
  alpha = c(startval_alpha0, rnorm(ncol(Xzoop[[1]]))),
  kappa_eta = 1,
  eta = foreach(i = 1:length(Xzoop), .combine = cbind) %do% {
    dummy = t(chol(exp(-distmat/phi_eta))) %*% rnorm(nrow(Xzoop[[1]]))
    dummy - mean(dummy)
  },
  lam0 = rnorm(1),
  lam1 = rnorm(1),
  sig2_sur = 1
)

source(path.r)

set.seed(2)

Rf_fitZoop_etat(
  dat.zoop, Xzoop,
  coord, cellarea, phi_eta, sig2_obl, tau2,
  n.iter = 1000, n.keep = 500, n.save = 500, n.adap = 200,
  new.run = T,
  startval, path.cpp, filename)


# n.iter = 100; n.keep = 50; n.save = 50; n.adap = 50; new.run = T


