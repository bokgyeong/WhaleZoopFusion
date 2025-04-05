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

fiti = 'Joint_etat_psit'
filename = paste0(path.fit, datai, '_', fiti, '.RData')

# =============================================================================-
# fit the model ----
# =============================================================================-

load(paste0(path.data, datai, '.RData'))

Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
distmat = geodist(coord) / 1000
dim(distmat)

dat.obl = dat.zoop %>%
  filter(type == 'obl') %>%
  dplyr::select(month, day, type, lon, lat, zoop)
logYobl = log(dat.obl$zoop)

numWobs = dat.dist %>%
  group_by(day) %>%
  summarise(num = n())


# fixed values ----
sig2_obl = trpar$sig2_obl
tau2 = trpar$tau2
shape_c = 9/4
scale_c = 10/7.5



## starting values ----
startval_alpha0 = rnorm(nrow(Tz), mean(logYobl), sqrt(tau2))
startvals_beta0 = log(numWobs$num / (cellarea*nrow(coord)))
startval_beta3 = 1.5

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
  sig2_sur = 1,
  beta = c(startvals_beta0, rnorm(ncol(Xwhale[[1]])), startval_beta3), ### fit ----
  kappa_psi = 1,
  psi = foreach(i = 1:length(Xwhale), .combine = cbind) %do% {
    dummy = t(chol(exp(-distmat/phi_psi))) %*% rnorm(nrow(Xwhale[[1]]))
    dummy - mean(dummy)
  },
  c = rep(2.46, nrow(Tw))
)


source(path.r)

set.seed(2)

Rf_fitJoint_etat_psit(
  dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
  dur_dist, f_dist_obs, f_dist, p_pam,
  coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, shape_c, scale_c,
  n.iter = 1000, n.keep = 500, n.save = 500, n.adap = 200,
  new.run = T,
  startval, path.cpp, filename)


