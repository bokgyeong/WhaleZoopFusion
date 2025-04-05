rm(list = ls())
require(tigris) # counties
require(tidyverse)
require(foreach)
require(xtable)
require(coda) # HPDinterval


# directory and file names ----
fold = 'github/'

path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
path.sum = paste0(fold, 'sum/')


datai = 'sim_etat_psi'
# datai = 'sim_etat_psit'

fiti = 'Joint_etat_psi'
# fiti = 'Joint_etat_psit'


# =============================================================================-
# load data ----
# =============================================================================-

load(paste0(path.data, datai, '.RData'))

## observation days ----

date.zoop = format(unique(dat.zoop$date), '%b %d')
date.whale = dat.dist %>% 
  dplyr::select(month, day) %>% 
  left_join(dat.zoop %>% dplyr::select(date, month, day)) %>% 
  dplyr::select(date) %>% 
  unique %>% 
  mutate(date = format(unique(date), '%b %d')) %>% 
  unlist() %>% 
  as.vector



## zoop parameters ----
cname.zoop = c(
  'alphaTilde_0', paste0('alpha_0,', date.zoop), paste0('alpha_', c('temp')), 
  paste0('lam', 0:1), 'sig2_sur', 'kappa_eta'
)
trp.zoop = data.frame(
  Parameter = factor(cname.zoop, levels = cname.zoop), 
  Value = c(
    trpar$alpha0tilde, trpar$alpha0 + colMeans(trpar$eta), 
    trpar$alpha, trpar$lam0, trpar$lam1, trpar$sig2_sur, trpar$kappa_eta
  )
)


## whale parameters ----
cname.whale = c(
  paste0('beta_0,', date.whale), paste0('beta_', c('bath', 'zoop')), 
  paste0('c_', date.whale), 'kappa_psi'
)
trp.whale = data.frame(
  Parameter = factor(cname.whale, levels = cname.whale), 
  Value = c(trpar$beta0 + mean(trpar$psi), trpar$beta, trpar$c, trpar$kappa_psi)
)



# =============================================================================-
# posterior samples ----
# =============================================================================-

load(paste0(path.fit, datai, '_', fiti, '.RData'))

curniter
nrow(postAlpha)
n.burn = curniter - nrow(postAlpha)

## zoop ----
postZoop = data.frame(
  cbind(
    n.burn+1:nrow(postAlpha),
    postAlpha0tilde, postAlpha, postLam, postSig2, postKappa[,1]
  )
)

## whale ----
postWhale = data.frame(
  cbind(n.burn+1:nrow(postAlpha), postBeta, postC, postKappa[,2])
)

colnames(postZoop) = c('Iteration', cname.zoop)
colnames(postWhale) = c('Iteration', cname.whale)

postSamp = cbind(postWhale, postZoop)




# =============================================================================-
# abundance estimates ----
# =============================================================================-

load(paste0(path.sum, datai, '_', fiti, '_sum.RData'))

## whale ----
postWabun = foreach(t = 1:ncol(post.n.whale), .combine = rbind) %do% {
  data.frame(
    Date = date.whale[t],
    Iteration = 1:nrow(post.n.whale),
    Value = post.n.whale[,t]
  )
}

## zoop ----
postZabun = foreach(t = 1:ncol(post.n.zoop), .combine = rbind) %do% {
  data.frame(
    Date = date.zoop[t],
    Iteration = 1:nrow(post.n.zoop),
    Value = post.n.zoop[,t]
  )
}

head(postZabun)
head(postWabun)

postZabun$Date = factor(postZabun$Date, levels = unique(postZabun$Date))
postWabun$Date = factor(postWabun$Date, levels = unique(postWabun$Date))



# =============================================================================-
# table for parameters ----
# =============================================================================-

# 
# require(rstan) # for Rhat
# require(xtable)
# require(batchmeans)


## whale ----
par.notations.w = c(
  paste0('$\\beta_{', c('bath', 'zoop'), '}$'), paste0('$c_{', date.whale, '}$'), 
  '$\\kappa_{\\psi}$'
)  
par.ind.w = c('beta_bath', 'beta_zoop','c_Feb 24', 'c_Mar 23', 'c_Apr 29', 'kappa_psi')

tab.whale.par = postSamp %>% 
  dplyr::select(
    beta_bath, beta_zoop,
    `c_Feb 24`, `c_Mar 23`, `c_Apr 29`,
    kappa_psi
  ) %>% 
  pivot_longer(
    cols = everything(),
    values_to = 'Value',
    names_to = 'Parameter'
  ) %>% 
  mutate(Parameter = factor(Parameter, levels = unique(Parameter))) %>% 
  group_by(Parameter) %>% 
  summarise(
    Median = format(round(median(Value), 2), nsmall = 2),
    `95% HPD` = paste0(
      '(', format(round(HPDinterval(as.mcmc(Value))[1], 2), nsmall = 2), ', ',
      format(round(HPDinterval(as.mcmc(Value))[2], 2), nsmall = 2), ')'
    )
  ) %>% 
  mutate(Parameter = par.notations.w) %>% 
  dplyr::select(
    Parameter, Median, `95% HPD`,
  ) 

tab.whale.par %>% 
  xtable() %>% 
  print(include.rownames = F, sanitize.text.function = function(x) {x})
  
  

## zoop ----
par.notations.z = c(
  '$\\alpha_{temp}$', paste0('$\\lambda^{sur}_{', 0:1, '}$'), '$\\sigma^2_{sur}$', '$\\kappa_{\\eta}$'
)  
par.ind.z = c('alpha_temp', 'lam0', 'lam1', 'sig2_sur', 'kappa_eta')


tab.zoop.par = postSamp %>% 
  dplyr::select(
    alpha_temp, lam0, lam1, sig2_sur, kappa_eta
  ) %>% 
  pivot_longer(
    cols = everything(),
    values_to = 'Value',
    names_to = 'Parameter'
  ) %>% 
  mutate(Parameter = factor(Parameter, levels = unique(Parameter))) %>% 
  group_by(Parameter) %>% 
  summarise(
    Median = format(round(median(Value), 2), nsmall = 2),
    `95% HPD` = paste0(
      '(', format(round(HPDinterval(as.mcmc(Value))[1], 2), nsmall = 2), ', ',
      format(round(HPDinterval(as.mcmc(Value))[2], 2), nsmall = 2), ')'
    )
  ) %>% 
  mutate(Parameter = par.notations.z) %>% 
  dplyr::select(
    Parameter, Median, `95% HPD`,
  ) 

tab.zoop.par %>% 
  xtable() %>% 
  print(include.rownames = F, sanitize.text.function = function(x) {x})



## combine----

cbind(
  tab.zoop.par %>%
    add_row(Parameter = NA, Median = NA, `95% HPD`  = NA),
  tab.whale.par
) %>% 
  xtable() %>% 
  print(include.rownames = F, sanitize.text.function = function(x) {x})



# =============================================================================-
# table for abundance ----
# =============================================================================-


## whale ----

tab.whale.abun = postWabun %>% 
  group_by(Date) %>% 
  summarise(
    Median = format(round(median(Value)), nsmall = 0),
    HPD = paste0(
      '(', format(round(HPDinterval(as.mcmc(Value))[1]), nsmall = 0),
      ', ', format(round(HPDinterval(as.mcmc(Value))[2]), nsmall = 0), ')') 
  ) 

tab.whale.abun %>% 
  xtable() %>%
  print(include.rownames = F)



## zoop ----

tab.zoop.abun = postZabun %>% 
  group_by(Date) %>% 
  summarise(
    Median = format(round(median(Value)), nsmall = 0),
    HPD = paste0(
      '(', format(round(HPDinterval(as.mcmc(Value))[1]), nsmall = 0),
      ', ', format(round(HPDinterval(as.mcmc(Value))[2]), nsmall = 0), ')') 
  ) 

tab.zoop.abun %>% 
  xtable() %>%
  print(include.rownames = F)


## combine ----

cbind(
  tab.zoop.abun %>%
    add_column(Empty = ''),
  rbind(
    data.frame(Median = NA, HPD = NA),
    data.frame(Median = NA, HPD = NA),
    tab.whale.abun %>% 
      ungroup %>% 
      select(Median, HPD),
    data.frame(Median = NA, HPD = NA)
  )
) %>% 
  xtable() %>% 
  print(include.rownames = F, sanitize.text.function = function(x) {x})

