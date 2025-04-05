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

sim_design = rbind(
  data.frame(
    Model = 'Joint_etat_psi',
    ModelName = '(i)'
  ),
  data.frame(
    Model = 'Joint_etat_psit',
    ModelName = '(ii)'
  )
)


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


## number of whales and calls ----
obs.n.whale.dist = dat.dist %>% 
  group_by(day) %>% 
  summarise(n = length(lon)) %>% 
  ungroup %>% select(n) %>% unlist %>% as.vector

obs.n.call = dat.pam %>% 
  group_by(day) %>% 
  summarise(n = sum(num)) %>% 
  ungroup %>% select(n) %>% unlist %>% as.vector



# =============================================================================-
# loglik and CRPS ----
# =============================================================================-

## loglik ----

dat.loglik = data.frame()
for(i in 1:nrow(sim_design)){
  fiti = sim_design$Model[i]
  
  load(paste0(path.sum, datai, '_', fiti, '_sum.RData'))
  
  dat.loglik = rbind(
    dat.loglik,
    data.frame(
      Model = sim_design$ModelName[i],
      Aerial = post.loglik[,1],
      PAM = post.loglik[,2],
      Oblique = post.loglik[,3],
      Surface = post.loglik[,4]
    )
  ) 
}



## RPS ----

dat.rps = data.frame()
for(i in 1:nrow(sim_design)){
  fiti = sim_design$Model[i]
  
  load(paste0(path.sum, datai, '_', fiti, '_sum.RData'))
  
  for(t in 1:length(date.whale)){
    
    dummy = density(post.n.whale.dist[,t])
    obs.cdf = ifelse(dummy$x >= obs.n.whale.dist[t], 1, 0)
    post.cdf = cumsum(dummy$y) / sum(dummy$y)
    rps.aerial = sum((post.cdf - obs.cdf)^2) * diff(dummy$x)[1]
    
    dummy.pam = density(post.n.call[,t])
    obs.cdf.pam = ifelse(dummy.pam$x >= obs.n.call[t], 1, 0)
    post.cdf.pam = cumsum(dummy.pam$y) / sum(dummy.pam$y)
    rps.pam = sum((post.cdf.pam - obs.cdf.pam)^2) * diff(dummy.pam$x)[1]
    
    dat.rps = rbind(
      dat.rps,
      data.frame(
        Model = sim_design$ModelName[i],
        Date = date.whale[t],
        Aerial = rps.aerial,
        PAM = rps.pam
      )  
    )
  } 
}





# =============================================================================-
# table ----
# =============================================================================-


## full loglik ----
tab.loglik.full = dat.loglik %>% 
  mutate(Full = Aerial + PAM + Oblique + Surface) %>% 
  select(Model, Full) %>% 
  group_by(Model) %>% 
  summarise(
    Median = format(round(median(Full)), nsmall = 0),
    HPD = paste0('(', 
      format(round(HPDinterval(as.mcmc(Full))[1]), nsmall = 0), ', ', 
      format(round(HPDinterval(as.mcmc(Full))[2]), nsmall = 0), ')'
    )
  )
  

## RPS for whale abundance ----
tab.rps.aerial = dat.rps %>% 
  mutate(
    Aerial = format(round(Aerial, 2), nsmall = 2),
    PAM = format(round(PAM, 2), nsmall = 2)
  ) %>% 
  select(-PAM) %>% 
  pivot_wider(
    id_cols = Model,
    names_from = Date,
    values_from = Aerial
  )


## combine ----
cbind(
  tab.loglik.full %>% add_column(Empty = ''),
  tab.rps.aerial %>% ungroup %>% select(-Model)
) %>% xtable %>% 
  print(include.rownames = F)



