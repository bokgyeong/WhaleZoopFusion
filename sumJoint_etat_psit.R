rm(list = ls())
require(Rcpp); require(RcppArmadillo); require(foreach)
require(tidyverse); require(geodist)

runID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# runID = 1

# directory and file names ----
fold = 'github/'

path.data = paste0(fold, 'data/')
path.fit = paste0(fold, 'fit/')
path.sum = paste0(fold, 'sum/')
ifelse(!dir.exists(path.sum), dir.create(path.sum, recursive = T), FALSE)

path.r = paste0(fold, '/src/RFtns.R')
path.cpp = paste0(fold, '/src/RcppFtns.cpp')

data = c('sim_etat_psi', 'sim_etat_psit')
datai = data[runID]

# datai = 'sim_etat_psi'
# datai = 'sim_etat_psit'

fiti = 'Joint_etat_psit'

filename = paste0(path.fit, datai, '_', fiti, '.RData')
filename.sum = paste0(path.sum, datai, '_', fiti, '_sum.RData')

# =============================================================================-
# summarize results ----
# =============================================================================-

load(paste0(path.data, datai, '.RData'))
load(filename)


## observation days ----
Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix # 6 days of Zoop collection
Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix

indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2]))

sig2_obl = trpar$sig2_obl



## summarize ----
mean.logZ = mean.Z = mean.loglam = mean.lam = 0 # posterior mean surfaces
post.n.zoop = post.n.whale = c() # posterior samples for abundance
post.loglik = c() # posterior samples for log likelihood
post.n.whale.dist = post.n.call = c() # posterior samples for the number of observed whales and calls 

start = 1
# load(filename.sum); start = nrow(post.n.whale) + 1

source(path.r)
sourceCpp(path.cpp)

for(iter in start:nrow(postAlpha)){
  logZ = cppf_logZ_a0t_etat_ja(Xzoop, postAlpha[iter,], sapply(1:length(postEta), function(t) postEta[[t]][iter,]))
  Z = exp(logZ)
  
  loglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, postBeta[iter,], sapply(1:length(postPsi), function(t) postPsi[[t]][iter,]), indtwhale-1)
  lam = exp(loglam)
  
  mean.logZ = ( mean.logZ * (iter - 1) + logZ ) / iter
  mean.Z = ( mean.Z * (iter - 1) + Z ) / iter
  mean.loglam = ( mean.loglam * (iter - 1) + loglam ) / iter
  mean.lam = ( mean.lam * (iter - 1) + lam ) / iter
  post.n.zoop = rbind(post.n.zoop, colMeans(Z))
  post.n.whale = rbind(post.n.whale, colSums(lam) * cellarea)
  
  post.loglik = rbind(
    post.loglik,
    c(cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist-1, indmdist-1, indmwhale-1, cellarea),
      cppf_loglik_pam_ct(Wpam, lam, dur_dist, postC[iter,], p_pam, cellarea),
      cppf_loglik_obl(logYobl, cppf_logZobl(logZ, indobl-1), sig2_obl),
      cppf_loglik_sur(logYsur, cppf_logZsur(logZ, indsur-1), postLam[iter,1], postLam[iter,2], postSig2[iter,1])
    )
  )
  
  post.n.whale.dist = rbind(
    post.n.whale.dist,  
    t(cppf_intlam_dist(lam, pii, f_dist, indmwhale-1, cellarea))
  )
  
  post.n.call = rbind(
    post.n.call,
    unlist(lapply(cppf_lam_pam_ct(lam, dur_dist, postC[iter,], p_pam, cellarea), sum))
  )
  
  if((iter %% 1000) == 0){
    save(
      mean.logZ, mean.Z, post.n.zoop,
      mean.loglam, mean.lam, post.n.whale,
      post.loglik,
      post.n.whale.dist, post.n.call,
      file = filename.sum
    )
    cat(paste0('computed ', iter, 'th iteration...\n'))
  }
}

save(
  mean.logZ, mean.Z, post.n.zoop,
  mean.loglam, mean.lam, post.n.whale,
  post.loglik,
  post.n.whale.dist, post.n.call,
  file = filename.sum
)
