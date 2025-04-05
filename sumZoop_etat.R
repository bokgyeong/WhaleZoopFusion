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

fiti = 'Zoop_etat'

filename = paste0(path.fit, datai, '_', fiti, '.RData')
filename.sum = paste0(path.sum, datai, '_', fiti, '_sum.RData')

# =============================================================================-
# summarize results ----
# =============================================================================-

load(paste0(path.data, datai, '.RData'))
load(filename)


## observation days ----
Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix # 6 days of Zoop collection

sig2_obl = trpar$sig2_obl


## summarize ----
mean.logZ = mean.Z = 0 # posterior mean surfaces
post.n.zoop = c() # posterior samples for abundance
post.loglik = c() # posterior samples for log likelihood


start = 1
# load(filename3); start = nrow(post.n.whale) + 1

source(path.r)
sourceCpp(path.cpp)


for(iter in start:nrow(postAlpha)){
  logZ = cppf_logZ_a0t_etat_ja(Xzoop, postAlpha[iter,], sapply(1:length(postEta), function(t) postEta[[t]][iter,]))
  Z = exp(logZ)
  
  mean.logZ = ( mean.logZ * (iter - 1) + logZ ) / iter
  mean.Z = ( mean.Z * (iter - 1) + exp(logZ) ) / iter
  post.n.zoop = rbind(post.n.zoop, colMeans(exp(logZ)))
  
  post.loglik = rbind(
    post.loglik,
    c(cppf_loglik_obl(logYobl, cppf_logZobl(logZ, indobl-1), sig2_obl),
      cppf_loglik_sur(logYsur, cppf_logZsur(logZ, indsur-1), postLam[iter,1], postLam[iter,2], postSig2[iter,1])
    )
  )
  
  if((iter %% 1000) == 0){
    save(
      mean.logZ, mean.Z, post.n.zoop,
      post.loglik,
      file = filename.sum
    )
    cat(paste0('computed ', iter, 'th iteration...\n'))
  }
}
save(
  mean.logZ, mean.Z, post.n.zoop,
  post.loglik,
  file = filename.sum
)



