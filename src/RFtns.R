# Simulate point patters over irregular spatial domain
Rf_simulateNHPP = function(lam, coord, x.range, y.range, mas){
  res = list()
  
  for(t in 1:ncol(lam)){
    
    maxintensity = max(lam[,t])
    dummy.x = data.frame(lon = x.range, lat = y.range[1])
    dummy.y = data.frame(lon = x.range[1], lat = y.range) 
    nn = rpois(1, maxintensity * 
                 geodist(dummy.x[1,], dummy.x[2,]) / 1000 * 
                 geodist(dummy.y[1,], dummy.y[2,]) / 1000)
    
    points = cbind(runif(nn, x.range[1], x.range[2]), 
                   runif(nn, y.range[1], y.range[2]))
    
    land = which(point.in.polygon(points[,1], points[,2], mas[,1], mas[,2], mode.checked = FALSE) == 1) # find points on land
    land = c(land, which(points[,1] > -70.08 & points[,2] > 41.95))
    xy = points[-land,] # remove points on land
    xy = xy[-c(which(xy[,1] < -70.6 & xy[,2] < 41.8),  # remove points on bottom left and top right
               which(xy[,1]> -70.15 & xy[,2] > 42.06), 
               which(xy[,1]> -70.13 & xy[,2] > 42.05),
               which(xy[,2] < 41.735), # remove points on where whales don't go
               which(xy[,1] < -70.635 & xy[,2] < 42.05),
               which(xy[,1] < -70.605 & xy[,2] < 42.01),
               which(xy[,1] < -70.63 & xy[,2] > 42.0 & xy[,2] < 42.04),
               which(xy[,1] < -70.61 & xy[,2] > 42.0 & xy[,2] < 42.02)),]
    
    xy = as.data.frame(xy)
    colnames(xy) = c('lon', 'lat')
    n.xy = nrow(xy)
    
    res[[t]] = foreach(i = 1:n.xy, .combine = rbind) %do% {
      inds = which.min( geodist(coord, xy[i,]) )
      if( runif(1) < (lam[inds,t] / maxintensity) ){
        xy[i,]
      }  
    }
  }
  
  return( res )
}


# acoustic detection function 
Rf_p_pam = function(noise, si, sj){
  thres = 26 + noise + 14.5 * log10( geodist(si, sj)[1,1] ) # unit is meters
  
  if(thres < 141) {
    return(1)
  } else if(thres <= 197){
    return(1 - (thres - 141) / (197 - 141))
  } else {
    return(0)
  }
  
  # 1-pnorm(noise + 26 + 14.5 * log10( geodist(si, sj)[1,1] ), mean = 192, sd = sqrt(10))
}


Rf_fitJoint_etat_psi = function(
    dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
    dur_dist, f_dist_obs, f_dist, p_pam,
    coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, shape_c, scale_c,
    n.iter, n.keep, n.save, n.adap, new.run,
    startval, path.cpp, filename){

  # ===========================================================================-
  # Zoop ----
  # ===========================================================================-
  nG = nrow(coord)
  distmat = geodist(coord) / 1000
  Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix


  # ---------------------------------------------------------------------------=
  ## oblique tow ----
  # ---------------------------------------------------------------------------=
  dat.obl = dat.zoop %>%
    filter(type == 'obl') %>%
    dplyr::select(month, day, type, lon, lat, zoop)

  indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
    indt = which(dat.obl$day[i] == Tz[,2])
    inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
    c(inds, indt)
  }
  rownames(indobl) = NULL
  logYobl = log(dat.obl$zoop)

  # ---------------------------------------------------------------------------=
  ## surface tow ----
  # ---------------------------------------------------------------------------=
  dat.sur = dat.zoop %>%
    filter(type == 'sfc') %>%
    dplyr::select(month, day, type, lon, lat, zoop)

  indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
    indt = which(dat.sur$day[i] == Tz[,2])
    inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
    c(inds, indt)
  }
  rownames(indsur) = NULL
  logYsur = log(dat.sur$zoop)


  # ===========================================================================-
  # whale ----
  # ===========================================================================-

  Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
  indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
  indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data


  # ---------------------------------------------------------------------------=
  ## Areal sightings ----
  # ---------------------------------------------------------------------------=
  Sdist = cbind(dat.dist$lon, dat.dist$lat)

  inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
    indt = which(dat.dist$day[i] == Tw[,2])
    inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
    c(inds, indt)
  }
  rownames(inddist) = NULL

  indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data


  # ---------------------------------------------------------------------------=
  ## PAM ----
  # ---------------------------------------------------------------------------=
  Wpam = foreach(i = 1:nrow(Tw)) %do% {
    dat.pam %>%
      filter(day == Tw[i,2]) %>%
      dplyr::select(num) %>%
      unlist %>% as.vector
  }


  # ===========================================================================-
  # MCMC ----
  # ===========================================================================-

  outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))

  updateCOV = TRUE
  adaptFactorExponent = 0.8
  adaptInterval = n.adap

  if(new.run){

    # initial parameter values for Zoop model
    alpha0tilde = startval$alpha0tilde
    alpha = startval$alpha
    kappa_eta = startval$kappa_eta
    eta = startval$eta
    lam0 = startval$lam0
    lam1 = startval$lam1
    sig2_sur = startval$sig2_sur

    # initial parameter values for whale model
    beta = c(startval$beta0, startval$beta)
    kappa_psi = startval$kappa_psi
    psi = startval$psi
    c = startval$c

    # MH proposal covariance matrices
    sigma2 = rep(0.05, 1+1)
    COValpha = diag(length(alpha))
    COVbeta = diag(length(beta))
    adapIter = rep(1, length(sigma2))

    start = 1; Accprob = 0; rtime = 0
    postAlpha0tilde = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = postC = c()
    postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)

  } else {
    load(filename)
    start = which(outers == curniter)
  }

  sourceCpp(path.cpp)

  for(i in (start+1):length(outers) ){
    outeri = outers[i]-outers[i-1]
    curniter = outers[i]

    ptm = proc.time()[3]
    dummy = cppf_fitJoint_etat_psi(
      outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
      Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
      indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha,
      eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur,
      beta, psi, pii, c, shape_c, scale_c, kappa_psi, phi_psi,
      cellarea, sigma2, COValpha, COVbeta, updateCOV, adaptInterval,
      adaptFactorExponent, adapIter)
    rtime = rtime + proc.time()[3] - ptm

    Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]

    postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
    postAlpha = rbind(postAlpha, dummy$postAlpha)
    for(t in 1:ncol(eta)){
      postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
    }
    postKappa = rbind(postKappa, dummy$postKappa)
    postLam = rbind(postLam, dummy$postLam)
    postSig2 = rbind(postSig2, dummy$postSig2)
    postBeta = rbind(postBeta, dummy$postBeta)
    postPsi = rbind(postPsi, dummy$postPsi)
    postC = rbind(postC, dummy$postC)

    nSamples = nrow(postAlpha)
    alpha0tilde = postAlpha0tilde[nSamples,]
    alpha = postAlpha[nSamples,]
    eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
    lam0 = postLam[nSamples,1]
    lam1 = postLam[nSamples,2]
    sig2_sur = postSig2[nSamples,1]
    kappa_eta = postKappa[nSamples,1]
    beta = postBeta[nSamples,]
    psi = postPsi[nSamples,]
    kappa_psi = postKappa[nSamples,2]
    c = postC[nSamples,]

    sigma2 = dummy$sigma2
    adapIter = dummy$adapIter
    COValpha = dummy$COValpha
    COVbeta = dummy$COVbeta


    if(nSamples == (n.keep + n.save)){
      postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
      postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
      for(t in 1:ncol(eta)){
        postEta[[t]] = postEta[[t]][-(1:n.save),]
      }
      postKappa = postKappa[-(1:n.save),]
      postLam = postLam[-(1:n.save),]
      postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
      postBeta = postBeta[-(1:n.save),]
      postPsi = postPsi[-(1:n.save),]
      postC = matrix(postC[-(1:n.save),], nrow = n.keep)
    }


    save(Accprob, rtime, curniter,
         postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2,
         postBeta, postPsi, postC,
         distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
         Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
         inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
         kappa_eta, phi_eta, beta, psi, shape_c, scale_c,
         pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta,
         adapIter, file = filename)
  }
}


Rf_fitJoint_etat_psit = function(
    dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
    dur_dist, f_dist_obs, f_dist, p_pam,
    coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, shape_c, scale_c,
    n.iter, n.keep, n.save, n.adap, new.run,
    startval, path.cpp, filename){

  # ===========================================================================-
  # Zoop ----
  # ===========================================================================-
  nG = nrow(coord)
  distmat = geodist(coord) / 1000
  Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix


  # ---------------------------------------------------------------------------=
  ## oblique tow ----
  # ---------------------------------------------------------------------------=
  dat.obl = dat.zoop %>%
    filter(type == 'obl') %>%
    dplyr::select(month, day, type, lon, lat, zoop)

  indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
    indt = which(dat.obl$day[i] == Tz[,2])
    inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
    c(inds, indt)
  }
  rownames(indobl) = NULL
  logYobl = log(dat.obl$zoop)

  # ---------------------------------------------------------------------------=
  ## surface tow ----
  # ---------------------------------------------------------------------------=
  dat.sur = dat.zoop %>%
    filter(type == 'sfc') %>%
    dplyr::select(month, day, type, lon, lat, zoop)

  indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
    indt = which(dat.sur$day[i] == Tz[,2])
    inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
    c(inds, indt)
  }
  rownames(indsur) = NULL
  logYsur = log(dat.sur$zoop)


  # ===========================================================================-
  # whale ----
  # ===========================================================================-

  Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
  indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
  indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data


  # ---------------------------------------------------------------------------=
  ## Areal sightings ----
  # ---------------------------------------------------------------------------=
  Sdist = cbind(dat.dist$lon, dat.dist$lat)

  inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
    indt = which(dat.dist$day[i] == Tw[,2])
    inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
    c(inds, indt)
  }
  rownames(inddist) = NULL

  indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data


  # ---------------------------------------------------------------------------=
  ## PAM ----
  # ---------------------------------------------------------------------------=
  Wpam = foreach(i = 1:nrow(Tw)) %do% {
    dat.pam %>%
      filter(day == Tw[i,2]) %>%
      dplyr::select(num) %>%
      unlist %>% as.vector
  }


  # ===========================================================================-
  # MCMC ----
  # ===========================================================================-

  outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))

  updateCOV = TRUE
  adaptFactorExponent = 0.8
  adaptInterval = n.adap

  if(new.run){

    # initial parameter values for Zoop model
    alpha0tilde = startval$alpha0tilde
    alpha = startval$alpha
    kappa_eta = startval$kappa_eta
    eta = startval$eta
    lam0 = startval$lam0
    lam1 = startval$lam1
    sig2_sur = startval$sig2_sur

    # initial parameter values for whale model
    beta = c(startval$beta0, startval$beta)
    kappa_psi = startval$kappa_psi
    psi = startval$psi
    c = startval$c

    # MH proposal covariance matrices
    sigma2 = rep(0.05, 1+1)
    COValpha = diag(length(alpha))
    COVbeta = diag(length(beta))
    adapIter = rep(1, length(sigma2))

    start = 1; Accprob = 0; rtime = 0
    postAlpha0tilde = postAlpha = postKappa = postLam = postSig2 = postBeta = postC = c()
    postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
    postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)

  } else {
    load(filename)
    start = which(outers == curniter)
  }

  sourceCpp(path.cpp)

  for(i in (start+1):length(outers) ){
    outeri = outers[i]-outers[i-1]
    curniter = outers[i]

    ptm = proc.time()[3]
    dummy = cppf_fitJoint_etat_psit(
      outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
      Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
      indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha,
      eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur,
      beta, psi, pii, c, shape_c, scale_c, kappa_psi, phi_psi,
      cellarea, sigma2, COValpha, COVbeta, updateCOV, adaptInterval,
      adaptFactorExponent, adapIter)
    rtime = rtime + proc.time()[3] - ptm

    Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]

    postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
    postAlpha = rbind(postAlpha, dummy$postAlpha)
    for(t in 1:ncol(eta)){
      postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
    }
    postKappa = rbind(postKappa, dummy$postKappa)
    postLam = rbind(postLam, dummy$postLam)
    postSig2 = rbind(postSig2, dummy$postSig2)
    postBeta = rbind(postBeta, dummy$postBeta)
    for(t in 1:ncol(psi)){
      postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
    }
    postC = rbind(postC, dummy$postC)

    nSamples = nrow(postAlpha)
    alpha0tilde = postAlpha0tilde[nSamples,]
    alpha = postAlpha[nSamples,]
    eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
    lam0 = postLam[nSamples,1]
    lam1 = postLam[nSamples,2]
    sig2_sur = postSig2[nSamples,1]
    kappa_eta = postKappa[nSamples,1]
    beta = postBeta[nSamples,]
    psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
    kappa_psi = postKappa[nSamples,2]
    c = postC[nSamples,]

    sigma2 = dummy$sigma2
    adapIter = dummy$adapIter
    COValpha = dummy$COValpha
    COVbeta = dummy$COVbeta


    if(nSamples == (n.keep + n.save)){
      postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
      postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
      for(t in 1:ncol(eta)){
        postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
      }
      postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
      postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
      postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
      postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
      for(t in 1:ncol(psi)){
        postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
      }
      postC = matrix(postC[-(1:n.save),], nrow = n.keep)
    }


    save(Accprob, rtime, curniter,
         postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2,
         postBeta, postPsi, postC,
         distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
         Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
         inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
         kappa_eta, phi_eta, beta, psi, shape_c, scale_c,
         pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta,
         adapIter, file = filename)
  }
}


Rf_fitWhale_psit = function(
    dat.dist, dat.pam, Xwhale,
    dur_dist, f_dist_obs, f_dist, p_pam,
    coord, cellarea, pii, phi_psi, shape_c, scale_c,
    n.iter, n.keep, n.save, n.adap, new.run,
    startval, path.cpp, filename){

  # ===========================================================================-
  # whale ----
  # ===========================================================================-
  nG = nrow(coord)
  distmat = geodist(coord) / 1000

  Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
  indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data


  # ---------------------------------------------------------------------------=
  ## Areal sightings ----
  # ---------------------------------------------------------------------------=
  Sdist = cbind(dat.dist$lon, dat.dist$lat)

  inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
    indt = which(dat.dist$day[i] == Tw[,2])
    inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
    c(inds, indt)
  }
  rownames(inddist) = NULL

  indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data


  # ---------------------------------------------------------------------------=
  ## PAM ----
  # ---------------------------------------------------------------------------=
  Wpam = foreach(i = 1:nrow(Tw)) %do% {
    dat.pam %>%
      filter(day == Tw[i,2]) %>%
      dplyr::select(num) %>%
      unlist %>% as.vector
  }


  # ===========================================================================-
  # MCMC ----
  # ===========================================================================-

  outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))

  updateCOV = TRUE
  adaptFactorExponent = 0.8
  adaptInterval = n.adap

  if(new.run){

    # initial parameter values for whale model
    beta = c(startval$beta0, startval$beta)
    kappa_psi = startval$kappa_psi
    psi = startval$psi
    c = startval$c

    # MH proposal covariance matrices
    sigma2 = rep(0.05, 1)
    COVbeta = diag(length(beta))
    adapIter = rep(1, length(sigma2))

    start = 1; Accprob = 0; rtime = 0
    postKappa = postBeta = postC = c()
    postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)

  } else {
    load(filename)
    start = which(outers == curniter)
  }

  sourceCpp(path.cpp)

  for(i in (start+1):length(outers) ){
    outeri = outers[i]-outers[i-1]
    curniter = outers[i]

    ptm = proc.time()[3]
    dummy = cppf_fitWhale_psit(
      outeri, distmat, Xwhale, Sdist, Wpam, dur_dist, f_dist_obs,
      f_dist, p_pam, inddist, indmdist, indmwhale, beta, psi, pii,
      c, shape_c, scale_c, kappa_psi, phi_psi, cellarea, sigma2,
      COVbeta, updateCOV, adaptInterval, adaptFactorExponent, adapIter)
    rtime = rtime + proc.time()[3] - ptm

    Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]

    postKappa = rbind(postKappa, dummy$postKappa)
    postBeta = rbind(postBeta, dummy$postBeta)
    for(t in 1:ncol(psi)){
      postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
    }
    postC = rbind(postC, dummy$postC)

    nSamples = nrow(postBeta)
    beta = postBeta[nSamples,]
    psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
    kappa_psi = postKappa[nSamples,1]
    c = postC[nSamples,]

    sigma2 = dummy$sigma2
    adapIter = dummy$adapIter
    COVbeta = dummy$COVbeta


    if(nSamples == (n.keep + n.save)){
      postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
      postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
      for(t in 1:ncol(psi)){
        postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
      }
      postC = matrix(postC[-(1:n.save),], nrow = n.keep)
    }


    save(Accprob, rtime, curniter,
         postKappa, postBeta, postPsi, postC,
         distmat, Xwhale, Sdist,
         Wpam, dur_dist, f_dist_obs, f_dist, p_pam,
         inddist, indmdist, indmwhale, beta, psi, shape_c, scale_c,
         pii, c, kappa_psi, phi_psi, cellarea, sigma2, COVbeta,
         adapIter, file = filename)
  }
}

# Rf_fitWhale_b0t_psi_ct = function(
#     dat.dist, dat.pam, Xwhale, 
#     dur_dist, f_dist_obs, f_dist, p_pam,
#     coord, cellarea, pii, phi_psi, shape_c, scale_c,
#     n.iter, n.keep, n.save, n.adap, new.run, 
#     startval, path.cpp, filename){
#   
#   # ===========================================================================-
#   # whale ----
#   # ===========================================================================-
#   nG = nrow(coord)
#   distmat = geodist(coord) / 1000
#   
#   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
#   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
#   
#   
#   # ---------------------------------------------------------------------------=
#   ## Areal sightings ----
#   # ---------------------------------------------------------------------------=
#   Sdist = cbind(dat.dist$lon, dat.dist$lat)
#   
#   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
#     indt = which(dat.dist$day[i] == Tw[,2])
#     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
#     c(inds, indt)
#   }
#   rownames(inddist) = NULL
#   
#   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
#   
#   
#   # ---------------------------------------------------------------------------=
#   ## PAM ----
#   # ---------------------------------------------------------------------------=
#   Wpam = foreach(i = 1:nrow(Tw)) %do% {
#     dat.pam %>% 
#       filter(day == Tw[i,2]) %>% 
#       dplyr::select(num) %>% 
#       unlist %>% as.vector
#   }
#   
#   
#   # ===========================================================================-
#   # MCMC ----
#   # ===========================================================================-
#   
#   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
#   
#   updateCOV = TRUE
#   adaptFactorExponent = 0.8
#   adaptInterval = n.adap
#   
#   if(new.run){
#     
#     # initial parameter values for whale model
#     beta = c(startval$beta0, startval$beta)
#     kappa_psi = startval$kappa_psi
#     psi = startval$psi
#     c = startval$c
#     
#     # MH proposal covariance matrices
#     sigma2 = rep(0.05, 1)
#     COVbeta = diag(length(beta))
#     adapIter = rep(1, length(sigma2))
#     
#     start = 1; Accprob = 0; rtime = 0
#     postKappa = postBeta = postC = postPsi = c()
#     
#   } else {
#     load(filename)
#     start = which(outers == curniter)
#   }
#   
#   sourceCpp(path.cpp)
#   
#   for(i in (start+1):length(outers) ){
#     outeri = outers[i]-outers[i-1]
#     curniter = outers[i]
#     
#     ptm = proc.time()[3]
#     dummy = cppf_fitWhale_b0t_psi_ct(
#       outeri, distmat, Xwhale, Sdist, Wpam, dur_dist, f_dist_obs, 
#       f_dist, p_pam, inddist, indmdist, indmwhale, beta, psi, pii, 
#       c, shape_c, scale_c, kappa_psi, phi_psi, cellarea, sigma2, 
#       COVbeta, updateCOV, adaptInterval, adaptFactorExponent, adapIter) 
#     rtime = rtime + proc.time()[3] - ptm
#     
#     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
#     
#     postKappa = rbind(postKappa, dummy$postKappa)
#     postBeta = rbind(postBeta, dummy$postBeta)
#     postPsi = rbind(postPsi, dummy$postPsi)
#     postC = rbind(postC, dummy$postC)
#     
#     nSamples = nrow(postBeta)
#     beta = postBeta[nSamples,]
#     psi = postPsi[nSamples,]
#     kappa_psi = postKappa[nSamples,1]
#     c = postC[nSamples,]
#     
#     sigma2 = dummy$sigma2
#     adapIter = dummy$adapIter
#     COVbeta = dummy$COVbeta
#     
#     
#     if(nSamples == (n.keep + n.save)){
#       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
#       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
#       postPsi = matrix(postPsi[-(1:n.save),], nrow = n.keep)
#       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
#     }
#     
#     
#     save(Accprob, rtime, curniter,
#          postKappa, postBeta, postPsi, postC,
#          distmat, Xwhale, Sdist, 
#          Wpam, dur_dist, f_dist_obs, f_dist, p_pam,
#          inddist, indmdist, indmwhale, beta, psi, shape_c, scale_c,
#          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COVbeta, 
#          adapIter, file = filename)
#   }
# }


Rf_fitZoop_etat = function(
    dat.zoop, Xzoop,
    coord, cellarea, phi_eta, sig2_obl, tau2,
    n.iter, n.keep, n.save, n.adap, new.run,
    startval, path.cpp, filename){

  # ===========================================================================-
  # Zoop ----
  # ===========================================================================-
  nG = nrow(coord)
  distmat = geodist(coord) / 1000
  Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix


  # ---------------------------------------------------------------------------=
  ## oblique tow ----
  # ---------------------------------------------------------------------------=
  dat.obl = dat.zoop %>%
    filter(type == 'obl') %>%
    dplyr::select(month, day, type, lon, lat, zoop)

  indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
    indt = which(dat.obl$day[i] == Tz[,2])
    inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
    c(inds, indt)
  }
  rownames(indobl) = NULL
  logYobl = log(dat.obl$zoop)

  # ---------------------------------------------------------------------------=
  ## surface tow ----
  # ---------------------------------------------------------------------------=
  dat.sur = dat.zoop %>%
    filter(type == 'sfc') %>%
    dplyr::select(month, day, type, lon, lat, zoop)

  indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
    indt = which(dat.sur$day[i] == Tz[,2])
    inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
    c(inds, indt)
  }
  rownames(indsur) = NULL
  logYsur = log(dat.sur$zoop)


  # ===========================================================================-
  # MCMC ----
  # ===========================================================================-

  outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))

  updateCOV = TRUE
  adaptFactorExponent = 0.8
  adaptInterval = n.adap

  if(new.run){

    # initial parameter values for Zoop model
    alpha0tilde = startval$alpha0tilde
    alpha = startval$alpha
    kappa_eta = startval$kappa_eta
    eta = startval$eta
    lam0 = startval$lam0
    lam1 = startval$lam1
    sig2_sur = startval$sig2_sur

    # MH proposal covariance matrices
    sigma2 = rep(0.05, 1)
    COValpha = diag(length(alpha))
    adapIter = rep(1, length(sigma2))

    start = 1; Accprob = 0; rtime = 0
    postAlpha0tilde = postAlpha = postKappa = postLam = postSig2 = c()
    postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)

  } else {
    load(filename)
    start = which(outers == curniter)
  }

  sourceCpp(path.cpp)

  for(i in (start+1):length(outers) ){
    outeri = outers[i]-outers[i-1]
    curniter = outers[i]

    ptm = proc.time()[3]
    dummy = cppf_fitZoop_etat(
      outeri, distmat, Xzoop, logYobl, logYsur, indobl, indsur,
      alpha0tilde, alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta,
      sig2_obl, sig2_sur, cellarea, sigma2, COValpha, updateCOV,
      adaptInterval, adaptFactorExponent, adapIter)
    rtime = rtime + proc.time()[3] - ptm

    Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]

    postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
    postAlpha = rbind(postAlpha, dummy$postAlpha)
    for(t in 1:ncol(eta)){
      postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
    }
    postKappa = rbind(postKappa, dummy$postKappa)
    postLam = rbind(postLam, dummy$postLam)
    postSig2 = rbind(postSig2, dummy$postSig2)


    nSamples = nrow(postAlpha)
    alpha0tilde = postAlpha0tilde[nSamples,]
    alpha = postAlpha[nSamples,]
    eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
    lam0 = postLam[nSamples,1]
    lam1 = postLam[nSamples,2]
    sig2_sur = postSig2[nSamples,1]
    kappa_eta = postKappa[nSamples,1]

    sigma2 = dummy$sigma2
    adapIter = dummy$adapIter
    COValpha = dummy$COValpha


    if(nSamples == (n.keep + n.save)){
      postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
      postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
      for(t in 1:ncol(eta)){
        postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
      }
      postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
      postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
      postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
    }


    save(Accprob, rtime, curniter,
         postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2,
         distmat, Xzoop, logYobl, logYsur, indobl, indsur,
         alpha0tilde, alpha, eta, lam0, lam1, sig2_sur, kappa_eta, phi_eta,
         cellarea, sigma2, COValpha,
         adapIter, file = filename)
  }
}



# # =============================================================================-
# # uniform prior for c ----
# # =============================================================================-
# 
# Rf_fitJoint_cZ_a0t_etat_b0t_psi_ct_ja_uc = function(
#     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
#     dur_dist, f_dist_obs, f_dist, p_pam,
#     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, lb_c,
#     n.iter, n.keep, n.save, n.adap, new.run, 
#     startval, path.cpp, filename){
#   
#   # ===========================================================================-
#   # Zoop ----
#   # ===========================================================================-
#   nG = nrow(coord)
#   distmat = geodist(coord) / 1000
#   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
#   
#   
#   # ---------------------------------------------------------------------------=
#   ## oblique tow ----
#   # ---------------------------------------------------------------------------=
#   dat.obl = dat.zoop %>% 
#     filter(type == 'obl') %>% 
#     dplyr::select(month, day, type, lon, lat, zoop)
#   
#   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
#     indt = which(dat.obl$day[i] == Tz[,2])
#     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
#     c(inds, indt)
#   }
#   rownames(indobl) = NULL
#   logYobl = log(dat.obl$zoop)
#   
#   # ---------------------------------------------------------------------------=
#   ## surface tow ----
#   # ---------------------------------------------------------------------------=
#   dat.sur = dat.zoop %>% 
#     filter(type == 'sfc') %>% 
#     dplyr::select(month, day, type, lon, lat, zoop)
#   
#   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
#     indt = which(dat.sur$day[i] == Tz[,2])
#     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
#     c(inds, indt)
#   }
#   rownames(indsur) = NULL
#   logYsur = log(dat.sur$zoop)
#   
#   
#   # ===========================================================================-
#   # whale ----
#   # ===========================================================================-
#   
#   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
#   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
#   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
#   
#   
#   # ---------------------------------------------------------------------------=
#   ## Areal sightings ----
#   # ---------------------------------------------------------------------------=
#   Sdist = cbind(dat.dist$lon, dat.dist$lat)
#   
#   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
#     indt = which(dat.dist$day[i] == Tw[,2])
#     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
#     c(inds, indt)
#   }
#   rownames(inddist) = NULL
#   
#   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
#   
#   
#   # ---------------------------------------------------------------------------=
#   ## PAM ----
#   # ---------------------------------------------------------------------------=
#   Wpam = foreach(i = 1:nrow(Tw)) %do% {
#     dat.pam %>% 
#       filter(day == Tw[i,2]) %>% 
#       dplyr::select(num) %>% 
#       unlist %>% as.vector
#   }
#   
#   
#   # ===========================================================================-
#   # MCMC ----
#   # ===========================================================================-
#   
#   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
#   
#   updateCOV = TRUE
#   adaptFactorExponent = 0.8
#   adaptInterval = n.adap
#   
#   if(new.run){
#     
#     # initial parameter values for Zoop model
#     alpha0tilde = startval$alpha0tilde
#     alpha = startval$alpha
#     kappa_eta = startval$kappa_eta
#     eta = startval$eta
#     lam0 = startval$lam0
#     lam1 = startval$lam1
#     sig2_sur = startval$sig2_sur
#     
#     # initial parameter values for whale model
#     beta = c(startval$beta0, startval$beta)
#     kappa_psi = startval$kappa_psi
#     psi = startval$psi
#     c = startval$c
#     
#     # MH proposal covariance matrices
#     sigma2 = rep(0.05, 1+1+1)
#     COValpha = diag(length(alpha))
#     COVbeta = diag(length(beta))
#     COVc = diag(length(c))
#     adapIter = rep(1, length(sigma2))
#     
#     start = 1; Accprob = 0; rtime = 0
#     postAlpha0tilde = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = postC = c()
#     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)  
#     
#   } else {
#     load(filename)
#     start = which(outers == curniter)
#   }
#   
#   sourceCpp(path.cpp)
#   
#   for(i in (start+1):length(outers) ){
#     outeri = outers[i]-outers[i-1]
#     curniter = outers[i]
#     
#     ptm = proc.time()[3]
#     dummy = cppf_fitJoint_cZ_a0t_etat_b0t_psi_ct_ja_uc(
#       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
#       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
#       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha, 
#       eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur, 
#       beta, psi, pii, c, lb_c, kappa_psi, phi_psi, cellarea, 
#       sigma2, COValpha, COVbeta, COVc, updateCOV, adaptInterval, 
#       adaptFactorExponent, adapIter) 
#     rtime = rtime + proc.time()[3] - ptm
#     
#     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
#     
#     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
#     postAlpha = rbind(postAlpha, dummy$postAlpha)
#     for(t in 1:ncol(eta)){
#       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
#     }
#     postKappa = rbind(postKappa, dummy$postKappa)
#     postLam = rbind(postLam, dummy$postLam)
#     postSig2 = rbind(postSig2, dummy$postSig2)
#     postBeta = rbind(postBeta, dummy$postBeta)
#     postPsi = rbind(postPsi, dummy$postPsi)
#     postC = rbind(postC, dummy$postC)
#     
#     nSamples = nrow(postAlpha)
#     alpha0tilde = postAlpha0tilde[nSamples,]
#     alpha = postAlpha[nSamples,]
#     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
#     lam0 = postLam[nSamples,1]
#     lam1 = postLam[nSamples,2]
#     sig2_sur = postSig2[nSamples,1]
#     kappa_eta = postKappa[nSamples,1]
#     beta = postBeta[nSamples,]
#     psi = postPsi[nSamples,]
#     kappa_psi = postKappa[nSamples,2]
#     c = postC[nSamples,]
#     
#     sigma2 = dummy$sigma2
#     adapIter = dummy$adapIter
#     COValpha = dummy$COValpha
#     COVbeta = dummy$COVbeta
#     COVc = dummy$COVc
#     
#     
#     if(nSamples == (n.keep + n.save)){
#       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
#       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
#       for(t in 1:ncol(eta)){
#         postEta[[t]] = postEta[[t]][-(1:n.save),]
#       }
#       postKappa = postKappa[-(1:n.save),]
#       postLam = postLam[-(1:n.save),]
#       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
#       postBeta = postBeta[-(1:n.save),]
#       postPsi = postPsi[-(1:n.save),]
#       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
#     }
#     
#     
#     save(Accprob, rtime, curniter,
#          postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2, 
#          postBeta, postPsi, postC,
#          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
#          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
#          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
#          kappa_eta, phi_eta, beta, psi, lb_c, 
#          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, COVc,
#          adapIter, file = filename)
#   }
# }
# 
# 
# # =============================================================================-
# # ****** old ***** ----
# # =============================================================================-
# 
# 
# # # =============================================================================-
# # # fix call rates ----
# # # =============================================================================-
# # ## not centering logZ ----
# # # -----------------------------------------------------------------------------=
# # Rf_fitJoint_eta_b0t_psi_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha = postEta = postKappa = postLam = postSig2 = postBeta = postPsi = c()
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_eta_b0t_psi_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha, eta, lam0, 
# #       lam1, kappa_eta, phi_eta, sig2_obl, sig2_sur, beta, psi, 
# #       pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #       updateCOV, adaptInterval, adaptFactorExponent, adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     postEta = rbind(postEta, dummy$postEta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha = postAlpha[nSamples,]
# #     eta = postEta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       postEta = matrix(postEta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postPsi = matrix(postPsi[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, 
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # Rf_fitJoint_eta_b0t_psit_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>% 
# #     arrange(MARU) %>% 
# #     dplyr::select(lon, lat) %>% 
# #     unique() %>% 
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha = postEta = postKappa = postLam = postSig2 = postBeta = c()
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_eta_b0t_psit_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha, eta, lam0, 
# #       lam1, kappa_eta, phi_eta, sig2_obl, sig2_sur, beta, psi, 
# #       pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #       updateCOV, adaptInterval, adaptFactorExponent, adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     postEta = rbind(postEta, dummy$postEta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha = postAlpha[nSamples,]
# #     eta = postEta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       postEta = matrix(postEta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, 
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # Rf_fitJoint_a0t_etat_b0t_psi_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_a0t_etat_b0t_psi_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha, 
# #       eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur, 
# #       beta, psi, pii, c, kappa_psi, phi_psi, cellarea, sigma2, 
# #       COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent, 
# #       adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = postEta[[t]][-(1:n.save),]
# #       }
# #       postKappa = postKappa[-(1:n.save),]
# #       postLam = postLam[-(1:n.save),]
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = postBeta[-(1:n.save),]
# #       postPsi = postPsi[-(1:n.save),]
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, 
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # Rf_fitJoint_a0t_etat_b0t_psit_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postKappa = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)  
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_a0t_etat_b0t_psit_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha, 
# #       eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur, 
# #       beta, psi, pii, c, kappa_psi, phi_psi, cellarea, sigma2, 
# #       COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent, 
# #       adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, 
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # Rf_fitJoint_a0t_eta_b0t_psit_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postEta = postKappa = postLam = postSig2 = postBeta = c()
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_a0t_eta_b0t_psit_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha, 
# #       eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur, 
# #       beta, psi, pii, c, kappa_psi, phi_psi, cellarea, sigma2, 
# #       COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent, 
# #       adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     postEta = rbind(postEta, dummy$postEta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = postEta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       postEta = matrix(postEta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, 
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # # -----------------------------------------------------------------------------=
# # ## centering logZ ----
# # # -----------------------------------------------------------------------------=
# # 
# # Rf_fitJoint_cZ_eta_b0t_psi_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha = postEta = postKappa = postLam = postSig2 = postBeta = postPsi = c()
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_cZ_eta_b0t_psi_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha, eta, lam0, 
# #       lam1, kappa_eta, phi_eta, sig2_obl, sig2_sur, beta, psi, 
# #       pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #       updateCOV, adaptInterval, adaptFactorExponent, adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     postEta = rbind(postEta, dummy$postEta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha = postAlpha[nSamples,]
# #     eta = postEta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       postEta = matrix(postEta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postPsi = matrix(postPsi[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, 
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # Rf_fitJoint_cZ_eta_b0t_psit_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha = postEta = postKappa = postLam = postSig2 = postBeta = c()
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_cZ_eta_b0t_psit_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha, eta, lam0, 
# #       lam1, kappa_eta, phi_eta, sig2_obl, sig2_sur, beta, psi, 
# #       pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #       updateCOV, adaptInterval, adaptFactorExponent, adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     postEta = rbind(postEta, dummy$postEta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha = postAlpha[nSamples,]
# #     eta = postEta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       postEta = matrix(postEta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, 
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # Rf_fitJoint_cZ_a0t_etat_b0t_psi_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_cZ_a0t_etat_b0t_psi_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha, 
# #       eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur, 
# #       beta, psi, pii, c, kappa_psi, phi_psi, cellarea, sigma2, 
# #       COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent, 
# #       adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = postEta[[t]][-(1:n.save),]
# #       }
# #       postKappa = postKappa[-(1:n.save),]
# #       postLam = postLam[-(1:n.save),]
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = postBeta[-(1:n.save),]
# #       postPsi = postPsi[-(1:n.save),]
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, 
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # # =============================================================================-
# # # estimate call rates ----
# # # =============================================================================-
# # 
# # 
# # # -----------------------------------------------------------------------------=
# # ## centering logZ ----
# # # -----------------------------------------------------------------------------=
# # 
# # Rf_fitJoint_cZ_eta_b0t_psi = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha = postEta = postKappa = postLam = postSig2 = postBeta = postPsi = postC = c()
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_cZ_eta_b0t_psi(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha, eta, lam0, 
# #       lam1, kappa_eta, phi_eta, sig2_obl, sig2_sur, beta, psi, 
# #       pii, c, shape_c, scale_c, kappa_psi, phi_psi, cellarea, sigma2, 
# #       COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent, 
# #       adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     postEta = rbind(postEta, dummy$postEta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha = postAlpha[nSamples,]
# #     eta = postEta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       postEta = matrix(postEta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postPsi = matrix(postPsi[-(1:n.save),], nrow = n.keep)
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, shape_c, scale_c,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # Rf_fitJoint_cZ_eta_b0t_psit = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>% 
# #     arrange(MARU) %>% 
# #     dplyr::select(lon, lat) %>% 
# #     unique() %>% 
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha = postEta = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_cZ_eta_b0t_psit(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha, eta, lam0, 
# #       lam1, kappa_eta, phi_eta, sig2_obl, sig2_sur, beta, psi, 
# #       pii, c, shape_c, scale_c, kappa_psi, phi_psi, cellarea, sigma2, 
# #       COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent, 
# #       adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     postEta = rbind(postEta, dummy$postEta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha = postAlpha[nSamples,]
# #     eta = postEta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       postEta = matrix(postEta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, shape_c, scale_c,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # 
# # ### fix kappa_eta ----
# # 
# # 
# # Rf_fitJoint_cZ_eta_b0t_psit_fixKeta = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam, coord, cellarea, 
# #     pii, kappa_eta, phi_eta, phi_psi, sig2_obl, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>% 
# #     arrange(MARU) %>% 
# #     dplyr::select(lon, lat) %>% 
# #     unique() %>% 
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha = startval$alpha
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha = postEta = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_cZ_eta_b0t_psit_fixKeta(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha, eta, lam0, 
# #       lam1, kappa_eta, phi_eta, sig2_obl, sig2_sur, beta, psi, 
# #       pii, c, shape_c, scale_c, kappa_psi, phi_psi, cellarea, sigma2, 
# #       COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent, 
# #       adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     postEta = rbind(postEta, dummy$postEta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha = postAlpha[nSamples,]
# #     eta = postEta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,1]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       postEta = matrix(postEta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, shape_c, scale_c,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # Rf_fitJoint_cZ_a0t_etat_b0t_psi_ja_fixKeta = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam, coord, cellarea, 
# #     pii, kappa_eta, phi_eta, phi_psi, sig2_obl, tau2, 
# #     shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>% 
# #     arrange(MARU) %>% 
# #     dplyr::select(lon, lat) %>% 
# #     unique() %>% 
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_cZ_a0t_etat_b0t_psi_ja_fixKeta(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha, 
# #       eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur, 
# #       beta, psi, pii, c, shape_c, scale_c, kappa_psi, phi_psi, 
# #       cellarea, sigma2, COValpha, COVbeta, updateCOV, adaptInterval, 
# #       adaptFactorExponent, adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,1]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postPsi = matrix(postPsi[-(1:n.save),], nrow = n.keep)
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, shape_c, scale_c,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, 
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # # -----------------------------------------------------------------------------=
# # ## not centering logZ ----
# # # -----------------------------------------------------------------------------=
# # Rf_fitJoint_a0t_etat_b0t_psi_ja = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale, 
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run, 
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>% 
# #     filter(type == 'obl') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>% 
# #     filter(type == 'sfc') %>% 
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>% 
# #     arrange(MARU) %>% 
# #     dplyr::select(lon, lat) %>% 
# #     unique() %>% 
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)  
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_a0t_etat_b0t_psi_ja(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, 
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha, 
# #       eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur, 
# #       beta, psi, pii, c, shape_c, scale_c, kappa_psi, phi_psi, 
# #       cellarea, sigma2, COValpha, COVbeta, updateCOV, adaptInterval, 
# #       adaptFactorExponent, adapIter) 
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = postEta[[t]][-(1:n.save),]
# #       }
# #       postKappa = postKappa[-(1:n.save),]
# #       postLam = postLam[-(1:n.save),]
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = postBeta[-(1:n.save),]
# #       postPsi = postPsi[-(1:n.save),]
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postKappa, postLam, postSig2, 
# #          postBeta, postPsi, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, 
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale, 
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi, pii, c, kappa_psi, phi_psi, 
# #          cellarea, sigma2, COValpha, COVbeta, adapIter, tau2, shape_c, scale_c,
# #          file = filename)
# #   }
# # }
# # 
# # 
# # Rf_fitJoint_b0 = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = startval$beta
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, pii, c, shape_c, scale_c, cellarea, sigma2, COValpha0, COValpha,
# #       COVbeta, updateCOV, adaptInterval, adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = postAlpha0[-(1:n.save),]
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = postEta[[t]][-(1:n.save),]
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = postLam[-(1:n.save),]
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = postBeta[-(1:n.save),]
# #       postC = postC[-(1:n.save),]
# #     }
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0, alpha0tilde, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta,
# #          pii, c, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # Rf_fitJoint_b0_gp = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = startval$beta
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0_gp(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, psi, pii, c, shape_c, scale_c, kappa_psi, phi_psi, cellarea,
# #       sigma2, COValpha0, COValpha, COVbeta, updateCOV, adaptInterval,
# #       adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = postAlpha0[-(1:n.save),]
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = postEta[[t]][-(1:n.save),]
# #       }
# #       postKappa = postKappa[-(1:n.save),]
# #       postLam = postLam[-(1:n.save),]
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = postBeta[-(1:n.save),]
# #       postPsi = postPsi[-(1:n.save),]
# #       postC = postC[-(1:n.save),]
# #     }
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postC, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # Rf_fitJoint_b0t = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run, startval,
# #     path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, pii, c, shape_c, scale_c, cellarea, sigma2, COValpha0, COValpha,
# #       COVbeta, updateCOV, adaptInterval, adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = matrix(postAlpha0[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta,
# #          pii, c, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # Rf_fitJoint_b0t_pref = function(
#     #     dat.zoop, dat.dist, dat.pam, dat.pref, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_delta, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run, startval,
# #     path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## irregular stations ----
# #   # ---------------------------------------------------------------------------=
# #   Spref = cbind(dat.pref$lon, dat.pref$lat)
# #   
# #   indpref = foreach(i = 1:nrow(dat.pref), .combine = rbind) %do% {
# #     which.min(geodist(coord, dat.pref %>% dplyr::select(lon, lat) %>% slice(i)))
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     mu0 = startval$mu0
# #     kappa_delta = startval$kappa_delta
# #     delta = startval$delta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     COVmu0 = 1
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postMu0 = postDelta = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t_pref(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Spref, bath,
# #       noise, Sdist, Wpam, dur_dist, f_dist_obs, f_dist, p_pam,
# #       indobl, indsur, indpref, indtwhale, inddist, indmdist, indmwhale,
# #       alpha0tilde, alpha0, alpha, eta, mu0, delta, lam0, lam1,
# #       tau2, kappa_eta, phi_eta, kappa_delta, phi_delta, sig2_obl,
# #       sig2_sur, beta, pii, c, shape_c, scale_c, cellarea, sigma2,
# #       COValpha0, COValpha, COVbeta, COVmu0, updateCOV, adaptInterval,
# #       adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postMu0 = rbind(postMu0, dummy$postMu0)
# #     postDelta = rbind(postDelta, dummy$postDelta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     mu0 = postMu0[nSamples,1]
# #     delta = postDelta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     kappa_delta = postKappa[nSamples,2]
# #     beta = postBeta[nSamples,]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     COVmu0 = dummy$COVmu0
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = matrix(postAlpha0[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postMu0 = matrix(postMu0[-(1:n.save),], nrow = n.keep)
# #       postDelta = matrix(postDelta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postMu0, postDelta,
# #          postKappa, postLam, postSig2, postBeta, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, Spref,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indpref, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, mu0, delta,
# #          lam0, lam1, sig2_sur, kappa_eta, phi_eta, kappa_delta, phi_delta, beta,
# #          pii, c, cellarea, sigma2, COValpha0, COValpha, COVbeta, COVmu0,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # Rf_fitJoint_b0t_psi = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t_psi(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, psi, pii, c, shape_c, scale_c, kappa_psi, phi_psi, cellarea,
# #       sigma2, COValpha0, COValpha, COVbeta, updateCOV, adaptInterval,
# #       adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = postAlpha0[-(1:n.save),]
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = postEta[[t]][-(1:n.save),]
# #       }
# #       postKappa = postKappa[-(1:n.save),]
# #       postLam = postLam[-(1:n.save),]
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = postBeta[-(1:n.save),]
# #       postPsi = postPsi[-(1:n.save),]
# #       postC = postC[-(1:n.save),]
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postC, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # Rf_fitJoint_b0t_psit = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t_psit(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, psi, pii, c, shape_c, scale_c, kappa_psi,
# #       phi_psi, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #       updateCOV, adaptInterval, adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = matrix(postAlpha0[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postC, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # Rf_fitJoint_b0t_psit_1day = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.02, 1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t_psit_1day(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha, eta, lam0,
# #       lam1, tau2, kappa_eta, phi_eta, sig2_obl, sig2_sur, beta,
# #       psi, pii, c, shape_c, scale_c, kappa_psi, phi_psi, cellarea,
# #       sigma2, COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent,
# #       adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postC, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # # fit the joint model of Zoop and whale
# # Rf_fitZoop = function(
#     #     dat.zoop, Xzoop,
# #     coord, cellarea, phi_eta, sig2_obl, tau2,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postSig2 = postLam = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitZoop(
# #       outeri, distmat, Xzoop, logYobl, logYsur, indobl, indsur,
# #       alpha0tilde, alpha0, alpha, eta, lam0, lam1, tau2, kappa_eta,
# #       phi_eta, sig2_obl, sig2_sur, sigma2, COValpha0, COValpha,
# #       updateCOV, adaptInterval, adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postLam = rbind(postLam, dummy$postLam)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     kappa_eta = postKappa[nSamples,1]
# #     sig2_sur = postSig2[nSamples,1]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = postAlpha0[-(1:n.save),]
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = postEta[[t]][-(1:n.save),]
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postLam = postLam[-(1:n.save),]
# #     }
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postSig2, postLam,
# #          distmat, Xzoop, logYobl, logYsur, indobl, indsur,
# #          alpha0tilde, alpha0, alpha, eta, lam0, lam1,
# #          kappa_eta, phi_eta, sig2_sur,
# #          sigma2, COValpha0, COValpha, adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # # fit the whale model
# # Rf_fitWhale_b0t_psi = function(
#     #     dat.dist, dat.pam, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_psi, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1)
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postKappa = postBeta = postC = c()
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i] - outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitWhale_b0t_psi(
# #       outeri, distmat, Xwhale, Sdist, Wpam, dur_dist, f_dist_obs,
# #       f_dist, p_pam, inddist, indmdist, indmwhale, beta, psi, pii,
# #       c, shape_c, scale_c, kappa_psi, phi_psi, cellarea, sigma2, COVbeta, updateCOV,
# #       adaptInterval, adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postBeta)
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     c = postC[nSamples,]
# #     kappa_psi = postKappa[nSamples,1]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COVbeta = dummy$COVbeta
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = postPsi[[t]][-(1:n.save),]
# #       }
# #       postC = postC[-(1:n.save),]
# #     }
# #     
# #     save(Accprob, rtime, curniter,
# #          postKappa, postBeta, postC, postPsi,
# #          distmat, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam,
# #          inddist, indmdist, indmwhale,
# #          beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # # =============================================================================-
# # # fit the joint model after some number of iterations ----
# # # =============================================================================-
# # 
# # Rf_fitSeparJoint_b0t = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, sig2_obl, tau2, shape_c, scale_c,
# #     n.iter, n.keep, n.save, n.adap, joiniter, new.run, startval,
# #     path.cpp, filename.separ, filename.joint){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     c = startval$c
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postLam = postSig2 = postBeta = postC = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     
# #     if(dir.exists(filename.joint)){
# #       filename = filename.joint
# #     } else {
# #       filename = filename.separ
# #     }
# #     
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitSeparJoint_b0t(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, pii, c, shape_c, scale_c, cellarea, sigma2,
# #       COValpha0, COValpha, COVbeta, updateCOV, adaptInterval, adaptFactorExponent,
# #       adapIter, curniter, joiniter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postC = rbind(postC, dummy$postC)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     c = postC[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = matrix(postAlpha0[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postC = matrix(postC[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     if(curniter > joiniter){
# #       filename = filename.joint
# #     } else {
# #       filename = filename.separ
# #     }
# #     
# #     save(Accprob, rtime, curniter, joiniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postC,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta,
# #          pii, c, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # # =============================================================================-
# # # fix call rates ----
# # # =============================================================================-
# # 
# # Rf_fitJoint_b0t_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, startval,
# #     path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == cur.n.iter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     cur.n.iter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, pii, c, cellarea, sigma2, COValpha0, COValpha,
# #       COVbeta, updateCOV, adaptInterval, adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = matrix(postAlpha0[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     save(Accprob, rtime, cur.n.iter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta,
# #          pii, c, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # Rf_fitJoint_b0t_psi_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postPsi = postKappa = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t_psi_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, psi, pii, c, kappa_psi, phi_psi, cellarea,
# #       sigma2, COValpha0, COValpha, COVbeta, updateCOV, adaptInterval,
# #       adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = postAlpha0[-(1:n.save),]
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = postEta[[t]][-(1:n.save),]
# #       }
# #       postKappa = postKappa[-(1:n.save),]
# #       postLam = postLam[-(1:n.save),]
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = postBeta[-(1:n.save),]
# #       postPsi = postPsi[-(1:n.save),]
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # Rf_fitJoint_b0t_psit_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t_psit_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, psi, pii, c, kappa_psi, phi_psi, cellarea,
# #       sigma2, COValpha0, COValpha, COVbeta, updateCOV, adaptInterval,
# #       adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = matrix(postAlpha0[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # 
# # # MH proposal covariance matrices
# # Rf_fitZoopJoint_b0t_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run, startval,
# #     sigma2, COValpha0, COValpha, COVbeta, adapIter,
# #     path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == cur.n.iter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     cur.n.iter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_b0t_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, pii, c, cellarea, sigma2, COValpha0, COValpha,
# #       COVbeta, updateCOV, adaptInterval, adaptFactorExponent, adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = matrix(postAlpha0[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     save(Accprob, rtime, cur.n.iter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta,
# #          pii, c, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # Rf_fitSeparJoint_b0t_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, joiniter, new.run, startval,
# #     path.cpp, filename.separ, filename.joint){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha0 = startval$alpha0
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, length(alpha0)+1+1)
# #     COValpha0 = rep(1, length(alpha0))
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha0 = postAlpha = postKappa = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     
# #     if(dir.exists(filename.joint)){
# #       filename = filename.joint
# #     } else {
# #       filename = filename.separ
# #     }
# #     
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitSeparJoint_b0t_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #       Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur,
# #       indtwhale, inddist, indmdist, indmwhale, alpha0tilde, alpha0,
# #       alpha, eta, lam0, lam1, tau2, kappa_eta, phi_eta, sig2_obl,
# #       sig2_sur, beta, pii, c, cellarea, sigma2, COValpha0, COValpha,
# #       COVbeta, updateCOV, adaptInterval, adaptFactorExponent, adapIter,
# #       curniter, joiniter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha0 = rbind(postAlpha0, dummy$postAlpha0)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha0 = postAlpha0[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha0 = dummy$COValpha0
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha0 = matrix(postAlpha0[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     if(curniter > joiniter){
# #       filename = filename.joint
# #     } else {
# #       filename = filename.separ
# #     }
# #     
# #     save(Accprob, rtime, curniter, joiniter,
# #          postAlpha0tilde, postAlpha0, postAlpha, postEta, postKappa, postLam, postSig2,
# #          postBeta,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha0, alpha, eta, lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta,
# #          pii, c, cellarea, sigma2, COValpha0, COValpha, COVbeta,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # # =============================================================================-
# # # Fixing call rates & pref sampling ----
# # # =============================================================================-
# # 
# # Rf_fitJoint_pref_b0t_psit_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, dat.pref, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, phi_delta, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## irregular stations ----
# #   # ---------------------------------------------------------------------------=
# #   Spref = cbind(dat.pref$lon, dat.pref$lat)
# #   
# #   indpref = foreach(i = 1:nrow(dat.pref), .combine = rbind) %do% {
# #     which.min(geodist(coord, dat.pref %>% dplyr::select(lon, lat) %>% slice(i)))
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     mu0 = startval$mu0
# #     kappa_delta = startval$kappa_delta
# #     delta = startval$delta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = startval$beta
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     COVmu0 = 1
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postKappa = postMu0 = postDelta = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_pref_b0t_psit_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Spref, Xwhale,
# #       Sdist, Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl,
# #       indsur, indpref, indtwhale, inddist, indmdist, indmwhale,
# #       alpha0tilde, alpha, eta, mu0, delta, lam0, lam1, tau2, kappa_eta,
# #       phi_eta, kappa_delta, phi_delta, sig2_obl, sig2_sur, beta,
# #       psi, pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha,
# #       COVbeta, COVmu0, updateCOV, adaptInterval, adaptFactorExponent,
# #       adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postMu0 = rbind(postMu0, dummy$postMu0)
# #     postDelta = rbind(postDelta, dummy$postDelta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     mu0 = postMu0[nSamples,1]
# #     delta = postDelta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     kappa_delta = postKappa[nSamples,3]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     COVmu0 = dummy$COVmu0
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postMu0 = matrix(postMu0[-(1:n.save),], nrow = n.keep)
# #       postDelta = matrix(postDelta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postMu0, postDelta,
# #          postKappa, postLam, postSig2,
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, Spref,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indpref, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta,
# #          mu0, delta, kappa_delta, phi_delta,
# #          lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, COVmu0,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # Rf_fitJoint_pref_b0t_psi_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, dat.pref, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, phi_delta, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## irregular stations ----
# #   # ---------------------------------------------------------------------------=
# #   Spref = cbind(dat.pref$lon, dat.pref$lat)
# #   
# #   indpref = foreach(i = 1:nrow(dat.pref), .combine = rbind) %do% {
# #     which.min(geodist(coord, dat.pref %>% dplyr::select(lon, lat) %>% slice(i)))
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     mu0 = startval$mu0
# #     kappa_delta = startval$kappa_delta
# #     delta = startval$delta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = startval$beta
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     COVmu0 = 1
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postKappa = postMu0 = postDelta = postLam = postSig2 = postBeta = postPsi = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_pref_b0t_psi_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Spref, Xwhale,
# #       Sdist, Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl,
# #       indsur, indpref, indtwhale, inddist, indmdist, indmwhale,
# #       alpha0tilde, alpha, eta, mu0, delta, lam0, lam1, tau2, kappa_eta,
# #       phi_eta, kappa_delta, phi_delta, sig2_obl, sig2_sur, beta,
# #       psi, pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha,
# #       COVbeta, COVmu0, updateCOV, adaptInterval, adaptFactorExponent,
# #       adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postMu0 = rbind(postMu0, dummy$postMu0)
# #     postDelta = rbind(postDelta, dummy$postDelta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     mu0 = postMu0[nSamples,1]
# #     delta = postDelta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     kappa_delta = postKappa[nSamples,3]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     COVmu0 = dummy$COVmu0
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postMu0 = matrix(postMu0[-(1:n.save),], nrow = n.keep)
# #       postDelta = matrix(postDelta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postPsi = matrix(postPsi[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postMu0, postDelta,
# #          postKappa, postLam, postSig2,
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, Spref,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indpref, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta,
# #          mu0, delta, kappa_delta, phi_delta,
# #          lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, COVmu0,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # 
# # Rf_fitJoint_pref_b0t_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, dat.pref, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_delta, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## irregular stations ----
# #   # ---------------------------------------------------------------------------=
# #   Spref = cbind(dat.pref$lon, dat.pref$lat)
# #   
# #   indpref = foreach(i = 1:nrow(dat.pref), .combine = rbind) %do% {
# #     which.min(geodist(coord, dat.pref %>% dplyr::select(lon, lat) %>% slice(i)))
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     mu0 = startval$mu0
# #     kappa_delta = startval$kappa_delta
# #     delta = startval$delta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = startval$beta
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     COVmu0 = 1
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postKappa = postMu0 = postDelta = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitJoint_pref_b0t_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Spref, Xwhale,
# #       Sdist, Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl,
# #       indsur, indpref, indtwhale, inddist, indmdist, indmwhale,
# #       alpha0tilde, alpha, eta, mu0, delta, lam0, lam1, tau2, kappa_eta,
# #       phi_eta, kappa_delta, phi_delta, sig2_obl, sig2_sur, beta,
# #       pii, c, cellarea, sigma2, COValpha,
# #       COVbeta, COVmu0, updateCOV, adaptInterval, adaptFactorExponent,
# #       adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postMu0 = rbind(postMu0, dummy$postMu0)
# #     postDelta = rbind(postDelta, dummy$postDelta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     mu0 = postMu0[nSamples,1]
# #     delta = postDelta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     kappa_delta = postKappa[nSamples,2]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     COVmu0 = dummy$COVmu0
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postMu0 = matrix(postMu0[-(1:n.save),], nrow = n.keep)
# #       postDelta = matrix(postDelta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postMu0, postDelta,
# #          postKappa, postLam, postSig2, postBeta,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, Spref,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indpref, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta,
# #          mu0, delta, kappa_delta, phi_delta,
# #          lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta,
# #          pii, c, cellarea, sigma2, COValpha, COVbeta, COVmu0,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # 
# # # =============================================================================-
# # # Mod model & Fixing call rates & pref sampling ----
# # # =============================================================================-
# # 
# # Rf_fitMod_pref_b0t_psit_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, dat.pref, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, phi_delta, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## irregular stations ----
# #   # ---------------------------------------------------------------------------=
# #   Spref = cbind(dat.pref$lon, dat.pref$lat)
# #   
# #   indpref = foreach(i = 1:nrow(dat.pref), .combine = rbind) %do% {
# #     which.min(geodist(coord, dat.pref %>% dplyr::select(lon, lat) %>% slice(i)))
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     mu0 = startval$mu0
# #     kappa_delta = startval$kappa_delta
# #     delta = startval$delta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = startval$beta
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     COVmu0 = 1
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postKappa = postMu0 = postDelta = postLam = postSig2 = postBeta = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     postPsi = sapply(1:ncol(psi), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitMod_pref_b0t_psit_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Spref, Xwhale,
# #       Sdist, Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl,
# #       indsur, indpref, indtwhale, inddist, indmdist, indmwhale,
# #       alpha0tilde, alpha, eta, mu0, delta, lam0, lam1, tau2, kappa_eta,
# #       phi_eta, kappa_delta, phi_delta, sig2_obl, sig2_sur, beta,
# #       psi, pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha,
# #       COVbeta, COVmu0, updateCOV, adaptInterval, adaptFactorExponent,
# #       adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postMu0 = rbind(postMu0, dummy$postMu0)
# #     postDelta = rbind(postDelta, dummy$postDelta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     for(t in 1:ncol(psi)){
# #       postPsi[[t]] = rbind(postPsi[[t]], dummy$postPsi[[t]])
# #     }
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     mu0 = postMu0[nSamples,1]
# #     delta = postDelta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = sapply(1:ncol(psi), function(t) postPsi[[t]][nSamples,])
# #     kappa_psi = postKappa[nSamples,2]
# #     kappa_delta = postKappa[nSamples,3]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     COVmu0 = dummy$COVmu0
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postMu0 = matrix(postMu0[-(1:n.save),], nrow = n.keep)
# #       postDelta = matrix(postDelta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(psi)){
# #         postPsi[[t]] = matrix(postPsi[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postMu0, postDelta,
# #          postKappa, postLam, postSig2,
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, Spref,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indpref, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta,
# #          mu0, delta, kappa_delta, phi_delta,
# #          lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, COVmu0,
# #          adapIter, file = filename)
# #   }
# # }
# # 
# # 
# # Rf_fitMod_pref_b0t_psi_ja_fixC = function(
#     #     dat.zoop, dat.dist, dat.pam, dat.pref, Xzoop, Xwhale,
# #     dur_dist, f_dist_obs, f_dist, p_pam,
# #     coord, cellarea, pii, phi_eta, phi_psi, phi_delta, sig2_obl, tau2, c,
# #     n.iter, n.keep, n.save, n.adap, new.run,
# #     startval, path.cpp, filename){
# #   
# #   # ===========================================================================-
# #   # Zoop ----
# #   # ===========================================================================-
# #   nG = nrow(coord)
# #   distmat = geodist(coord) / 1000
# #   Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## oblique tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.obl = dat.zoop %>%
# #     filter(type == 'obl') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
# #     indt = which(dat.obl$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indobl) = NULL
# #   logYobl = log(dat.obl$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## surface tow ----
# #   # ---------------------------------------------------------------------------=
# #   dat.sur = dat.zoop %>%
# #     filter(type == 'sfc') %>%
# #     dplyr::select(month, day, type, lon, lat, zoop)
# #   
# #   indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
# #     indt = which(dat.sur$day[i] == Tz[,2])
# #     inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
# #     c(inds, indt)
# #   }
# #   rownames(indsur) = NULL
# #   logYsur = log(dat.sur$zoop)
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## irregular stations ----
# #   # ---------------------------------------------------------------------------=
# #   Spref = cbind(dat.pref$lon, dat.pref$lat)
# #   
# #   indpref = foreach(i = 1:nrow(dat.pref), .combine = rbind) %do% {
# #     which.min(geodist(coord, dat.pref %>% dplyr::select(lon, lat) %>% slice(i)))
# #   }
# #   
# #   
# #   # ===========================================================================-
# #   # whale ----
# #   # ===========================================================================-
# #   
# #   Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix
# #   indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # day indicator for whale data
# #   indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## Areal sightings ----
# #   # ---------------------------------------------------------------------------=
# #   Sdist = cbind(dat.dist$lon, dat.dist$lat)
# #   
# #   inddist = foreach(i = 1:nrow(dat.dist), .combine = rbind) %do% {
# #     indt = which(dat.dist$day[i] == Tw[,2])
# #     inds = which.min( geodist(coord, data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])) / 1000 )
# #     c(inds, indt)
# #   }
# #   rownames(inddist) = NULL
# #   
# #   indmdist = sapply(1:nrow(dat.dist), function(i) which(dat.dist$month[i] == unique(Tw[,1]))) # month for areal data
# #   
# #   
# #   # ---------------------------------------------------------------------------=
# #   ## PAM ----
# #   # ---------------------------------------------------------------------------=
# #   Wpam = foreach(i = 1:nrow(Tw)) %do% {
# #     dat.pam %>% 
# #       filter(day == Tw[i,2]) %>% 
# #       dplyr::select(num) %>% 
# #       unlist %>% as.vector
# #   }
# #   
# #   locMARU = dat.pam %>%
# #     arrange(MARU) %>%
# #     dplyr::select(lon, lat) %>%
# #     unique() %>%
# #     as.matrix()
# #   
# #   
# #   # ===========================================================================-
# #   # MCMC ----
# #   # ===========================================================================-
# #   
# #   outers = unique(c(0, seq(n.save, n.iter, by = n.save), n.iter))
# #   
# #   updateCOV = TRUE
# #   adaptFactorExponent = 0.8
# #   adaptInterval = n.adap
# #   
# #   if(new.run){
# #     
# #     # initial parameter values for Zoop model
# #     alpha0tilde = startval$alpha0tilde
# #     alpha = startval$alpha
# #     kappa_eta = startval$kappa_eta
# #     eta = startval$eta
# #     mu0 = startval$mu0
# #     kappa_delta = startval$kappa_delta
# #     delta = startval$delta
# #     lam0 = startval$lam0
# #     lam1 = startval$lam1
# #     sig2_sur = startval$sig2_sur
# #     
# #     # initial parameter values for whale model
# #     beta = c(startval$beta0, startval$beta)
# #     kappa_psi = startval$kappa_psi
# #     psi = startval$psi
# #     
# #     # MH proposal covariance matrices
# #     sigma2 = rep(0.05, 1+1+1)
# #     COValpha = diag(length(alpha))
# #     COVbeta = diag(length(beta))
# #     COVmu0 = 1
# #     adapIter = rep(1, length(sigma2))
# #     
# #     start = 1; Accprob = 0; rtime = 0
# #     postAlpha0tilde = postAlpha = postKappa = postMu0 = postDelta = postLam = postSig2 = postBeta = postPsi = c()
# #     postEta = sapply(1:ncol(eta), function(i) c(), simplify = F)
# #     
# #   } else {
# #     load(filename)
# #     start = which(outers == curniter)
# #   }
# #   
# #   sourceCpp(path.cpp)
# #   
# #   for(i in (start+1):length(outers) ){
# #     outeri = outers[i]-outers[i-1]
# #     curniter = outers[i]
# #     
# #     ptm = proc.time()[3]
# #     dummy = cppf_fitMod_pref_b0t_psi_ja_fixC(
# #       outeri, distmat, Xzoop, logYobl, logYsur, Spref, Xwhale,
# #       Sdist, Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl,
# #       indsur, indpref, indtwhale, inddist, indmdist, indmwhale,
# #       alpha0tilde, alpha, eta, mu0, delta, lam0, lam1, tau2, kappa_eta,
# #       phi_eta, kappa_delta, phi_delta, sig2_obl, sig2_sur, beta,
# #       psi, pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha,
# #       COVbeta, COVmu0, updateCOV, adaptInterval, adaptFactorExponent,
# #       adapIter)
# #     rtime = rtime + proc.time()[3] - ptm
# #     
# #     Accprob = ( Accprob * outers[i-1] + colSums(dummy$Accprob) ) / outers[i]
# #     
# #     postAlpha0tilde = rbind(postAlpha0tilde, dummy$postAlpha0tilde)
# #     postAlpha = rbind(postAlpha, dummy$postAlpha)
# #     for(t in 1:ncol(eta)){
# #       postEta[[t]] = rbind(postEta[[t]], dummy$postEta[[t]])
# #     }
# #     postMu0 = rbind(postMu0, dummy$postMu0)
# #     postDelta = rbind(postDelta, dummy$postDelta)
# #     postKappa = rbind(postKappa, dummy$postKappa)
# #     postLam = rbind(postLam, dummy$postLam)
# #     postSig2 = rbind(postSig2, dummy$postSig2)
# #     postBeta = rbind(postBeta, dummy$postBeta)
# #     postPsi = rbind(postPsi, dummy$postPsi)
# #     
# #     nSamples = nrow(postAlpha)
# #     alpha0tilde = postAlpha0tilde[nSamples,]
# #     alpha = postAlpha[nSamples,]
# #     eta = sapply(1:ncol(eta), function(t) postEta[[t]][nSamples,])
# #     mu0 = postMu0[nSamples,1]
# #     delta = postDelta[nSamples,]
# #     lam0 = postLam[nSamples,1]
# #     lam1 = postLam[nSamples,2]
# #     sig2_sur = postSig2[nSamples,1]
# #     kappa_eta = postKappa[nSamples,1]
# #     beta = postBeta[nSamples,]
# #     psi = postPsi[nSamples,]
# #     kappa_psi = postKappa[nSamples,2]
# #     kappa_delta = postKappa[nSamples,3]
# #     
# #     sigma2 = dummy$sigma2
# #     adapIter = dummy$adapIter
# #     COValpha = dummy$COValpha
# #     COVbeta = dummy$COVbeta
# #     COVmu0 = dummy$COVmu0
# #     
# #     if(nSamples == (n.keep + n.save)){
# #       postAlpha0tilde = matrix(postAlpha0tilde[-(1:n.save),], nrow = n.keep)
# #       postAlpha = matrix(postAlpha[-(1:n.save),], nrow = n.keep)
# #       for(t in 1:ncol(eta)){
# #         postEta[[t]] = matrix(postEta[[t]][-(1:n.save),], nrow = n.keep)
# #       }
# #       postMu0 = matrix(postMu0[-(1:n.save),], nrow = n.keep)
# #       postDelta = matrix(postDelta[-(1:n.save),], nrow = n.keep)
# #       postKappa = matrix(postKappa[-(1:n.save),], nrow = n.keep)
# #       postLam = matrix(postLam[-(1:n.save),], nrow = n.keep)
# #       postSig2 = matrix(postSig2[-(1:n.save),], nrow = n.keep)
# #       postBeta = matrix(postBeta[-(1:n.save),], nrow = n.keep)
# #       postPsi = matrix(postPsi[-(1:n.save),], nrow = n.keep)
# #     }
# #     
# #     
# #     save(Accprob, rtime, curniter,
# #          postAlpha0tilde, postAlpha, postEta, postMu0, postDelta,
# #          postKappa, postLam, postSig2,
# #          postBeta, postPsi,
# #          distmat, Xzoop, logYobl, logYsur, Xwhale, Sdist, Spref,
# #          Wpam, dur_dist, f_dist_obs, f_dist, p_pam, indobl, indsur, indpref, indtwhale,
# #          inddist, indmdist, indmwhale, alpha0tilde, alpha, eta,
# #          mu0, delta, kappa_delta, phi_delta,
# #          lam0, lam1, sig2_sur,
# #          kappa_eta, phi_eta, beta, psi,
# #          pii, c, kappa_psi, phi_psi, cellarea, sigma2, COValpha, COVbeta, COVmu0,
# #          adapIter, file = filename)
# #   }
# # }




