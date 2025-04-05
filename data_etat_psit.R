rm(list=ls())
require(tigris) # counties
require(tidyverse)
require(sf) # st_coordinates
require(sp) # point.in.polygon
require(geodist)
require(foreach)

fold = 'github/' # name of folder

path.data = paste0(fold, 'data/')
ifelse(!dir.exists(path.data), dir.create(path.data, recursive = T), FALSE)

filename = paste0(path.data, 'sim_etat_psit.RData')

path.r = paste0(fold, 'src/RFtns.R')
source(path.r)
# =============================================================================-
# Spatial domain ----
# =============================================================================-

## Create an sf POINTS object ----
x.range = c(-70.73, -70.08)
y.range = c(41.71, 42.13)
x = seq(x.range[1], x.range[2], by = .01)
y = seq(y.range[1], y.range[2], by = .01)
points = expand.grid(x, y)

options(tigris_use_cache = TRUE)
counties_sf = counties(cb = TRUE)
mas.sf = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')
mas = st_coordinates(mas.sf)[,1:2]

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


coord = data.frame(lon = xy[,1], lat = xy[,2])

ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_point(aes(x = lon, y = lat), data = coord) +
  theme_bw()



## Distance matrix ----
nG = nrow(coord)
distmat = geodist(coord) / 1000
dim(distmat)
max(distmat)


## Cell area ----
cellarea = as.numeric(
  geodist( data.frame(lon = x[1], lat = y[1]), data.frame(lon = x[2], lat = y[1]) ) / 1000 *
    geodist( data.frame(lon = x[1], lat = y[1]), data.frame(lon = x[1], lat = y[2]) ) / 1000
)


## Load data containing collection dates, covariates, locations of collection sites, etc ----
load(paste0(path.data, 'collection.RData'))


# =============================================================================-
# Zooplankton data ----
# =============================================================================-

### Collection dates
Tz = date.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix


### Temperature
temp.scaled = (temp - mean(temp)) / sd(temp)
colMeans(temp)
colMeans(temp.scaled)


### Design matrix
Xzoop = foreach(i = 1:nrow(Tz)) %do% {
  cbind(temp.scaled[,i])
}


# -----------------------------------------------------------------------------=
## True parameter values ----
# -----------------------------------------------------------------------------=

alpha0tilde = 5.5
tau2 = 0.04
alpha0 = alpha0tilde + sqrt(tau2) * rnorm(nrow(Tz))
alpha = c(0.2)
kappa_eta = 1
phi_eta = max(distmat)/12
eta = sapply(1:nrow(Tz), function(i) t(chol( kappa_eta * exp( - distmat / phi_eta ) )) %*% rnorm(nG))
sig2_sur = 0.5
lam0 = -0.7
lam1 = 1
sig2_obl = 1

trpar = list(alpha0tilde = alpha0tilde, tau2 = tau2, alpha0 = alpha0, alpha = alpha, 
             kappa_eta = kappa_eta, phi_eta = phi_eta, eta = eta, 
             sig2_obl = sig2_obl, sig2_sur = sig2_sur, lam0 = lam0, lam1 = lam1)




# -----------------------------------------------------------------------------=
## True Zoop surface ----
# -----------------------------------------------------------------------------=
logZ = foreach(i = 1:nrow(Tz), .combine = cbind) %do% {
  alpha0[i] + cbind(Xzoop[[i]]) %*% alpha + eta[,i]
}

range(exp(logZ))
range(logZ)
colMeans(logZ)

# range of real oblique tow observations is (6.5 19265.5)
# range of log of real oblique tow observations is (1.88 9.87)

trpar$logZ = logZ

foreach(i = 1:nrow(Tz), .combine = rbind) %do% {
  data.frame(date = as.Date(date.zoop$date[i]), x = coord[,1], y = coord[,2], value = logZ[,i])
} %>%
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  facet_wrap(~date) +
  scale_fill_viridis_c() +
  theme_bw()




# -----------------------------------------------------------------------------=
## Collection locations ----
# -----------------------------------------------------------------------------=

n.points = 30

dat.zoop = c()
for(i in 1:nrow(Tz)){
  
  sites_loc = Rf_simulateNHPP(cbind(rep(0.1, nG)), coord, x.range, y.range, mas)[[1]] %>% 
    slice(1:n.points)
  
  dat.zoop = rbind(
    dat.zoop, 
    data.frame(
      date = date.zoop$date[date.zoop$day == Tz[i,2]],
      month = Tz[i,1],
      day = Tz[i,2],
      type = 'sfc', 
      lon = sites_loc$lon, 
      lat = sites_loc$lat
    ),
    data.frame(
      date = date.zoop$date[date.zoop$day == Tz[i,2]],
      month = Tz[i,1],
      day = Tz[i,2],
      type = 'obl', 
      lon = sites_loc$lon, 
      lat = sites_loc$lat
    )
  )
}

dat.zoop %>% 
  group_by(date, type) %>% 
  summarise(n = n()) %>% 
  group_by(date) %>% 
  summarise(max = max(n)) %>% 
  ungroup %>% select(max) %>% unlist()



# -----------------------------------------------------------------------------=
## Surface tow obs ----
# -----------------------------------------------------------------------------=

dat.sur = dat.zoop %>% 
  filter(type == 'sfc') %>% 
  dplyr::select(date, month, day, type, lon, lat)  

indsur = foreach(i = 1:nrow(dat.sur), .combine = rbind) %do% {
  indt = which(dat.sur$day[i] == Tz[,2])
  inds = which.min(geodist(coord, dat.sur %>% dplyr::select(lon, lat) %>% slice(i)))
  c(inds, indt)
}
rownames(indsur) = NULL
any(is.na(indsur))


logYsur = foreach(i = 1:nrow(indsur), .combine = c) %do% {
  trpar$lam0 + trpar$lam1 * trpar$logZ[indsur[i,1],indsur[i,2]] + rnorm(1, mean = 0, sd = sqrt(trpar$sig2_sur))
}
dat.sur$zoop = exp(logYsur)


dat.sur %>%
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_point(aes(x = lon, y = lat, color = log(zoop))) +
  facet_wrap(~day) +
  labs(x = 'Longitude', y = 'Latitude', color = expression('logY'^'sur')) +
  scale_colour_viridis_c() +
  theme_bw()




# -----------------------------------------------------------------------------=
## Oblique tow obs ----
# -----------------------------------------------------------------------------=

dat.obl = dat.zoop %>% 
  filter(type == 'obl') %>% 
  dplyr::select(date, month, day, type, lon, lat) 

indobl = foreach(i = 1:nrow(dat.obl), .combine = rbind) %do% {
  indt = which(dat.obl$day[i] == Tz[,2])
  inds = which.min(geodist(coord, dat.obl %>% dplyr::select(lon, lat) %>% slice(i)))
  c(inds, indt)
}
rownames(indobl) = NULL
any(is.na(indobl))

logYobl = foreach(i = 1:nrow(indobl), .combine = c) %do% {
  trpar$logZ[indobl[i,1],indobl[i,2]] + rnorm(1, mean = 0, sd = sqrt(trpar$sig2_obl))
}
dat.obl$zoop = exp(logYobl)


dat.obl %>%
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_point(aes(x = lon, y = lat, color = log(zoop))) +
  facet_wrap(~day) +
  labs(x = 'Longitude', y = 'Latitude', color = expression('logY'^'obl')) +
  scale_colour_viridis_c() +
  theme_bw()


# -----------------------------------------------------------------------------=
## final Zoop data ----
# -----------------------------------------------------------------------------=
dat.zoop = rbind(dat.sur, dat.obl) %>% 
  arrange(month, day)





# =============================================================================-
# Whale data ----
# =============================================================================-

### Collection dates
Tw = date.whale %>% dplyr::select(month, day) %>% unique %>% as.matrix

indmwhale = sapply(1:length(Tw[,1]), function(i) which(Tw[i,1] == unique(Tw[,1]))) # month indicator for whale data
indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2])) # date indicator associated with zoop data


### Design matrix
Xwhale = foreach(i = 1:nrow(Tw)) %do% {
  cbind(scale(bath))
}




# -----------------------------------------------------------------------------=
## True parameter values ----
# -----------------------------------------------------------------------------=
beta0 = rep(-4, nrow(Tw))
beta = c(-0.26, 1.5) # coefficients of bathy, noise, and zoop
kappa_psi = 0.2
phi_psi = max(distmat)/12
psi = sapply(1:nrow(Tw), function(i) t(chol( kappa_psi * exp( - distmat / phi_psi ) )) %*% rnorm(nG))
pii = c(0.34, 0.31, 0.55) # prob of surface time for Feb, Mar, Apr from Ganley et al. (2019)
c = rep(2.46, nrow(Tw)) # daily call rates

trpar$beta0 = beta0
trpar$beta = beta
trpar$kappa_psi = kappa_psi
trpar$phi_psi = phi_psi
trpar$psi = psi
trpar$pii = pii
trpar$c = c


# -----------------------------------------------------------------------------=
## True whale abundance intensity ----
# -----------------------------------------------------------------------------=

loglam = foreach(i = 1:nrow(Tw), .combine = cbind) %do% {
  beta0[i] + Xwhale[[i]] %*% beta[1] + beta[2] * (logZ[,indtwhale[i]] - mean(logZ)) + psi[,i]
}
trpar$loglam = loglam
lam = exp(loglam)
# range(lam)
round(sapply(1:ncol(loglam), function(i) sum(lam[,i]) * cellarea))


foreach(i = 1:nrow(Tw), .combine = rbind) %do% {
  data.frame(
    date = date.whale$date[i],
    lon = coord$lon,
    lat = coord$lat,
    value = loglam[,i]
  )
} %>%
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_tile(aes(x = lon, y = lat, fill = value)) +
  facet_wrap(~date) + 
  scale_fill_viridis_c() +
  theme_bw()



# -----------------------------------------------------------------------------=
## True whale locations ----
# -----------------------------------------------------------------------------=
Strue = Rf_simulateNHPP(lam, coord, x.range, y.range, mas)
# plot(Strue[[6]])
sapply(1:length(Strue), function(i) nrow(Strue[[i]]))
trpar$Strue = Strue


p.dat.Strue = data.frame()
for(t in 1:nrow(Tw)){
  p.dat.Strue = rbind(
    p.dat.Strue,
    data.frame(
      day = Tw[t,2],
      lon = trpar$Strue[[t]][,1],
      lat = trpar$Strue[[t]][,2]
    ))
}

p.dat.Strue %>% 
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_point(aes(x = lon, y = lat), size = 0.5) +
  facet_wrap(~day, ncol = 4) +
  theme_bw()




# -----------------------------------------------------------------------------=
## Distance sampling ----
# -----------------------------------------------------------------------------=

### transects for each day ----
dat.line = list()
for(t in 1:nrow(Tw)){
  dat.line[[t]] = dat.line.can # canonical transects
}


### Aerial sightings ----
dat.dist = foreach(t = 1:nrow(Tw), .combine = rbind) %do% {
  
  dat.distl = foreach(l = 1:nrow(dat.line[[t]]), .combine = rbind) %do% {
    
    distsl = rep(0, nrow(Strue[[t]]))
    for(i in 1:nrow(Strue[[t]])){
      
      if(dat.line[[t]]$west[l] <= Strue[[t]]$lon[i] & Strue[[t]]$lon[i] <= dat.line[[t]]$east[l]){
        distsl[i] = geodist(
          data.frame(lon = Strue[[t]]$lon[i], lat = dat.line[[t]]$lat[l]),
          data.frame(lon = Strue[[t]]$lon[i], lat = Strue[[t]]$lat[i])
        ) / 1000  
        
      } else {
        distsl[i] = Inf
      }
    }
    
    f_d = rep(1, length(distsl))
    f_d[ distsl == Inf ] = 0
    f_d[ distsl > 0.75 ] = exp( - (distsl[distsl > 0.75] - 0.75)^2 )
    p_d = pii[indmwhale[t]] * f_d
    
    dummy = rbinom(length(p_d), 1, p_d)
    
    if(sum(dummy) != 0){
      cbind(Tw[t,1], Tw[t,2], dat.line[[t]][l,], Strue[[t]][(dummy == 1),]) # (month, day, line indicator, line lat, line west, line east, whale lon, whale lat) 
    }
  }
  dat.distl
}
dat.dist = data.frame(dat.dist)
colnames(dat.dist) = c('month', 'day', 'line', 'line_lat', 'line_west', 'line_east', 'lon', 'lat')
rownames(dat.dist) = NULL
Sdist = cbind(dat.dist$lon, dat.dist$lat)


t = 3
ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_point(aes(x = lon, y = lat), data = Strue[[t]]) +
  geom_point(aes(x = lon, y = lat), data = dat.dist %>% filter(day == Tw[t,2]), shape = 1, color = 'red') +
  geom_line(aes(x = lon, y = lat, group = line), linetype = 'dotted',
            data = dat.line[[t]] %>% pivot_longer(cols = c(west, east), values_to = 'lon', names_to = 'point')) +
  theme_bw()



### Detection function at observed locations ----
f_dist_obs = foreach(i = 1:nrow(dat.dist), .combine = c) %do% {
  
  dists = geodist(
    data.frame(lon = dat.dist$lon[i], lat = dat.dist$line_lat[i]),
    data.frame(lon = dat.dist$lon[i], lat = dat.dist$lat[i])
  ) / 1000  
  
  if(dists > 0.75){
    exp( - (dists - 0.75)^2 )
  } else {
    1
  }
}


### detection function across spatial domain ----
f_dist = foreach(t = 1:nrow(Tw)) %do% {
  res = foreach(l = 1:nrow(dat.line[[t]]), .combine = cbind) %do% {
    distsl = rep(0, nG)
    for(i in 1:nG){
      
      if(dat.line[[t]]$west[l] <= coord$lon[i] & coord$lon[i] <= dat.line[[t]]$east[l]){
        distsl[i] = geodist(
          data.frame(lon = coord$lon[i], lat = dat.line[[t]]$lat[l]),
          data.frame(lon = coord$lon[i], lat = coord$lat[i])
        ) / 1000  
        
      } else {
        distsl[i] = Inf
      }
    }
    f_d = rep(1, length(distsl))
    f_d[ distsl > 0.75 ] = exp( - (distsl[distsl > 0.75] - 0.75)^2 )
    f_d
  }
  res
}


data.frame(lon = coord[,1], lat = coord[,2], value = 1 - apply(1-f_dist[[1]], 1, prod)) %>% 
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_tile(aes(x = lon, y = lat, fill = value)) +
  scale_fill_viridis_c() +
  theme_bw()



# -----------------------------------------------------------------------------=
## PAM ----
# -----------------------------------------------------------------------------=

### aerial survey time ----
dur_dist = rep(5, nrow(Tw))



### simulate the number of calls per whale ----
Wtrue = foreach(t = 1:nrow(Tw)) %do% {
  cbind(rpois(nrow(trpar$Strue[[t]]), dur_dist[t] * trpar$c[t]), trpar$Strue[[t]][,1], trpar$Strue[[t]][,2]) # (number of calls, lon, lat)
}


dat.pam = foreach(t = 1:nrow(Tw), .combine = rbind) %do% {
  
  recordk = foreach(
    k = 1:nrow(locMARU), # locMARU: locations of hydrophones
    .combine = rbind) %do% {
    record = 0
    for(i in 1:nrow(Wtrue[[t]])){
      prob_pam = Rf_p_pam(
        median(noise.daily$value), # noise.daily: median ambient noise level per day
        data.frame(lon = locMARU[k,1], lat = locMARU[k,2]),
        data.frame(lon = Wtrue[[t]][i,2], lat = Wtrue[[t]][i,3])
      )
      record = record + sum(rbinom(Wtrue[[t]][i,1], 1, prob_pam))
    }
    c(Tw[t,1], Tw[t,2], k, locMARU[k,1], locMARU[k,2], record) # (month, day, MARU, lon, lat, number of detected calls)
  }
  recordk
}
dat.pam = data.frame(dat.pam)
colnames(dat.pam) = c('month', 'day', 'MARU', 'lon', 'lat', 'num')
rownames(dat.pam) = NULL


Wpam = foreach(i = 1:nrow(Tw)) %do% {
  dat.pam %>% 
    filter(day == Tw[i,2]) %>% 
    dplyr::select(num) %>% 
    unlist %>% as.vector
}


### detection probability and intensity over spatial domain ----
p_pam = foreach(i = 1:nrow(Tw)) %:% 
  foreach(k = 1:nrow(locMARU), 
          .combine = cbind) %do% {
            sapply(1:nrow(coord), function(j) 
              Rf_p_pam(
                median(noise.daily$value),
                data.frame(lon = locMARU[k,1], lat = locMARU[k,2]), 
                data.frame(lon = coord[j,1], lat = coord[j,2])
              )
            )
          }



data.frame(
  lon = coord[,1],
  lat = coord[,2],
  value = 1 - apply(1-p_pam[[1]], 1, prod)
) %>%
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_tile(aes(x = lon, y = lat, fill = value)) +
  scale_fill_viridis_c() +
  theme_bw()



# Save data ----
save(Xzoop, Xwhale, 
     dat.zoop, dat.dist, dat.pam, Wtrue, 
     dat.line.can, dat.line,
     f_dist_obs, f_dist, p_pam,
     dur_dist, coord, cellarea, pii, phi_eta, phi_psi, trpar,
     file = filename)


