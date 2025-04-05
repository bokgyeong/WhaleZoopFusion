rm(list = ls())
require(tigris) # counties
require(tidyverse)
require(foreach)
require(xtable)
require(egg) # ggarrange


# directory and file names ----
fold = 'github/'

path.data = paste0(fold, 'data/')
path.sum = paste0(fold, 'sum/')
path.fig = paste0(fold, 'fig/')
ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)


datai = 'sim_etat_psi'
# datai = 'sim_etat_psit'

fiti = 'Joint_etat_psi'
# fiti = 'Joint_etat_psit'


# =============================================================================-
# load data ----
# =============================================================================-

load(paste0(path.data, datai, '.RData'))


## spatial domain ----
options(tigris_use_cache = TRUE)
counties_sf <- counties(cb = TRUE)
mas.sf = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')


## zoop and whale days ----
zoop.date = data.frame(
  date = format(unique(dat.zoop$date), '%b %d, %Y'),
  date2 = format(unique(dat.zoop$date), '%b %d'),
  day = unique(dat.zoop$day)
)

dat.zoop = dat.zoop %>% 
  dplyr::select(-date) %>% 
  left_join(zoop.date)

Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix # 17 days of Zoop collection
Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix

indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2]))

dat.zoop %>% 
  group_by(type) %>% 
  summarise(n = format(n(), nsmall = 0)) %>% 
  t() %>% 
  xtable() %>% 
  print(booktabs = T, include.colnames = F)





## whale data ----

# p.dat.dist = p.dat.line = p.dat.pam = data.frame()
# for(t in 1:nrow(Tw)){
#   p.dat.dist = rbind(
#     p.dat.dist,
#     dat.dist %>%
#       filter(day == Tw[t,2]) %>%
#       mutate(
#         date = zoop.date$date[zoop.date$day == Tw[t,2]],
#         date2 = zoop.date$date2[zoop.date$day == Tw[t,2]]
#       ) %>% 
#       dplyr::select(date, date2, lon, lat) %>% 
#       unique()
#   )
#   
#   p.dat.line = rbind(
#     p.dat.line,
#     dat.line[[t]] %>%
#       pivot_longer(cols = c(west, east), values_to = 'lon', names_to = 'point') %>%
#       add_column(
#         date = zoop.date$date[zoop.date$day == Tw[t,2]],
#         date2 = zoop.date$date2[zoop.date$day == Tw[t,2]]
#       )
#   )
#   p.dat.pam = rbind(
#     p.dat.pam,
#     dat.pam %>%
#       filter(day == Tw[t,2]) %>%
#       mutate(
#         date = zoop.date$date[zoop.date$day == Tw[t,2]], 
#         date2 = zoop.date$date2[zoop.date$day == Tw[t,2]]
#       )
#   )
# }
# 
# p.dat.dist$date = factor(p.dat.dist$date, levels = zoop.date$date)
# p.dat.line$date = factor(p.dat.line$date, levels = zoop.date$date)
# p.dat.pam$date = factor(p.dat.pam$date, levels = zoop.date$date)
# 
# p.dat.dist$date2 = factor(p.dat.dist$date2, levels = zoop.date$date2)
# p.dat.line$date2 = factor(p.dat.line$date2, levels = zoop.date$date2)
# p.dat.pam$date2 = factor(p.dat.pam$date2, levels = zoop.date$date2)




# =============================================================================-
# summarize inference results ----
# =============================================================================-

load(paste0(path.sum, datai, '_', fiti, '_sum.RData'))


## zoop surface ----
dat.zoop.est = data.frame()
for(i in 1:ncol(mean.logZ)){
  dummy = rbind(
    data.frame(date = zoop.date$date2[i], x = coord[,1], y = coord[,2], value = mean.Z[,i])
  )
  dat.zoop.est = rbind(dat.zoop.est, dummy)
}

dat.zoop.est$date = factor(dat.zoop.est$date, levels = zoop.date$date2)


p.zoop.est = dat.zoop.est %>% 
  # filter(date == zoop.date$date2[6]) %>% 
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(trans = scales::log_trans(base = 10), labels = label_comma(), limits = c(6, 51000)) +
  facet_wrap(~date, nrow = 1) +
  labs(title = '(b) Average zooplankton abundance per cubic meter') +
  theme_bw(base_size = 10) +
  theme(
    # legend.position = 'right',
    # legend.key.width = unit(0.1, 'cm'),
    # legend.key.height = unit(0.4, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.1, 'cm'),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.1, 'cm'),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    plot.title = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill=NA, color=NA),
    strip.text=element_text(size=9)
  )

p.zoop.est




## whale surface ----
dat.whale.intensity.est = data.frame()
for(i in 1:ncol(mean.loglam)){
  dummy = rbind(
    data.frame(date = zoop.date$date2[zoop.date$day == Tw[i,2]], x = coord[,1], y = coord[,2], value = mean.lam[,i])
  )
  dat.whale.intensity.est = rbind(dat.whale.intensity.est, dummy)
}

dat.whale.intensity.est$date = factor(dat.whale.intensity.est$date, levels = zoop.date$date2)


p.whale.est = dat.whale.intensity.est %>%
  # mutate(date = factor(date, levels = unique(dat.zoop.est$date))) %>% 
  add_row(date = unique(dat.zoop.est$date)[c(1, 2, 6)], x = NA, y = NA, value = NA) %>% 
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c(trans = scales::log_trans(base = 10), labels = scales::number_format(accuracy = 0.0001)) +
  facet_wrap(~date, nrow = 1) +
  labs(title = '(c) Whale abundance intensity') +
  theme_bw(base_size = 10) +
  theme(
    # legend.position = 'right',
    # legend.key.width = unit(0.1, 'cm'),
    # legend.key.height = unit(0.4, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.1, 'cm'),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.1, 'cm'),
    legend.title = element_blank(),
    axis.title=element_blank(),
    plot.title = element_text(size = 10),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill=NA, color=NA),
    strip.text=element_text(size=9)
  )

p.whale.est

