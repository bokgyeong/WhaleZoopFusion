rm(list = ls())
require(tigris) # counties
require(tidyverse)
require(sf) # st_bbox
# require(sp) # point.in.polygon
require(foreach)
require(xtable)
require(rnaturalearth)
require(egg) # ggarrange
require(scales) # label_comma


# directory and file names ----
fold = 'github/'

path.data = paste0(fold, 'data/')
path.fig = paste0(fold, 'fig/')
ifelse(!dir.exists(path.fig), dir.create(path.fig, recursive = T), FALSE)


datai = 'sim_etat_psi'
# datai = 'sim_etat_psit'


# =============================================================================-
# load data ----
# =============================================================================-

load(paste0(path.data, datai, '.RData'))

## spatial domain ----
options(tigris_use_cache = TRUE)
counties_sf = counties(cb = TRUE)
mas.sf = counties_sf %>% filter(NAME %in% c('Plymouth', 'Barnstable'), STATE_NAME == 'Massachusetts')


## observation days ----
zoop.date = data.frame(
  date = format(unique(dat.zoop$date), '%b %d, %Y'),
  date2 = format(unique(dat.zoop$date), '%b %d'),
  day = unique(dat.zoop$day)
)

dat.zoop = dat.zoop %>% 
  dplyr::select(-date) %>% 
  left_join(zoop.date)

dat.zoop$date = factor(dat.zoop$date, levels = unique(dat.zoop$date))
dat.zoop$date2 = factor(dat.zoop$date2, levels = unique(dat.zoop$date2))

Tz = dat.zoop %>% dplyr::select(month, day) %>% unique %>% as.matrix # 17 days of Zoop collection
Tw = dat.dist %>% dplyr::select(month, day) %>% unique %>% as.matrix

indtwhale = sapply(1:nrow(Tw), function(i) which(Tw[i,2] == Tz[,2]))

dat.zoop %>% 
  group_by(type) %>% 
  summarise(n = format(n(), nsmall = 0)) %>% 
  t() %>% 
  xtable() %>% 
  print(booktabs = T, include.colnames = F)



# =============================================================================-
# CCB and zoop sites ----
# =============================================================================-

## CCB ----

# load world coastlines and US states
world <- ne_countries(scale = "medium", returnclass = "sf")
us_states <- ne_states(country = "united states of america", returnclass = "sf")
proj_albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"


# filter for the Northwest Atlantic region
# You might need to adjust the bounding box coordinates to better fit your area of interest
nw_atlantic_bbox <- st_bbox(c(xmin = -80, xmax = -60, ymin = 35, ymax = 50), crs = st_crs(world))
nw_atlantic_prj <- st_bbox(c(xmin = -80, xmax = -60, ymin = 35, ymax = 50), crs = st_crs(proj_albers))


# transform spatial data to Albers equal-area conic projection
# note: you may need to adjust the parameters for your specific region of interest
world_proj <- st_transform(world, proj_albers)
us_states_proj <- st_transform(us_states, proj_albers)


# plot
ploc <- ggplot() +
  geom_sf(data = world_proj, fill = "antiquewhite", color = "gray") + # World coastline
  geom_sf(data = us_states_proj, fill = NA, color = "black") + # US state boundaries
  # annotate("text", x = 2100000, y = 850002, label = "CCB", size = 3, hjust = "left") + # Annotation for Cape Cod Bay
  annotate("text", x = 2100000, y = 850002, label = "CCB", size = 2, hjust = "left") + # Annotation for Cape Cod Bay
  coord_sf(xlim = c(1100000, 2359794), 
           ylim = c(-900000, 1300201)) + # Adjusting view to bounding box
  # labs(x = '', y = '', title = '(a)')+
  labs(x = '', y = '', title = '')+
  theme_bw() +
  theme(
    plot.margin = margin(r = -1, b = -10, t = -10, l = -10),
    # plot.margin = margin(0, 0, 0, 0),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )


## zoop sites ----

p.zoop.site = dat.zoop %>%
  ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_point(aes(x = lon, y = lat, shape = type, color = type)) +
  scale_shape_manual(values = c(1, 4), labels = c('Oblique', 'Surface')) +
  scale_color_manual(values = c('red', 'blue'), labels = c('Oblique', 'Surface')) +
  facet_wrap(~date2, ncol = 3) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_blank(),
    legend.position = 'bottom',
    legend.key.height = unit(0.1, 'cm'),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0, 'cm'),
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill=NA, color=NA),
    strip.text=element_text(size=9)
  )
# p.zoop.site



## combine ----

p.ccb.zoop.site = ggarrange(ploc, p.zoop.site, nrow = 1, widths = c(1, 3.4))



# =============================================================================-
# whale ----
# =============================================================================-

p.dat.dist = p.dat.line = p.dat.pam = data.frame()
for(t in 1:nrow(Tw)){
  p.dat.dist = rbind(
    p.dat.dist,
    dat.dist %>%
      filter(day == Tw[t,2]) %>%
      mutate(
        date = zoop.date$date[zoop.date$day == Tw[t,2]],
        date2 = zoop.date$date2[zoop.date$day == Tw[t,2]]
      ) %>% 
      dplyr::select(date, date2, lon, lat) %>% 
      unique()
  )
  
  p.dat.line = rbind(
    p.dat.line,
    dat.line[[t]] %>%
      pivot_longer(cols = c(west, east), values_to = 'lon', names_to = 'point') %>%
      add_column(
        date = zoop.date$date[zoop.date$day == Tw[t,2]],
        date2 = zoop.date$date2[zoop.date$day == Tw[t,2]]
      )
  )
  p.dat.pam = rbind(
    p.dat.pam,
    dat.pam %>%
      filter(day == Tw[t,2]) %>%
      mutate(
        date = zoop.date$date[zoop.date$day == Tw[t,2]], 
        date2 = zoop.date$date2[zoop.date$day == Tw[t,2]]
      )
  )
}

p.dat.dist$date = factor(p.dat.dist$date, levels = zoop.date$date)
p.dat.line$date = factor(p.dat.line$date, levels = zoop.date$date)
p.dat.pam$date = factor(p.dat.pam$date, levels = zoop.date$date)

p.dat.dist$date2 = factor(p.dat.dist$date2, levels = zoop.date$date2)
p.dat.line$date2 = factor(p.dat.line$date2, levels = zoop.date$date2)
p.dat.pam$date2 = factor(p.dat.pam$date2, levels = zoop.date$date2)


dat.n.whale = rbind(
  p.dat.dist %>% add_column(Name = 'Observed') %>% dplyr::select(date2, lon, lat, Name)  
) %>% 
  group_by(Name, date2) %>% 
  summarise(n = length(lon)) 

dat.n.whale %>% 
  pivot_wider(
    id_cols = Name,
    names_from = date2,
    values_from = n
  ) %>% 
  xtable() %>% 
  print(include.rownames = F)


# whale locations and the number of calls
p.whale.obs = ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_point(aes(x = lon, y = lat), data = p.dat.dist, size = 0.8, color = 'blue', shape = 4, alpha = 0.8) +
  geom_line(aes(x = lon, y = lat, group = line), data = p.dat.line, alpha = 0.5, linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat, size = num), data = p.dat.pam, shape = 1, color = 'red') +
  facet_wrap(~date2, nrow = 1) +
  labs(size = 'Number of calls') +
  # scale_size_continuous(breaks = c(min(dat.pam$num), 10, 20, 40, max(dat.pam$num)), range = c(0.5, 10)) +
  scale_size_continuous(
    breaks = as.vector(round(quantile(dat.pam$num, probs = c(0, 0.25, 0.50, 0.75, 1)))), 
    range = c(0.2, 7)) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_blank(),
    legend.position = 'bottom',
    legend.key.height = unit(0.1, 'cm'),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0, 'cm'),
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill=NA, color=NA),
    strip.text=element_text(size=9)
  )

p.whale.obs


# =============================================================================-
# whale and Zoop ----
# =============================================================================-

dat.zoop %>% filter(type == 'obl') %>% select(zoop) %>% range

p.obs.obl.whale = ggplot() +
  geom_sf(data = mas.sf) +
  coord_sf(xlim = c(-70.8, -69.95), ylim = c(41.7, 42.13)) +
  geom_point(aes(x = lon, y = lat), data = p.dat.dist, size = 0.4, color = 'blue', shape = 4, alpha = 0.6) +
  geom_point(aes(x = lon, y = lat, size = num), data = p.dat.pam, shape = 1, color = 'red', alpha = 0.6) +
  geom_point(aes(x = lon, y = lat, color = zoop), data = dat.zoop %>% filter(type == 'obl'), size = 0.9) +
  # scale_color_viridis_c() +
  scale_color_viridis_c(trans = scales::log_trans(base = 10), labels = label_comma(), limits = c(6, 51000)) +
  facet_wrap(~date2, nrow = 1) +
  labs(size = 'Number of calls') +
  # scale_size_continuous(breaks = c(min(dat.pam$num), 10, 20, 40, max(dat.pam$num)), range = c(0.5, 10)) +
  scale_size_continuous(
    breaks = as.vector(round(quantile(dat.pam$num, probs = c(0, 0.25, 0.50, 0.75, 1)))), 
    range = c(0.1, 4)) +
  labs(x = '', title = '(a) Observations', color = expression(Y^'obl')) +
  theme_bw(base_size = 10) +
  theme(
    # plot.title = element_blank(),
    plot.title = element_text(size = 10),
    legend.position = 'bottom',
    legend.key.height = unit(0.1, 'cm'),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0, 'cm'),
    # legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill=NA, color=NA),
    strip.text=element_text(size=9)
  )

p.obs.obl.whale

