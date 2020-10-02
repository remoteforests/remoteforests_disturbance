
# 0. Libraries and konections ---------------------------------------------
library(tidyverse); library(DBI); library(pool); library(zoo); library(pracma); library(ggrepel); library(cowplot); library(gridExtra); library(sp)
source('disturbance/0_dist_functions.R')

source('disturbance/pw.R')
# function to get tables
tbl_k <- function(...){tbl(KELuser, ...)}

# 1. Make the connection to the data and filter the required --------------

tbl_k('plot') %>%
  select(plot_id = id, foresttype, country, stand) %>%
  inner_join( tbl_k('dist_plot'), by = 'plot_id') %>%
  collect() ->
  dist_plot.df

tbl_k('plot') %>%
  select(plot_id = id, foresttype, country, stand, lng, lat) %>%
  inner_join( tbl_k('dist_plot'), by = 'plot_id') %>%
  inner_join(tbl_k('dist_plot_event') %>% select(id = dist_plot_id), by = 'id') %>%
  collect() ->
  dist_plot_event.df
  

# 2. Calculate the kde with events and other figures ----------------------
dist_plot.df %>%
  group_by(foresttype, country, stand, year) %>%
  summarise(ca = mean(ca_per)) %>%
  group_by(foresttype, country, stand) %>%
  complete(year = c(1650:2020), fill = list(ca = 0)) %>%
  mutate(kde = mdsFun(ca, k = 30, bw = 5, st = 7)) %>%
  ungroup() ->
  data.mds

#- Detect the peacks in disturbances
data.mds %>%
  group_by(foresttype, country, stand) %>%
  do( filter(., row_number() %in% peakDetection(x = kde, threshold = 10, nups = 5,  mindist = 10))) %>%
  ungroup() ->
  data.peaks


# predict the disturbance
dist_plot_event.df %>%
  filter(!is.na(lng)) %>%
  #filter(foresttype == 'beech', country == 'Slovakia', stand == 'Sramkova') %>%
  group_by(foresttype,country, stand) %>%
  filter(length(unique(plot_id)) > 2) %>%
  nest() %>%
  mutate(dist.plot = map(data, ~predictDisStand(.x))) %>%
  select(-data) %>% unnest() ->
  data.dist.pred


# fit the kde based on the plot level peaks
dist_plot_event.df %>%
  group_by(foresttype, country, stand, year) %>%
  summarise(ca = sum(kde)) %>%
  mutate(ca = ca * 100 / sum(ca)) %>%
  group_by(foresttype, country, stand) %>%
  complete(year = c(1650:2020), fill = list(ca = 0)) %>%
  mutate(kde = mdsFun(ca, k = 30, bw = 5, st = 7)) %>%
  ungroup()->
  plot_per_mds

plot_per_mds %>%
  group_by(foresttype, country, stand) %>%
  do( filter(., row_number() %in% peakDetection(x = kde, threshold = 10, nups = 5,  mindist = 10))) %>%
  ungroup() ->
  plot_per_mds.peaks

# fit the kde based on number of plots
dist_plot_event.df %>%
  group_by(foresttype, country, stand, year) %>%
  summarise(ca = n()) %>%
  mutate(ca = ca * 100 / sum(ca)) %>%
  group_by(foresttype, country, stand) %>%
  complete(year = c(1650:2020), fill = list(ca = 0)) %>%
  mutate(kde = mdsFun(ca, k = 30, bw = 5, st = 7)) %>%
  ungroup()->
  plot_n_mds

plot_n_mds %>%
  group_by(foresttype, country, stand) %>%
  do( filter(., row_number() %in% peakDetection(x = kde, threshold = 10, nups = 5,  mindist = 10))) %>%
  ungroup() ->
  plot_n_mds.peaks


# merge all data together
left_join(data.mds %>% group_by(foresttype,country, stand) %>% nest(.key = "mds"),
           data.peaks %>% group_by(foresttype,country, stand) %>% nest(.key = "peaks"), by = c('foresttype','country', 'stand')) %>%
  left_join(dist_plot.df %>% group_by(foresttype,country, stand) %>% nest(.key = "mds.plot"), by = c('foresttype','country', 'stand')) %>%
  left_join(dist_plot_event.df %>% group_by(foresttype,country, stand) %>% nest(.key = "peaks.plot"), by = c('foresttype','country', 'stand')) %>%
  left_join(data.dist.pred %>% group_by(foresttype,country, stand) %>% nest(.key = "pred.dist"), by = c('foresttype','country', 'stand')) %>%
  left_join(plot_per_mds %>% group_by(foresttype,country, stand) %>% nest(.key = "plot_per_mds"), by = c('foresttype','country', 'stand')) %>%
  left_join(plot_per_mds.peaks %>% group_by(foresttype,country, stand) %>% nest(.key = "plot_per_mds.peaks"), by = c('foresttype','country', 'stand')) %>%
  left_join(plot_n_mds %>% group_by(foresttype,country, stand) %>% nest(.key = "plot_n_mds"), by = c('foresttype','country', 'stand')) %>%
  left_join(plot_n_mds.peaks %>% group_by(foresttype,country, stand) %>% nest(.key = "plot_n_mds.peaks"), by = c('foresttype','country', 'stand')) %>%
  unite(id, c('foresttype','country', 'stand'), sep = '/') %>%
  group_by(id) %>%
  nest(.key = "all_d") ->
  data.all.nest

  
# 3. Visualize the results ------------------------------------------------
#-- map the pltos
data.all.nest %>%
  mutate(plot = map2(all_d, id, ~plotDistStand(.x, .y))) ->
  plots.all

# 3. Plot all the data ----------------------------------------------------

pdf("KEL stand level disturbances.pdf", width = 9, height = 7,  pointsize = 12, onefile = T)
for(id in 1:nrow(plots.all)){
  cat(id)
  print(plots.all$plot[id])
}
dev.off()


# 9. Close the connection -------------------------------------------------
poolClose(KELuser)