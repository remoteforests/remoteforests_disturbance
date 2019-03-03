
# 0. Libraries and konections ---------------------------------------------
library(tidyverse); library(DBI); library(pool); library(zoo); library(pracma); library(ggrepel); library(cowplot); library(gridExtra); library(sp)
source('0_dist_functions.R')

source('pw.R')

# function to get tables
tbl_k <- function(...){tbl(KELuser, ...)}

# 1. Make the connection to the data and filter the required --------------

# All data
tbl_k('dist_tree') %>% select(dist_param, ring_id, event, dbh_event = dbh_mm, age) %>%
  inner_join( tbl_k('dist_param') %>% select(dist_param = id, dbh_th = dbh_mm, dbh_ca_f), by = 'dist_param') %>%
  inner_join( tbl_k('ring') %>% select(ring_id = id, core_id, year), by = 'ring_id') %>%
  inner_join( tbl_k('ring') %>% group_by(core_id) %>% summarise(year_min = min(year)), by = 'core_id') %>%
  inner_join( tbl_k('core') %>% select(core_id = id, tree_id, missing_mm, missing_years, corestatus, crossdated), by = 'core_id') %>%
  full_join( tbl_k('tree') %>% select(tree_id = id, plot_id, treeid, species, x_m, y_m, dbh_mm, growth, treetype, status), by = 'tree_id') %>%
  inner_join( tbl_k('plot') %>% select(plot_id = id, plotid, country, location, stand, foresttype), by = 'plot_id') %>%
  inner_join( tbl_k('species_fk') %>% rename(species = id), by = 'species') %>%
  filter(foresttype == 'spruce' & dbh_mm >= 100 | foresttype == 'beech' & dbh_mm >= 60) %>%
  mutate( year = if_else(event != 'release', year - age, year),
          year_min = year_min - missing_years) %>%
  collect() %>%
  filter(dbh_event <= dbh_th | is.na(dbh_event)) %>%
  mutate(missing_mm = ifelse(missing_mm %in% NA, -1, missing_mm),
         dist_use = case_when(
           growth %in% c(1,99,-1) & treetype %in% c('x', 'm', '0') & missing_mm < 30 & !corestatus %in% c(2,3) & !crossdated %in% c(20:30) & !is.na(event) ~ 'yes',
           TRUE ~ 'no'),
         species = factor(sp_group_dist, levels = c('Picea', 'Fagus', 'Abies', 'Acer', 'Pinus', 'Others'))) ->
  data.all

# Disturbance history data
data.all %>% 
  filter(dist_use == 'yes') %>%
  rowwise() %>%
  mutate(ca = eval(parse(text = dbh_ca_f))) %>% 
  ungroup() %>%
  #select(country, location, stand,foresttype, plotid, plot_id, tree_id, species = sp_group_dist, event, year, dbh_mm,  dbh_ca_f) %>%
  do({
    x <- .
    inner_join(
      x %>% group_by(country, location, stand, foresttype, plot_id, plotid, species, event, year) %>% summarise(ca = sum(ca)),
      x %>% distinct(tree_id, .keep_all = T) %>% group_by(plot_id) %>% summarise(ca_f = sum(ca), n = n()) %>% filter(n >= 5),
      by = 'plot_id'
    ) 
  }) %>%
  ungroup() %>%
  mutate(ca = ca * 100 / ca_f) %>%
  unite(plotid, c('foresttype','country', 'location', 'stand', 'plotid'), sep = "/") %>%
  arrange(plotid, year) %>%
  filter(year %in% c(1600:2010)) %>%
  select(plotid, species, event, year, ca) %>%
  gather(plot_type, value, -year, -ca, -plotid) %>%
  mutate(plot_type = factor(plot_type, levels = c('species', 'event')),
         Species = factor(value, levels = c('Picea', 'Fagus', 'Abies', 'Acer', 'Pinus', 'Others', 'gap', 'release', 'no event')))->
  data.dist

#- Calculate the MDS and the moving sums
data.dist %>%
  filter(plot_type == 'species') %>%
  group_by(plotid, year) %>%
  summarise(ca = sum(ca)) %>%
  group_by(plotid) %>%
  complete(year = 1600:2030, fill = list(ca = 0)) %>%
  mutate(value = mdsFun(ca, k = 30, bw = 5, st = 7)) %>%
  #select(-ca) %>%
  ungroup() %>%
  filter(year %in% c(1600:2010)) ->
  data.mds

#- Detect the peacks in disturbances
data.mds %>%
  group_by(plotid) %>%
  do({
    x <- .
    bind_rows(
      x %>% filter(row_number() %in% peakDetection(x = value, threshold = 10, nups = 10,  mindist = 20)) %>% mutate(method = '10_20_10'),
      x %>% filter(row_number() %in% peakDetection(x = value, threshold = 10, nups = 5,  mindist = 10)) %>% mutate(method = '10_10_5')
    )
  }) %>%
  ungroup() ->
  data.peaks

#- Growth trend data
tbl_k('ring') %>%
  inner_join(.,
             data.all %>% unite(plotid, c('foresttype','country', 'location', 'stand', 'plotid'), sep = "/") %>% filter(!is.na(year_min)) %>% select(plotid, core_id, species, dist_use),
             by = 'core_id', copy = T) %>%
  filter(year %in% c(1700:2000)) %>%
  group_by(plotid, species, dist_use, year) %>%
  summarise(incr_mm = mean(incr_mm)) %>%
  collect() ->
  growth.trend

#- Data for the ploting the tree position map
data.all %>%
  unite(plotid, c('foresttype','country', 'location', 'stand', 'plotid'), sep = "/") %>%
  mutate(status = cut(status, c(-Inf, 0, 9, Inf), c('stump', 'alive', 'dead'))) %>%
  select(plotid, tree_id, x_m, y_m, status, species, dbh_mm, event, year, year_min, dist_use) %>%
  mutate(Species = factor(species, levels = c('Picea', 'Fagus', 'Abies', 'Acer', 'Pinus', 'Others', 'gap', 'release', 'no event'))) ->
  data.position

data.position %>%
  distinct(tree_id, .keep_all = T) ->
  data.age_dbh

#- Data for the ploting predicted year of the event
data.position %>%
  filter(dist_use == 'yes',
         !is.na(x_m))%>%
  #filter(plotid %in% c('spruce_Slovakia_Oravske Beskydy_Pilsko_SLO_PIL_146', 'beech_Slovakia_High Fatra_Kornietova_SLO_KOR_006_1'))%>%
  group_by(plotid) %>%
  nest() %>%
  mutate(dist.plot = map(data, ~predictDisYear(.x))) %>%
  select(-data) ->
  data.dist.pred

#- Information about the disturbances
data.peaks %>%
  filter(method == '10_10_5') %>%
  group_by(plotid) %>%
  summarise(dis_year = year[which.max(value)],
            dis_sev = max(value),
            dist_n = n(),
            dist_mean_sev = mean(value),
            dist_first_year = min(year),
            dist_first_sev = value[which.min(year)],
            dist_last_year = max(year),
            dist_last_sev = value[which.max(year)]) %>%
  mutate_at(vars(dis_sev, dist_mean_sev, dist_first_sev, dist_last_sev), funs(round(.,0))) %>%
  group_by(plotid) %>%
  nest(.key = data.dist.info) ->
  data.dist.info

#-- Calculate the current disturbance rate (up to date)
data.dist.curr <- tbl_k('plot') %>%
  filter(census %in% 1) %>%
  select(plot_id = id, year = date, plotid, country, foresttype, location, stand) %>%
  inner_join(., tbl_k('tree') %>% select(tree_id=id, plot_id, species, onplot, dbh_mm, growth, decay), by = 'plot_id') %>%
  filter(onplot %in% c(1:3), growth %in% c(1,-1,99)) %>%              
  inner_join(.,tbl_k('species_fk') %>% select( species=id, sp_group_dist), by="species") %>%
  inner_join(., tbl_k('dist_group'), by=c("country","foresttype","location","sp_group_dist"))%>% 
  inner_join(., tbl_k('dist_param')  %>% select (dist_param=id, Tdbh_mm=dbh_mm,dbh_ca_f), by="dist_param") %>% 
  filter(dbh_mm >= Tdbh_mm) %>%
  select(year, plot_id, plotid, country, foresttype, location, stand, dbh_mm, dbh_ca_f, decay) %>%
  collect() %>%
  unite(plotid, c('foresttype','country', 'location', 'stand', 'plotid'), sep = "/") %>%
  rowwise() %>% 
  mutate ( ca = eval(parse(text = dbh_ca_f))) %>%
  group_by (year, plotid) %>%
  summarise(disTot=sum (ca),
            rDis1_3N = sum (ca[decay%in%c(1:3)])) %>%
  mutate(ca = rDis1_3N/disTot * 100) %>%
  select(plotid, year, ca)


#-- Join all the data together
inner_join(data.dist %>% group_by(plotid) %>% nest(.key = dist),
           data.mds %>% group_by(plotid) %>% nest(.key = mds), by = 'plotid') %>%
  inner_join(data.peaks %>% group_by(plotid) %>% nest(.key = peaks), by = 'plotid') %>%
  inner_join(growth.trend %>% group_by(plotid) %>% nest(.key = growth.trend), by = 'plotid') %>%
  inner_join(data.position %>% group_by(plotid) %>% nest(.key = position), by = 'plotid') %>%
  inner_join(data.age_dbh %>% group_by(plotid) %>% nest(.key = age_dbh), by = 'plotid') %>%
  left_join(data.dist.pred, by = 'plotid') %>%
  left_join(data.dist.info, by = 'plotid') %>%
  left_join(data.dist.curr %>% group_by(plotid) %>% nest(.key = data.dist.curr), by = 'plotid') %>%
  group_by(plotid) %>%
  nest(.key = all_d) ->
  data.all.nest

#-- map the pltos
data.all.nest %>%
  mutate(plot = map2(all_d, plotid, ~plotDistPos(.x, .y))) ->
  plots.all

# 3. Plot all the data ----------------------------------------------------

pdf("disturbance/KEL plot level disturbances.pdf", width = 18, height = 10,  pointsize = 12, onefile = T)
for(id in 1:nrow(plots.all)){
  cat(id)
  print(plots.all$plot[id])
}
dev.off()


# 9. Close the connection -------------------------------------------------
poolClose(KELuser)
