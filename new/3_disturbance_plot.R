# 0. setup ----------------------------------------------------------------

library(tidyverse);library(pool);library(zoo);library(pracma);library(lme4)

source("new/pw.R")

source("new/0_disturbance_functions.R")

plot.id <- tbl(KELuser, "plot") %>%
  filter(ownership %in% 1,
         !country %in% c("Czech Republic", "France", "Germany"),
         foresttype %in% c("beech", "spruce"),
         plottype %in% c(3, 4),
         census %in% 1,
         !is.na(lng),
         !is.na(lat)) %>%
  filter(!date %in% 2022) %>%
  pull(id)

tree.id <- tbl(KELuser, "tree") %>%
  filter(dbh_mm >= 100,
         treetype %in% "0" & onplot %in% c(1, 2) | treetype %in% c("m", "x"),
         !species %in% c("Lians", "99")) %>%
  pull(id)

# 3. PLOT-LEVEL -----------------------------------------------------------

# 3. 1. crown area model --------------------------------------------------

data.ca <- tbl(KELuser, "dist_data_ca") %>% select(stand, subplot, plotid, sp_type, dbh_mm, ca_m2) %>% collect()

data.conif <- data.ca %>% filter(sp_type %in% "coniferous")

conif <- lmer(sqrt(ca_m2) ~ dbh_mm + (1|stand) + (1|subplot) + (1|plotid), data = data.conif)

data.broad <- data.ca %>% filter(sp_type %in% "broadleaved")

broad <- lmer(sqrt(ca_m2) ~ dbh_mm + (1|stand) + (1|subplot) + (1|plotid), data = data.broad)

# 3. 2. exclude recently disturbed plots ----------------------------------

X <- tbl(KELuser, "core") %>% filter(coretype %in% 1) %>%
  inner_join(., tbl(KELuser, "tree") %>% filter(id %in% tree.id), by = c("tree_id" = "id")) %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  inner_join(., tbl(KELuser, "plot") %>% filter(id %in% plot.id) %>% select(-stand, -subplot), by = c("plot_id" = "id")) %>%
  inner_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid") %>%
  select(stand, subplot, plotid, sp_type, status, growth, dbh_mm) %>%
  collect() %>%
  mutate(ca = case_when(
          sp_type %in% "broadleaved" ~ predict(object = broad, newdata = ., allow.new.levels = T),
          sp_type %in% "coniferous" ~  predict(object = conif, newdata = ., allow.new.levels = T)),
         ca = ca^2) %>%
  group_by(plotid) %>%
  summarise(ca_alive = sum(ca[status %in% 1 & growth %in% 1]),
            ca_dead = sum(ca[status %in% c(2:4, 11:14, 16:23)])) %>%
  mutate(pct_dead = (ca_dead / (ca_alive + ca_dead)) * 100) %>%
  filter(pct_dead > 5) %>%
  pull(plotid)

# 3. 3. data --------------------------------------------------------------

# 3. 3. 1. sampled trees --------------------------------------------------

data.sampled <- tbl(KELuser, "dist_tree") %>%
  inner_join(., tbl(KELuser, "core") %>% select(core_id = id, tree_id), by = "tree_id") %>%
  inner_join(., tbl(KELuser, "ring") %>% 
               group_by(core_id) %>% summarise(year_min = min(year), year_max = max(year) - 17),
             by = "core_id") %>%
  filter(!(year_min > year_max & year > year_max)) %>%
  inner_join(., tbl(KELuser, "tree"), by = c("tree_id" = "id")) %>%
  distinct(., plot_id, tree_id, year_min, year_max) %>%
  group_by(plot_id) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  inner_join(., tbl(KELuser, "plot") %>% 
               filter(!plotid %in% X) %>% select(-stand, -subplot), 
             by = c("plot_id" = "id")) %>%
  inner_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid") %>%
  inner_join(., tbl(KELuser, "tree"), by = c("tree_id" = "id", "plot_id")) %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  inner_join(., tbl(KELuser, "dist_tree"), by = "tree_id") %>%
  select(stand, subplot, plotid, tree_id, sp_type, dbh_mm, year_min, year_max, event, year) %>%
  collect()
 
# 3. 3. 2. unsampled trees ------------------------------------------------

data.unsampled <- tbl(KELuser, "core") %>%
  filter(coretype %in% 1) %>%
  right_join(., tbl(KELuser, "tree") %>% 
               filter(!id %in% local(data.sampled$tree_id),
                      id %in% tree.id,
                      status %in% 1,
                      growth %in% 1),
             by = c("tree_id" = "id")) %>%
  inner_join(., tbl(KELuser, "tree") %>% filter(id %in% local(data.sampled$tree_id)) %>% distinct(., plot_id), by = "plot_id") %>%
  inner_join(., tbl(KELuser, "plot"), by = c("plot_id" = "id")) %>%
  filter((foresttype %in% "spruce" & !stand %in% "Polana" & !is.na(subcore)) | (foresttype %in% "beech" | stand %in% "Polana")) %>%
  select(-stand, -subplot) %>%
  inner_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid") %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  mutate(year_min = 0,
         year_max = 0,
         event = "unsampled",
         year = 0) %>%
  select(stand, subplot, plotid, tree_id, sp_type, dbh_mm, year_min, year_max, event, year) %>%
  collect()

# 3. 3. 3. sampled + unsampled trees --------------------------------------

data.all <- bind_rows(data.sampled, data.unsampled) %>%
  mutate(ca = case_when(
          sp_type %in% "broadleaved" ~ predict(object = broad, newdata = ., allow.new.levels = T),
          sp_type %in% "coniferous" ~  predict(object = conif, newdata = ., allow.new.levels = T)),
         ca = ca^2)

# 3. 4. bootstrapping -----------------------------------------------------

# # minimum number of available tree records (valid + unsampled) = 12
# data.use %>% 
#   distinct(., plotid, treeid) %>%
#   group_by(plotid) %>%
#   summarise(n = n()) %>%
#   summarise(n = min(n))

data.boot <- data.use %>% 
  distinct(., plotid, treeid) %>%
  slice(rep(1:nrow(.), each = 1000)) %>%
  mutate(rep = rep(1:1000, times = nrow(.) / 1000)) 

set.seed(1)

data.dist.boot <- data.boot %>%
  group_by(plotid, rep) %>%
  # number of samples for bootstrapping = 12 - the minimum number of available tree records (valid + unsampled)
  slice_sample(., n = 12, replace = TRUE) %>%
  ungroup() %>%
  left_join(., data.use %>% select(plotid, treeid, year, event, ca, year_min), by = c("plotid", "treeid")) %>%
  do({
    x <- .
    y <- x %>% group_by(plotid, rep, treeid) %>% filter(year %in% first(year), !event %in% "unsampled") %>% group_by(plotid, rep) %>% 
      # filter(n() >= 10) %>% slice_min(., order_by = year_min, n = 10, with_ties = F) %>% summarise(year_min = max(year_min))
      filter(n() >= 1) %>% summarise(year_min = min(year_min), n_valid = n())
    inner_join(
      x %>% filter(!event %in% "unsampled") %>% group_by(plotid, rep, year) %>% summarise(ca = sum(ca)),
      x %>% group_by(plotid, rep, treeid) %>% filter(year %in% first(year)) %>% group_by(plotid, rep) %>% 
        summarise(ca_total = sum(ca), ca_min = min(ca[!event %in% "unsampled"])) %>% 
        inner_join(., y, by = c("plotid", "rep")),
      by = c("plotid", "rep")
    ) 
  }) %>%
  ungroup() %>%
  mutate(ca = ca * 100 / ca_total,
         ca_min = ca_min * 100 / ca_total)

beg <- min(data.dist.boot$year)
end <- max(data.dist.boot$year)

data.dist <- tibble()

for (p in unique(data.dist.boot$plotid)) {
  
  x <- data.dist.boot %>% filter(plotid %in% p) 
  
  out <- inner_join(
    x %>% select(plotid, rep, year, ca) %>% group_by(plotid, rep) %>% complete(year = beg:end, fill = list(ca = 0)) %>% 
      group_by(plotid, year) %>% summarise(ca = mean(ca)),
    x %>% distinct(., plotid, rep, year_min, ca_min, n_valid) %>% group_by(plotid) %>%
      summarise(year_min = round(mean(year_min), 0), ca_min = mean(ca_min), n_valid = round(mean(n_valid), 0)),
    by = "plotid") %>% ungroup()
  
  data.dist <- bind_rows(data.dist, out)
  
  remove(x, out)
}

# 3. 5. kernel density estimation & peak detection ------------------------

data.kde <- tbl(KELuser, "tree") %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  inner_join(., tbl(KELuser, "dp"), by = "sp_group_dist") %>%
  distinct(., plot_id, gapmaker_age) %>%
  inner_join(., tbl(KELuser, "plot") %>% filter(plotid %in% local(data.dist$plotid) & census %in% 1), by = c("plot_id" = "id")) %>%
  # cutting year according to the average age when tree reaches the minimum gapmaker DBH threshold [251 mm] (Frelich 2002)
  # maximum year due to the limits of release detection (follow growth[10] + sustained release[7])
  mutate(year_cut = date - gapmaker_age + 1,
         year_max = date - 17) %>%
  select(plotid, year_cut, year_max) %>%
  collect() %>%
  right_join(., data.dist, by = "plotid") %>%
  # cut the chronology by year_cut, year_min, and year_max
  filter(year > year_cut, year >= year_min, year <= year_max) %>%
  group_by(plotid) %>%
  complete(year = (min(year)-30):(max(year)+15), fill = list(ca = 0)) %>%
  mutate(kde = kdeFun(ca, k = 30, bw = 5, st = 7)) %>%
  filter(year %in% c((min(year)+15):(max(year)-15))) %>%
  ungroup()

data.peaks <- data.kde %>%
  group_by(plotid) %>%
  filter(row_number() %in% peakDetection(x = kde, threshold = 10, nups = 5, mindist = 10)) %>%
  ungroup()