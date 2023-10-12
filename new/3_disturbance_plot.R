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
  pull(id)

tree.id <- tbl(KELuser, "tree") %>%
  filter(dbh_mm >= 100,
         !status %in% c(0, 10, 15, 99),
         growth %in% c(-1, 1),
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

# 3. 2. data --------------------------------------------------------------

# 3. 2. 1. sampled trees --------------------------------------------------

data.sampled <- tbl(KELuser, "dist_tree") %>%
  inner_join(., tbl(KELuser, "core") %>% select(core_id = id, tree_id), by = "tree_id") %>%
  inner_join(., tbl(KELuser, "ring") %>% 
               group_by(core_id) %>% summarise(year_min = min(year), year_max = max(year) - 17),
             by = "core_id") %>%
  filter(!(year_min > year_max & year > year_max)) %>%
  inner_join(., tbl(KELuser, "tree") %>% filter(plot_id %in% plot.id), by = c("tree_id" = "id")) %>%
  distinct(., plot_id, tree_id, year_min, year_max) %>%
  group_by(plot_id) %>%
  filter(n() >= 8) %>%
  ungroup() %>%
  inner_join(., tbl(KELuser, "plot") %>% select(-stand, -subplot), by = c("plot_id" = "id")) %>%
  inner_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid") %>%
  inner_join(., tbl(KELuser, "tree"), by = c("tree_id" = "id", "plot_id")) %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  inner_join(., tbl(KELuser, "dist_tree"), by = "tree_id") %>%
  select(stand, subplot, plotid, tree_id, sp_type, dbh_mm, year_min, year_max, event, year) %>%
  collect()
 
# 3. 2. 2. unsampled trees ------------------------------------------------

data.unsampled <- tbl(KELuser, "core") %>%
  filter(coretype %in% 1) %>%
  right_join(., tbl(KELuser, "tree") %>% 
               filter(!id %in% local(data.sampled$tree_id),
                      id %in% tree.id),
             by = c("tree_id" = "id")) %>%
  inner_join(., tbl(KELuser, "tree") %>% filter(id %in% local(data.sampled$tree_id)) %>% distinct(., plot_id), by = "plot_id") %>%
  inner_join(., tbl(KELuser, "plot"), by = c("plot_id" = "id")) %>%
  filter((foresttype %in% "spruce" & !stand %in% "Polana" & !is.na(subcore)) | ((foresttype %in% "beech" | stand %in% "Polana") & (status %in% c(1:4) | !is.na(subcore)))) %>%
  select(-stand, -subplot) %>%
  inner_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid") %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  mutate(year_min = 0,
         year_max = 0,
         event = "unsampled",
         year = 0) %>%
  select(stand, subplot, plotid, tree_id, sp_type, dbh_mm, year_min, year_max, event, year) %>%
  collect()

# 3. 2. 3. sampled + unsampled trees --------------------------------------

data.all <- bind_rows(data.sampled, data.unsampled) %>%
  mutate(ca_m2 = case_when(
          sp_type %in% "broadleaved" ~ predict(object = broad, newdata = ., allow.new.levels = T),
          sp_type %in% "coniferous" ~  predict(object = conif, newdata = ., allow.new.levels = T)),
         ca_m2 = ca_m2^2)

# 3. 3. bootstrapping -----------------------------------------------------

data.boot <- data.all %>% 
  distinct(., plotid, tree_id) %>%
  slice(rep(1:nrow(.), each = 1000)) %>%
  mutate(rep = rep(1:1000, times = nrow(.) / 1000)) 

set.seed(1)

data.dist.boot <- data.boot %>%
  group_by(plotid, rep) %>%
  slice_sample(., n = 8, replace = TRUE) %>%
  ungroup() %>%
  left_join(., data.all %>% select(plotid, tree_id, year_min, year_max, event, year, ca_m2), by = c("plotid", "tree_id")) %>%
  do({
    
    x <- .
    
    x1 <- x %>% 
      filter(!event %in% "unsampled") %>% 
      group_by(plotid, rep, year) %>% 
      summarise(ca_event = sum(ca_m2)) %>%
      ungroup()
    
    x2 <- x %>% 
      group_by(plotid, rep, tree_id) %>% 
      filter(year %in% first(year)) %>% 
      group_by(plotid, rep) %>% 
      summarise(ca_total = sum(ca_m2)) %>%
      ungroup()
    
    x3 <- x %>% 
      group_by(plotid, rep, tree_id) %>% 
      filter(!event %in% "unsampled", 
             year %in% first(year)) %>% 
      group_by(plotid, rep) %>% 
      filter(n() >= 1) %>% 
      summarise(year_min = min(year_min),
                year_max = max(year_max),
                ncores = n()) %>%
      ungroup()
    
    inner_join(x1, x2 %>% inner_join(., x3, by = c("plotid", "rep")), by = c("plotid", "rep"))
    
  }) %>%
  mutate(ca_pct = ca_event * 100 / ca_total)

beg <- min(data.dist.boot$year)
end <- max(data.dist.boot$year)

data.dist <- tibble()

for (p in unique(data.dist.boot$plotid)) {
  
  x <- data.dist.boot %>% filter(plotid %in% p) 
  
  out <- inner_join(
    x %>% 
      select(plotid, rep, year, ca_pct) %>% 
      group_by(plotid, rep) %>% 
      complete(year = beg:end, fill = list(ca_pct = 0)) %>% 
      group_by(plotid, year) %>% 
      summarise(ca_pct = mean(ca_pct)) %>%
      ungroup(),
    x %>% 
      distinct(., plotid, rep, year_min, year_max, ncores) %>% 
      group_by(plotid) %>%
      summarise(year_min = round(mean(year_min), 0),
                year_max = round(mean(year_max), 0),
                ncores = round(mean(ncores), 0)) %>%
      ungroup(),
    by = "plotid")
  
  data.dist <- bind_rows(data.dist, out)
  
  remove(x, out)
}

# 3. 4. kernel density estimation & peak detection ------------------------

# 3. 4. 1. full chronologies ----------------------------------------------

data.kde.full <- data.dist %>%
  filter(year >= year_min, year <= year_max) %>%
  group_by(plotid) %>%
  complete(year = (min(year)-30):(max(year)+30), fill = list(ca_pct = 0)) %>%
  mutate(kde = kdeFun(ca_pct, k = 30, bw = 5, st = 7)) %>%
  filter(year %in% c((min(year)+15):(max(year)-15))) %>%
  ungroup() %>%
  mutate(type = "full")

data.peaks.full <- data.kde.full %>%
  group_by(plotid) %>%
  filter(row_number() %in% peakDetection(x = kde, threshold = 10, nups = 5, mindist = 10)) %>%
  ungroup() %>%
  mutate(type = "full")

# 3. 4. 2. cut chronologies -----------------------------------------------

data.kde.cut <- tbl(KELuser, "tree") %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  inner_join(., tbl(KELuser, "dist_param"), by = "sp_group_dist") %>%
  distinct(., plot_id, gapmaker_age) %>%
  inner_join(., tbl(KELuser, "plot") %>% filter(plotid %in% local(data.dist$plotid) & census %in% 1), by = c("plot_id" = "id")) %>%
  mutate(year_cut = date - gapmaker_age + 1) %>%
  select(plotid, year_cut) %>%
  collect() %>%
  right_join(., data.dist, by = "plotid") %>%
  filter(year >= year_min, year <= year_max, year > year_cut) %>%
  group_by(plotid) %>%
  mutate(year_min = min(year)) %>%
  select(-year_cut) %>%
  complete(year = (min(year)-30):(max(year)+30), fill = list(ca_pct = 0)) %>%
  mutate(kde = kdeFun(ca_pct, k = 30, bw = 5, st = 7)) %>%
  filter(year %in% c((min(year)+15):(max(year)-15))) %>%
  ungroup() %>%
  mutate(type = "cut")

data.peaks.cut <- data.kde.cut %>%
  group_by(plotid) %>%
  filter(row_number() %in% peakDetection(x = kde, threshold = 10, nups = 5, mindist = 10)) %>%
  ungroup() %>%
  mutate(type = "cut")

# 3. 4. 3. full + cut chronologies ----------------------------------------

data.kde.all <- bind_rows(data.kde.full, data.kde.cut) %>% filter(!is.na(ncores))

data.peaks.all <- bind_rows(data.peaks.full, data.peaks.cut) %>% filter(!is.na(ncores))

# ! close database connection ---------------------------------------------

poolClose(KELuser)
