# 0. setup ----------------------------------------------------------------

library(tidyverse);library(lme4);library(zoo);library(pracma);library(cutpointr);library(sf);library(glmmTMB)

source("C:/Users/Ondrej_Vostarek/Desktop/MVP/DB/scripts/conn.R")

source("1_article_fc.R")

plot.id <- tbl(KELuser, "plot") %>%
  filter(ownership %in% 1,
         !country %in% c("Czech Republic", "Germany"),
         foresttype %in% c("beech", "spruce"),
         plottype %in% c(3, 4),
         census %in% 1,
         !is.na(lng),
         !is.na(lat)) %>%
  pull(id)

tree.id <- tbl(KELuser, "tree") %>%
  filter(dbh_mm >= 100,
         status %in% 1,
         !treetype %in% c("g", "r", "t"),
         !species %in% c("Lians", "99")) %>%
  pull(id)

core.id <- tbl(KELuser, "core") %>%
  filter(coretype %in% 1,
         !crossdated %in% c(12, 20:22, 99),
         corestatus %in% c(0, 1)) %>%
  pull(id)

# 2. DISTURBANCE HISTORY --------------------------------------------------

data.ca <- read.table("parameters/data/ca_data.csv", sep = ",", header = T, stringsAsFactors = F)

data.conif <- data.ca %>% filter(sp_type %in% "coniferous")

conif <- lmer(sqrt(ca) ~ dbh_mm + (1|stand) + (1|subplot) + (1|plotid), data = data.conif)

data.broad <- data.ca %>% filter(sp_type %in% "broadleaved")

broad <- lmer(sqrt(ca) ~ dbh_mm + (1|stand) + (1|subplot) + (1|plotid), data = data.broad)

# exclude plots where dead or damaged cored trees represent more than 5 % of total cored canopy area
EX <- tbl(KELuser, "core") %>% 
  filter(coretype %in% 1) %>%
  inner_join(., tbl(KELuser, "tree") %>% 
               filter(dbh_mm >= 100,
                      treetype %in% "0" & onplot %in% c(1, 2) | treetype %in% c("m", "x"),
                      !species %in% c("Lians", "99")),
             by = c("tree_id" = "id")) %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  inner_join(., tbl(KELuser, "plot") %>% 
               filter(id %in% plot.id,
                      !plotsize %in% 500) %>%
               select(-stand, -subplot),
             by = c("plot_id" = "id")) %>%
  left_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid") %>%
  select(stand, subplot, plotid, treeid, sp_type, status, growth, dbh_mm) %>%
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

# 2. 2. unsampled trees ---------------------------------------------------

unsampled <- tbl(KELuser, "core") %>%
  filter(coretype %in% 1) %>%
  right_join(., tbl(KELuser, "tree") %>% 
               filter(!id %in% local(data.event$tree_id),
                      id %in% tree.id,
                      growth %in% 1,
                      treetype %in% "0" & onplot %in% c(1, 2) | treetype %in% c("m", "x")),
             by = c("tree_id" = "id")) %>%
  inner_join(., tbl(KELuser, "tree") %>% filter(id %in% local(data.event$tree_id)) %>% distinct(., plot_id), by = "plot_id") %>%
  inner_join(., tbl(KELuser, "plot"), by = c("plot_id" = "id")) %>%
  filter((foresttype %in% "spruce" & !stand %in% "Polana" & !is.na(subcore)) | (foresttype %in% "beech" | stand %in% "Polana")) %>%
  select(-stand, -subplot) %>%
  left_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid") %>%
  inner_join(., tbl(KELuser, "species_fk") %>% select(species = id, sp_group_dist, sp_type), by = "species") %>%
  mutate(year = 0,
         event = "unsampled",
         year_min = 0) %>%
  select(sp_group_dist, year, event, year_min, treeid, dbh_mm, species, sp_type, plotid, subplot, stand) %>%
  collect()

# 2. 3. plot-level disturbance --------------------------------------------

data.use <- bind_rows(
  data.event %>%
    inner_join(., tbl(KELuser, "core") %>% select(tree_id, core_id = id), by = "tree_id", copy = T) %>%
    inner_join(., tbl(KELuser, "ring") %>% group_by(core_id) %>% summarise(year_min = min(year)), by = "core_id", copy = T) %>%
    inner_join(., tbl(KELuser, "tree") %>% select(tree_id = id, treeid, plot_id, dbh_mm, species), by = "tree_id", copy = T) %>%
    inner_join(., tbl(KELuser, "species_fk") %>% select(species = id, sp_type), by = "species", copy = T) %>%
    inner_join(., tbl(KELuser, "plot") %>% select(plot_id = id, plotid), by = "plot_id", copy = T) %>%
    left_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid", copy = T) %>%
    select(-tree_id, -core_id, -plot_id),
  unsampled
) %>%
  mutate(ca = case_when(
    sp_type %in% "broadleaved" ~ predict(object = broad, newdata = ., allow.new.levels = T),
    sp_type %in% "coniferous" ~  predict(object = conif, newdata = ., allow.new.levels = T)),
    ca = ca^2)

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

write.table(data.dist, "data_dist.csv", sep = ",", row.names = F, na = "")

# -------------------------------------------------------------------------

# data.dist <- read.table("data_dist.csv", sep = ",", header = T, stringsAsFactors = F)
# 
# # disturbance detection sensitivity threshold
# data.dist %>% distinct(., plotid, ca_min) %>% summarise(ca_min = min(ca_min))
# 
# data.mds <- tbl(KELuser, "tree") %>%
#   inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
#   inner_join(., tbl(KELuser, "dp"), by = "sp_group_dist") %>%
#   distinct(., plot_id, gapmaker_age) %>%
#   inner_join(., tbl(KELuser, "plot") %>% filter(plotid %in% local(data.dist$plotid) & census %in% 1), by = c("plot_id" = "id")) %>%
#   # cutting year according to the average age when tree reaches the minimum gapmaker DBH threshold [251 mm] (Frelich 2002)
#   # maximum year due to the limits of release detection (follow growth[10] + sustained release[7])
#   mutate(year_cut = date - gapmaker_age + 1,
#          year_max = date - 17) %>%
#   select(plotid, year_cut, year_max) %>%
#   collect() %>%
#   right_join(., data.dist, by = "plotid") %>%
#   # cut the chronology by year_cut, year_min, and year_max
#   filter(year > year_cut, year >= year_min, year <= year_max) %>%
#   group_by(plotid) %>%
#   complete(year = (min(year)-30):(max(year)+15), fill = list(ca = 0)) %>%
#   mutate(kde = mdsFun(ca, k = 30, bw = 5, st = 7)) %>% # st = 7.5
#   filter(year %in% c((min(year)+15):(max(year)-15))) %>%
#   ungroup()
# 
# data.peaks <- data.mds %>%
#   group_by(plotid) %>%
#   filter(row_number() %in% peakDetection(x = kde, threshold = 10, nups = 5, mindist = 10)) %>% 
#   ungroup()
# 
# GG <- data.mds %>% distinct(., plotid) %>% pull(plotid)
# 
# mds.gg <- bind_rows(
#   tbl(KELuser, "dist_plot") %>% inner_join(., tbl(KELuser, "plot") %>% filter(plotid %in% GG), by = c("plot_id" = "id")) %>%
#     select(plotid, year, ca = ca_per, kde) %>% mutate(type = "DB") %>% collect(),
#   data.mds %>% select(plotid, year, ca, kde) %>% mutate(type = "NEW")
# )
# 
# peaks.gg <- bind_rows(
#   tbl(KELuser, "dist_plot_event") %>% inner_join(., tbl(KELuser, "dist_plot"), by = c("dist_plot_id" = "id")) %>%
#     inner_join(., tbl(KELuser, "plot") %>% filter(plotid %in% GG), by = c("plot_id" = "id")) %>%
#     select(plotid, year, kde) %>% mutate(type = "DB") %>% collect(),
#   data.peaks %>% select(plotid, year, kde) %>% mutate(type = "NEW")
# )
# 
# library(ggrepel)
# 
# pdf("plots.pdf", width = 18, height = 10, pointsize = 12, onefile = T)
# 
# for (i in GG) {
# 
#   mds <- mds.gg %>% filter(plotid %in% i)
#   peaks <- peaks.gg %>% filter(plotid %in% i)
# 
#   print(
#     ggplot(mds) +
#       geom_histogram(aes(year, weight = ca), binwidth = 10, fill = "grey80", breaks = seq(1600, 2010, 10)) +
#       geom_histogram(aes(year, weight = ca), binwidth = 1, fill = "grey20") +
#       geom_line(aes(year, kde), size = 1, colour = "grey20") +
#       geom_point(data = peaks, aes(year, kde), colour = "#78c2ef", size = 3) +
#       geom_text_repel(data = peaks, aes(year, kde, label = year), colour = "#78c2ef", size = 3) +
#       geom_hline(aes(yintercept = 10), linetype = 2, colour = "grey80") +
#       coord_cartesian(xlim = c(1600, 2020)) + coord_cartesian(ylim = c(0, 100)) +
#       ylab("Canopy area (%)") + xlab("Year") +
#       ggtitle(i) +
#       theme_bw() +
#       facet_wrap(~type)
#     )
# 
#   cat(i)
# 
# }
# 
# dev.off()