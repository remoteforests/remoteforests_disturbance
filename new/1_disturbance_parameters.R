# 0. setup ----------------------------------------------------------------

library(tidyverse);library(pool);library(lme4);library(zoo);library(pracma);library(cutpointr)

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
         status %in% 1,
         !treetype %in% c("g", "r", "t"),
         !species %in% c("Lians", "99")) %>%
  pull(id)

core.id <- tbl(KELuser, "core") %>%
  filter(coretype %in% 1,
         !crossdated %in% c(12, 20:22, 99),
         corestatus %in% c(0, 1)) %>%
  pull(id)

# 1. DISTURBANCE PARAMETERS -----------------------------------------------

# 1. 1. CA ----------------------------------------------------------------

# 1. 1. 1. data -----------------------------------------------------------

data.ca <- tbl(KELuser, "tree") %>%
  filter(id %in% tree.id,
         !is.na(crowndiam1_m),
         !is.na(crowndiam2_m)) %>%
  inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
  inner_join(., tbl(KELuser, "plot") %>% filter(id %in% plot.id) %>% select(-stand, -subplot), by = c("plot_id" = "id")) %>%
  left_join(., tbl(KELuser, "spatial_hierarchy"), by = "plotid") %>%
  select(stand, subplot, plotid, lng, lat, date, treeid, species, sp_type, dbh_mm, crowndiam1_m, crowndiam2_m) %>%
  collect() %>%
  mutate(ca_m2 = pi * (crowndiam1_m/2) * (crowndiam2_m/2))

write.table(data.ca, "parameters/data/ca_data.csv", sep = ",", row.names = F, na = "")

data.ca <- read.table("parameters/data/ca_data.csv", sep = ",", header = T, stringsAsFactors = F)

# 1. 1. 2. LMM ------------------------------------------------------------

# 1. 1. 2. 1. coniferous --------------------------------------------------

data.conif <- data.ca %>% filter(sp_type %in% "coniferous")

conif <- lmer(sqrt(ca_m2) ~ dbh_mm + (1|stand) + (1|subplot) + (1|plotid), data = data.conif)

summary(conif)

plot(conif)

plot(DHARMa::simulateResiduals(conif))

plot(effects::allEffects(conif, partial.residuals = TRUE))

performance::model_performance(conif) # R2 (cond.) = 0.718, R2 (marg.) = 0.516, RMSE = 0.888
performance::icc(conif, by_group = TRUE) # plotid = 0.153, subplot = 0.010, stand = 0.254

plot(ncf::spline.correlog(x = data.conif$lng, y = data.conif$lat, z = residuals(conif), resamp = 10, latlon = T))

# 1. 1. 2. 2. broadleaved -------------------------------------------------

data.broad <- data.ca %>% filter(sp_type %in% "broadleaved")

broad <- lmer(sqrt(ca_m2) ~ dbh_mm + (1|stand) + (1|subplot) + (1|plotid), data = data.broad)

summary(broad)

plot(broad)

plot(DHARMa::simulateResiduals(broad))

plot(effects::allEffects(broad, partial.residuals = TRUE))

performance::model_performance(broad) # R2 (cond.) = 0.692, R2 (marg.) = 0.544, RMSE = 1.642
performance::icc(broad, by_group = TRUE) # plotid = 0.133, subplot = 0.074, stand = 0.119 

plot(ncf::spline.correlog(x = data.broad$lng, y = data.broad$lat, z = residuals(broad), resamp = 10, latlon = T))

# 1. 2. AI ----------------------------------------------------------------

# 1. 2. 1. data -----------------------------------------------------------

# data.ai <- tbl(KELuser, "ring") %>%
#   inner_join(., tbl(KELuser, "core") %>% filter(id %in% core.id), by = c("core_id" = "id")) %>%
#   inner_join(., tbl(KELuser, "tree") %>% filter(id %in% tree.id), by = c("tree_id" = "id")) %>%
#   inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
#   inner_join(., tbl(KELuser, "plot") %>% filter(id %in% plot.id), by = c("plot_id" = "id")) %>%
#   select(date, treeid, species, sp_group_dist, core_id, year, incr_mm) %>%
#   collect() %>%
#   arrange(core_id, year) %>%
#   group_by(core_id) %>%
#   mutate(incr_mm = if_else(incr_mm %in% 0, NA_real_, incr_mm),
#          incr_mm = na.approx(incr_mm),
#          pg = priorGrowth(incr_mm, windowLength = 10),
#          fg = followGrowth(incr_mm, windowLength = 10),
#          ai = fg - pg) %>%
#   ungroup()
# 
# write.table(data.ai, "parameters/data/ai_data.csv", sep = ",", row.names = F, na = "")

data.ai <- read.table("parameters/data/ai_data.csv", sep = ",", header = T, stringsAsFactors = F) %>% filter(!ai %in% c("NaN", NA))

# 1. 2. 2. calculation ----------------------------------------------------

ai <- data.ai %>% 
  group_by(sp_group_dist) %>%
  mutate(ai_mm = sd(ai) * 1.25) %>%
  ungroup() %>%
  distinct(., core_id, .keep_all = T) %>%
  group_by(sp_group_dist, ai_mm) %>%
  summarise(n = n()) %>%
  ungroup()

openxlsx::write.xlsx(ai, "parameters/ai.xlsx")

# 1. 3. GAP ---------------------------------------------------------------

# 1. 3. 1. data -----------------------------------------------------------

# release <- tbl(KELuser, "ring") %>%
#   inner_join(., tbl(KELuser, "core"), by = c("core_id" = "id")) %>%
#   inner_join(., tbl(KELuser, "tree"), by = c("tree_id" = "id")) %>%
#   inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
#   inner_join(., tbl(KELuser, "plot"), by = c("plot_id" = "id")) %>%
#   select(sp_group_dist, date, treeid, core_id, year, incr_mm) %>%
#   collect() %>%
#   inner_join(., openxlsx::read.xlsx("parameters/ai.xlsx") %>% select(sp_group_dist, ai_mm), by = "sp_group_dist") %>%
#   arrange(core_id, year) %>%
#   group_by(core_id) %>%
#   mutate(incr_mm = if_else(incr_mm %in% 0, NA_real_, incr_mm),
#          incr_mm = na.approx(incr_mm),
#          pg = priorGrowth(incr_mm, windowLength = 10),
#          fg = followGrowth(incr_mm, windowLength = 10),
#          ai = fg - pg,
#          event = ifelse(row_number() %in% peakDetection(x = ai, threshold = first(ai_mm), nups = 1, mindist = 30), "release", NA),
#          event = ifelse(lead(fg, 7) <= pg, NA, event),
#          event = ifelse(lag(pg, 7) >= fg, NA, event)) %>%
#   ungroup() %>%
#   filter(!is.na(event)) %>%
#   distinct(., date, treeid, event)
# 
# data.gap <- tbl(KELuser, "ring") %>%
#   inner_join(., tbl(KELuser, "core") %>%
#                filter(id %in% core.id,
#                       !is.na(missing_mm),
#                       !is.na(missing_years)),
#              by = c("core_id" = "id")) %>%
#   inner_join(., tbl(KELuser, "tree") %>%
#                filter(id %in% tree.id,
#                       growth %in% c(0, 1)),
#              by = c("tree_id" = "id")) %>%
#   inner_join(., tbl(KELuser, "species_fk"), by = c("species" = "id")) %>%
#   inner_join(., tbl(KELuser, "plot") %>% filter(id %in% plot.id), by = c("plot_id" = "id")) %>%
#   select(date, treeid, species, sp_group_dist, growth, core_id, missing_years, year, incr_mm) %>%
#   collect() %>%
#   arrange(core_id, year) %>%
#   group_by(date, treeid, species, sp_group_dist, growth, core_id) %>%
#   mutate(age = year - min(year) + missing_years + 1,
#          incr_mm = if_else(incr_mm %in% 0, NA_real_, incr_mm),
#          incr_mm = na.approx(incr_mm)) %>%
#   filter(age %in% c(5:14)) %>%
#   summarise(gapGrowth = mean(incr_mm, na.rm = T),
#             n = length(incr_mm[!is.na(incr_mm)])) %>%
#   filter(n %in% 10) %>%
#   ungroup() %>%
#   left_join(., release, by = c("date", "treeid")) %>%
#   mutate(use = ifelse(growth %in% 1 & event %in% "release", 0 , 1))
# 
# write.table(data.gap, "parameters/data/gap_data.csv", sep = ",", row.names = F, na = "")

data.gap <- read.table("parameters/data/gap_data.csv", sep = ",", header = T, stringsAsFactors = F)

# 1. 3. 2. calculation ----------------------------------------------------

set.seed(10)

gap <- data.gap %>%
  filter(use %in% 1) %>%
  cutpointr(., x = gapGrowth, class = growth, subgroup = sp_group_dist,
            pos_class = 1, neg_class = 0, direction = ">=",
            method = minimize_boot_metric, metric = abs_d_sens_spec,
            boot_cut = 10000, boot_stratify = T, use_midpoints = T)

plot(gap) # export to .pdf - 4.72 x 7.48 in

# jpeg("parameters/GAP.jpg", width = 19, height = 12, units = "cm", res = 1000, quality = 100)
# 
# plot(gap)
# 
# dev.off()

openxlsx::write.xlsx(gap %>% select(-data, -roc_curve, -boot), "parameters/gap.xlsx")

# 1. 4. DBH ---------------------------------------------------------------

# 1. 4. 1. data -----------------------------------------------------------

# data.dbh <- tbl(KELuser, "tree") %>%
#   filter(id %in% tree.id,
#          growth %in% c(0, 1),
#          layer %in% c(11, 12, 13)) %>%
#   inner_join(., tbl(KELuser, "plot") %>% filter(id %in% plot.id), by = c("plot_id" = "id")) %>%
#   select(date, treeid, growth, layer, dbh_mm) %>%
#   collect() %>%
#   mutate(growth_use = case_when(
#     growth %in% 1 & layer %in% 11 ~ 1,
#     growth %in% 0 & layer %in% 13 ~ 0))
# 
# write.table(data.dbh, "parameters/data/dbh_data.csv", sep = ",", row.names = F, na = "")

data.dbh <- read.table("parameters/data/dbh_data.csv", sep = ",", header = T, stringsAsFactors = F)

# 1. 4. 2. calculation ----------------------------------------------------

set.seed(10)

dbh <- data.dbh %>%
  filter(!is.na(growth_use)) %>%
  cutpointr(., x = dbh_mm, class = growth_use,
            pos_class = 1, neg_class = 0, direction = ">=",
            method = minimize_boot_metric, metric = abs_d_sens_spec,
            boot_cut = 10000, boot_stratify = T, use_midpoints = T)

plot(dbh) # export to .pdf - 4.72 x 7.48 in

# jpeg("parameters/DBH.jpg", width = 19, height = 12, units = "cm", res = 1000, quality = 100)
# 
# plot(dbh)
# 
# dev.off()

openxlsx::write.xlsx(dbh %>% select(-data, -roc_curve, -boot), "parameters/dbh.xlsx")

# 1. 5. GAPMAKER AGE ------------------------------------------------------

# average age when tree reaches the minimum gapmaker DBH threshold [251 mm] (Frelich 2002)

# 1. 5. 1. data -----------------------------------------------------------

# data.age <- tbl(KELuser, "ring") %>%
#   inner_join(., tbl(KELuser, "core") %>%
#                filter(id %in% core.id,
#                       !is.na(missing_mm),
#                       !is.na(missing_years),
#                       missing_mm < 30,
#                       missing_years < 30),
#              by = c("core_id" = "id")) %>%
#   inner_join(., tbl(KELuser, "tree") %>% filter(id %in% tree.id), by = c("tree_id" = "id")) %>%
#   inner_join(., tbl(KELuser, "plot") %>% filter(id %in% plot.id), by = c("plot_id" = "id")) %>%
#   select(date, treeid, dbh_mm, core_id, missing_mm, missing_years, year, incr_mm) %>%
#   collect() %>%
#   arrange(core_id, year) %>%
#   group_by(core_id) %>%
#   mutate(incr_mm = if_else(incr_mm %in% 0, NA_real_, incr_mm),
#          incr_mm = na.approx(incr_mm),
#          dbh_growth = ifelse(row_number() == 1, incr_mm + missing_mm, incr_mm),
#          dbh_growth = cumsum(dbh_growth) * 2,
#          dbh_mm = ifelse(is.na(dbh_mm), max(dbh_growth), dbh_mm),
#          dbh_coef = max(dbh_mm) / max(dbh_growth),
#          dbh_growth = round(dbh_growth * dbh_coef, 0),
#          age = year - min(year) + missing_years + 1) %>%
#   ungroup()
# 
# write.table(data.age, "parameters/data/age_data.csv", sep = ",", row.names = F, na = "")

data.age <- read.table("parameters/data/age_data.csv", sep = ",", header = T, stringsAsFactors = F)

# 1. 5. 2. calculation ----------------------------------------------------

data.age %>%
  filter(dbh_growth %in% 251) %>%
  arrange(core_id, year) %>%
  group_by(core_id) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  summarise(age_mean = mean(age))

# 115