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

# 2. 1. tree-level disturbance --------------------------------------------

ID <- tbl(KELuser, "core") %>%
  filter(id %in% core.id,
         !is.na(missing_mm),
         !is.na(missing_years),
         missing_mm < 30,
         missing_years < 30) %>%
  inner_join(., tbl(KELuser, "tree") %>% 
               filter(id %in% tree.id,
                      growth %in% 1,
                      treetype %in% "0" & onplot %in% c(1, 2) | treetype %in% c("m", "x")),
             by = c("tree_id" = "id")) %>%
  inner_join(., tbl(KELuser, "plot") %>% 
               filter(id %in% plot.id,
                      !plotsize %in% 500,
                      !plotid %in% EX), 
             by = c("plot_id" = "id")) %>%
  # # minimum number of valid cores = 10
  # # - trade-off between number of cores and number of plots
  # # - set to balance number of plots per forest type
  # distinct(., foresttype, plotid, treeid) %>%
  # group_by(foresttype, plotid) %>%
  # summarise(n = n()) %>%
  # filter(n >= 10) %>%
  # group_by(foresttype) %>%
  # summarise(n = n())
  do({
    x <- .
    inner_join(
      x, x %>% distinct(., plotid, treeid) %>% group_by(plotid) %>% filter(n() >= 10) %>% distinct(., plotid),
      by = "plotid"
    )
  }) %>%
  pull(id)

data.list <- distGetData(ID = ID)

data.growth <- growthCalculate(data = data.list, windowLength = 10)

data.event <- eventCalculate(data = data.growth, gapAge = c(5:14), nprol = 7)