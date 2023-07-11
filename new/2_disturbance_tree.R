# 0. setup ----------------------------------------------------------------

library(tidyverse);library(pool);library(zoo);library(pracma)

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

core.id <- tbl(KELuser, "core") %>%
  filter(coretype %in% 1,
         !crossdated %in% c(12, 20:22, 99),
         corestatus %in% c(0, 1),
         !is.na(missing_mm),
         !is.na(missing_years),
         missing_mm < 30,
         missing_years < 30) %>%
  pull(id)

# 2. TREE-LEVEL -----------------------------------------------------------

# 2. 1. data --------------------------------------------------------------

ID <- tbl(KELuser, "core") %>% filter(id %in% core.id) %>%
  inner_join(., tbl(KELuser, "tree") %>% filter(id %in% tree.id), by = c("tree_id" = "id")) %>%
  inner_join(., tbl(KELuser, "plot") %>% filter(id %in% plot.id), by = c("plot_id" = "id")) %>%
  pull(id)

data.list <- distGetData(ID = ID)

# 2. 2. growth reconstruction ---------------------------------------------

data.growth <- growthCalculate(data = data.list, windowLength = 10)

# 2. 3. event detection ---------------------------------------------------

data.event <- eventCalculate(data = data.growth, gapAge = c(5:14), nprol = 7)

# ! close database connection ---------------------------------------------

poolClose(KELuser)
