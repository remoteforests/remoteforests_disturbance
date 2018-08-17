
# 0. Libraries and settings -----------------------------------------------
library(tidyverse); library(DBI); library(pool); library(zoo); library(pracma)
source('0_dist_functions.R')

source('pw.R')


# 1. Select core id’s for which disturbnace shall be calculated -----------
# e.g. select all trees from the slovakia and calculate for it
# tree.id <- tbl(KELuser, 'plot') %>%
#   filter(country %in% 'Slovakia') %>% select(plot_id = id) %>%
#   inner_join(., tbl(KELuser, 'tree'), by ='plot_id') %>%
#   select(tree_id = id) %>%
#   inner_join(., tbl(KELuser, 'core'), by = 'tree_id') %>%
#   pull(id)

tree.id <- tbl(KELuser, 'core') %>% pull(id)


# 2. Download and prepare the data ----------------------------------------
data.list <- distGetData(tree.id = tree.id)

# 3. Calculate the growth change ------------------------------------------
data.growth <- growthCalculate(data = data.list, windowLength = 10)


# 4. Calculate the releases -----------------------------------------------
data.release <- releaseCalculate(data = data.growth,  gapAge = c(5:15), nprol = 7)

# 5. Visualize ------------------------------------------------------------

# growth example
TR <- sample(data.release$core$core_id, 9)

data.release$ring %>%
  filter(core_id %in% TR) %>%
  left_join(., data.release$event %>% select(ring_id, event), by = 'ring_id') %>%
  ggplot()+
  geom_line(aes(year, incr_mm)) +
  geom_point(aes(year, incr_mm, col = event), size = 3) +
  theme_bw() +
  facet_wrap(~core_id) +
  scale_color_discrete(na.value="NA")

# dbh and age of event
data.release$event %>%
  filter(event %in% 'release') %>%
  left_join(., data.release$ring %>% select(ring_id, core_id), by = 'ring_id') %>%
  arrange(age) %>%
  group_by(core_id) %>%
  mutate(nevent = row_number()) %>%
  ungroup() %>%
  select(dist_param, nevent, age, dbh_mm) %>%
  gather(variable, value, -dist_param, -nevent) %>%
  ggplot()+
  geom_boxplot(aes(as.factor(nevent), value, fill = as.factor(dist_param))) +
  facet_wrap(~variable)
