# 0. setup ----------------------------------------------------------------

# R 3.6.3 (2020-02-29)

library(ggrepel) # 0.8.2
library(pool) # 0.1.4.3
library(tidyverse) # 1.3.0 (dplyr 1.0.7, forcats 0.5.0, ggplot2 3.3.5, purr 0.3.4, readr 1.3.1, stringr 1.4.0, tibble 3.0.0, tidyr 1.0.2)
library(RPostgreSQL) # 0.6-2 (DBI 1.1.0)

source("new/pw.R")

# 5. GRAPHS ---------------------------------------------------------------

# 5. 1. data --------------------------------------------------------------

data.all <- tbl(KELuser, "dist_event") %>% mutate(peak = "yes") %>%
  right_join(., tbl(KELuser, "dist_chrono"), by = c("dist_chrono_id" = "id")) %>%
  inner_join(., tbl(KELuser, "dist_plot"), by = c("dist_plot_id" = "id")) %>%
  inner_join(., tbl(KELuser, "plot"), by = c("plot_id" = "id")) %>%
  select(plotid, ncores, type, year, ca_pct, kde, peak) %>%
  collect()

# 5. 2. plotting ----------------------------------------------------------

pdf("new/plots.pdf", width = 18, height = 10, pointsize = 12, onefile = T)

for (p in unique(data.all$plotid)) {

  n <- data.all %>% filter(plotid %in% p) %>% distinct(., ncores) %>% pull()
  
  data.gg <- data.all %>% filter(plotid %in% p)
  
  print(
    ggplot(data.gg) +
      geom_histogram(aes(year, weight = ca_pct), breaks = seq(1590, 2010, 10), fill = "grey80") +
      geom_histogram(aes(year, weight = ca_pct), binwidth = 1, fill = "grey20") +
      geom_line(aes(year, kde), size = 1, colour = "grey20") +
      geom_point(data = data.gg %>% filter(peak %in% "yes"), aes(year, kde), shape = 21, colour = "grey20", fill = "#78c2ef", size = 3) +
      geom_text_repel(data = data.gg %>% filter(peak %in% "yes"), aes(year, kde, label = year), colour = "#78c2ef", size = 3) +
      geom_hline(aes(yintercept = 10), linetype = 2, colour = "grey80") +
      coord_cartesian(xlim = c(1590, 2010)) + coord_cartesian(ylim = c(0, 100)) +
      xlab("Year") + ylab("Canopy area (%)") + 
      ggtitle(paste(p, "(number of cores:", n, ")", sep = " ")) +
      theme_bw() +
      facet_wrap(~type)
    )
  
  cat(p)
  
  remove(n, data.gg, p)
}

dev.off()

# ! close database connection ---------------------------------------------

poolClose(KELuser)
