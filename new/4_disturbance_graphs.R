# 0. setup ----------------------------------------------------------------

library(tidyverse);library(pool);library(ggrepel)

source("new/pw.R")

# 4. GRAPHS ---------------------------------------------------------------

# 4. 1. data --------------------------------------------------------------

data.all <- tbl(KELuser, "dist_event") %>% mutate(peak = "yes") %>%
  right_join(., tbl(KELuser, "dist_chrono"), by = c("dist_chrono_id" = "id")) %>%
  inner_join(., tbl(KELuser, "dist_plot"), by = c("dist_plot_id" = "id")) %>%
  inner_join(., tbl(KELuser, "plot"), by = c("plot_id" = "id")) %>%
  select(plotid, ncores, type, year, ca_pct, kde, peak) %>%
  collect()

# 4. 2. plotting ----------------------------------------------------------

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
