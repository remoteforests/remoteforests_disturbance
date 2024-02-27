# 0. setup ----------------------------------------------------------------

# R 3.6.3 (2020-02-29)

# data.table 1.12.8
library(pool) # 0.1.4.3
library(tidyverse) # 1.3.0 (dplyr 1.0.7, forcats 0.5.0, ggplot2 3.3.5, purr 0.3.4, readr 1.3.1, stringr 1.4.0, tibble 3.0.0, tidyr 1.0.2)
library(pracma) # 2.2.9
library(RPostgreSQL) # 0.6-2 (DBI 1.1.0)
library(sf) # 0.9-5 (GEOS 3.8.0, GDAL 3.0.4, PROJ 6.3.1)
library(sp)
library(zoo) # 1.8-7

source("new/pw.R")

source("new/0_disturbance_functions.R")

# 4. STAND-LEVEL ----------------------------------------------------------

# 4. 1. data --------------------------------------------------------------

dist.events <- tbl(KELuser, "dist_event") %>%
  inner_join(., tbl(KELuser, "dist_chrono"), by = c("dist_chrono_id" = "id")) %>%
  right_join(., tbl(KELuser, "dist_plot") %>% filter(type %in% "full"), by = c("dist_plot_id" = "id")) %>%
  inner_join(., tbl(KELuser, "plot") %>% select(plot_id = id, plotid), by = "plot_id") %>%
  inner_join(., tbl(KELuser, "dist_polygons") %>% 
               group_by(stand) %>%
               mutate(nplots = n()) %>%
               ungroup(),
             by = "plotid") %>%
  select(stand, nplots, plotid, plot_id, year, year_min, year_max) %>%
  collect()
  
# 4. 2. bootstrapping -----------------------------------------------------

boot <- dist.events %>%
  distinct(., stand, plotid) %>%
  slice(rep(1:nrow(.), each = 1000)) %>%
  mutate(rep = rep(1:1000, times = nrow(.) / 1000)) 

set.seed(1)

dist.boot <- boot %>%
  group_by(stand, rep) %>%
  slice_sample(., n = 10, replace = TRUE) %>%
  ungroup() %>%
  left_join(., dist.events, by = c("stand", "plotid")) %>%
  do({
    
    x <- .
    
    inner_join(
      x %>% 
        filter(!is.na(year)) %>%
        group_by(stand, rep, year) %>%
        summarise(freq = n() / first(nplots)) %>%
        ungroup(),
      x %>% 
        group_by(stand, rep, plotid) %>% 
        filter(year %in% first(year)) %>% 
        group_by(stand, rep) %>% 
        summarise(year_min = max(year_min),
                  year_max = min(year_max)) %>%
        group_by(stand) %>%
        summarise(year_min = round(mean(year_min), 0),
                  year_max = round(mean(year_max), 0)) %>%
        ungroup(),
      by = "stand")
    
  }) %>%
  filter(year >= year_min, year <= year_max) %>%
  group_by(stand, rep) %>%
  complete(year = (min(year_min)-30):(max(year_max)+30), fill = list(freq = 0)) %>%
  mutate(density_pre = kdeFun(freq, k = 30, bw = 5, st = 7),
         density = rollapply(density_pre, width = 5, FUN = mean, fill = 0)) %>% 
  filter(year %in% c((min(year)+15):(max(year)-15))) %>%
  group_by(stand) %>%
  mutate(year_min = first(year_min[!is.na(year_min)]), 
         year_max = first(year_max[!is.na(year_max)])) %>%
  ungroup()

# write.table(dist.boot, "dist_boot.txt", sep = ",", row.names = F, na = "")
# dist.boot <- read.table("dist_boot.txt", sep = ",", header = T, stringsAsFactors = F)

# 4. 3. peak detection ----------------------------------------------------

dist.peaks <- dist.boot %>%
  group_by(stand, rep) %>%
  filter(row_number() %in% peakDetection(x = density, threshold = 0.00001, nups = 5, mindist = 10),
         year %in% c((min(year_min)+1):(max(year_max)-1))) %>%
  group_by(stand, year_min, year_max, year) %>%
  summarise(freq = n() / 100) %>%
  group_by(stand) %>%
  complete(year = (min(year_min)-10):(max(year_max)+10), fill = list(freq = 0)) %>%
  mutate(freqsmooth = kdeFun(freq, k = 11, bw = 1, st = 7)/10) %>%
  filter(row_number() %in% peakDetection(x = freqsmooth, threshold = 0.00001, nups = 0, mindist = 10)) %>%
  ungroup() %>%
  select(stand, year_min, year_max, peakyear = year)

# write.table(dist.peaks, "dist_peaks.csv", sep = ",", row.names = F, na = "")
# dist.peaks <- read.table("dist_peaks.csv", sep = ",", header = T, stringsAsFactors = F)

# 4. 4. join plot-level events to stand-level peaks -----------------------

dist.events.dt <- data.table::data.table(
  dist.events %>% 
    select(-year_min, -year_max) %>%
    inner_join(., dist.peaks %>% distinct(., stand, year_min, year_max), by = "stand") %>%
    filter(year >= year_min, year <= year_max),
  key = c("stand", "year"))

dist.peaks.dt <- data.table::data.table(dist.peaks %>% mutate(year = peakyear), key = c("stand", "year"))

dist.events.joined <-  data.frame(
  dist.peaks.dt[dist.events.dt,
                list(stand, nplots, year_min, year_max, peakyear, year, plot_id, plotid),
                roll = "nearest"])

# write.table(dist.events.joined, "dist_events_joined.csv", sep = ",", row.names = F, na = "")
# dist.events.joined <- read.table("dist_events_joined.csv", sep = ",", header = T, stringsAsFactors = F)

# 4. 5. patch area & stand size -------------------------------------------

polygons <- st_as_sf(tbl(KELuser, "dist_polygons") %>% collect(), wkt = "geometry") %>%
  left_join(., dist.events.joined %>% distinct(., plotid, plot_id, peakyear), by = "plotid")
  
polygons$area <- as.numeric(st_area(polygons) / 10000)
  
patches <- tibble()

for (s in unique(polygons$stand)) {
  
  data.in <- polygons %>% filter(stand %in% s & !is.na(peakyear)) %>% group_by(peakyear) %>% summarise()
  
  data.out <- tibble()
  
  for (p in data.in$peakyear) {
    
    x <- data.in %>% filter(peakyear %in% p)
    
    l <- list()
    
    for (i in 1:length(x[[2]][[1]])) {
      
      if(length(x[[2]][[1]][[i]]) > 1){
        
        print(paste(s,p,i, sep = "_"))
        
        # a <- Polygon(x[[2]][[1]][[i]][[1]])
        # b <- Polygon(x[[2]][[1]][[i]][[2]])
        # 
        # if(a@hole == T){
        #   
        #   print(paste(s,p,i,1, sep = "_"))
        #   
        #   y <- b
        #   
        # } else {
        #   
        #   print(paste(s,p,i,2, sep = "_"))
        #   
        #   y <- a
        #   
        # }
        # 
        # y <- Polygons(list(y), i); y <- SpatialPolygons(list(y))
        # 
        # l[[i]] <- y

        } else {
        
        y <- Polygon(x[[2]][[1]][[i]])
        
        if(y@hole == T) {
          
          print(paste(s,p,i, sep = "_"))
          
        } else {
          
          y <- Polygons(list(y), i); y <- SpatialPolygons(list(y))
          
          l[[i]] <- y
          
        }  
        
      }
      
    }
    
    sp <- SpatialPolygons(lapply(l, function(x){x@polygons[[1]]}))
    
    sf <- st_set_crs(st_as_sf(sp), 3035) %>% rownames_to_column("npatch") %>% mutate(stand = s, peakyear = p)
    sf$patch_area <- as.numeric(st_area(sf) / 10000)
    
    sf <- st_join(sf, polygons %>% filter(stand %in% s & peakyear %in% p) %>% group_by(plot_id) %>% summarise())
    
    sf <- as.data.frame(sf) %>% select(-geometry)
    
    data.out <- bind_rows(data.out, sf)
    
  }
  
  patches <- bind_rows(patches, data.out)
  
  remove(s, data.in, data.out, p, x, l, i, y, sp, sf)
  
}

# # plot can contribute only once within one peakyear!
# patches %>% group_by(peakyear, plot_id) %>% filter(n() > 1)

dist.patches <- patches %>%
  group_by(stand, peakyear, npatch) %>% 
  mutate(nplots_dist = n()) %>%
  inner_join(., as.data.frame(polygons) %>% 
              distinct(., stand, plotid, area) %>% 
              group_by(stand) %>% 
              summarise(stand_size = sum(area)),
            by = "stand") %>%
  select(stand, stand_size, peakyear, npatch, patch_area, nplots_dist, plot_id)

# write.table(dist.patches, "dist_patches.csv", sep = ",", row.names = F, na = "")
# dist.patches <- read.table("dist_patches.csv", sep = ",", header = T, stringsAsFactors = F)

# 4. 6. proportion of plots disturbed & finalisation ----------------------

dist.stand <- tbl(KELuser, "dist_event") %>%
  inner_join(., tbl(KELuser, "dist_chrono"), by = c("dist_chrono_id" = "id")) %>%
  inner_join(., tbl(KELuser, "dist_plot") %>% filter(type %in% "full"), by = c("dist_plot_id" = "id")) %>%
  select(dist_event_id = id, year, plot_id) %>%
  collect() %>%
  inner_join(., dist.events.joined, by = c("plot_id", "year")) %>%
  inner_join(., dist.patches, by = c("stand", "peakyear", "plot_id")) %>%
  inner_join(., dist.events.joined %>%
               distinct(., stand, nplots, peakyear, plotid) %>%
               group_by(stand, peakyear) %>%
               summarise(plotsprop_dist = n() / first(nplots)) %>%
               ungroup(),
             by = c("stand", "peakyear")) %>%
  select(stand, stand_size, nplots, year_min, year_max, peakyear, npatch, 
         patch_area, plotsprop_dist, nplots_dist, plot_id, dist_event_id) %>%
  arrange(stand, peakyear, npatch) %>%
  mutate(stand_size = round(stand_size, 3),
         patch_area = round(patch_area, 3),
         plotsprop_dist = round(plotsprop_dist, 2))

# ! close database connection ---------------------------------------------

poolClose(KELuser)
