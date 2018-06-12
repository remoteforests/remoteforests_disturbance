distGetData <- function(tree.id = NULL){
  #' @description function download the data from the database for the specified tree id
  #' @param tree.id the vector of core id to be used for calcuation
  
  core.tbl <- tbl(KELuser, 'core') %>%
    filter(id %in% tree.id) %>%
    inner_join(., tbl(KELuser, 'tree') %>% select(tree_id = id, dbh_mm, plot_id, species), by = 'tree_id') %>%
    inner_join(., tbl(KELuser, 'plot') %>% select(plot_id = id, country, foresttype, location), by = 'plot_id') %>%
    inner_join(., tbl(KELuser, 'species_fk') %>% select(species = id, sp_group_dist), by = 'species') %>%
    inner_join(., tbl(KELuser, 'dist_group') %>% select(-id), by = c('country', 'foresttype', 'location', 'sp_group_dist')) %>%
    select(dist_param, plot_id, tree_id, core_id = id, missing_mm, missing_years, dbh_mm) %>%
    collect()
  
  ring.tbl <- tbl(KELuser, 'ring') %>%
    filter(core_id %in% core.tbl$core_id) %>%
    select(core_id, ring_id = id, year, incr_mm) %>%
    collect() %>%
    group_by(core_id) %>%
    mutate(incr_mm = if_else(incr_mm %in% 0, NA_real_, incr_mm),
           incr_mm = na.approx(incr_mm)) %>%
    ungroup()
  
  dist_param.id <- c(unique(core.tbl$dist_param))
  dist_param.tbl <- tbl(KELuser, 'dist_param') %>%
    collect() %>%
    select(dist_param = id, ai_mm, gap_mm, dbh_mm, dbh_ca_f) %>%
    inner_join(., core.tbl %>% distinct(dist_param), by = 'dist_param')
  
  data <- list(dist_param = dist_param.tbl, core = core.tbl, ring = ring.tbl)
  
  return(data)
}

priorGrowth <- function(x, windowLength = 10){
  rollapply( x, 
             width = windowLength,
             FUN = mean,
             fill = NA,
             align = "right",
             na.rm = T,
             partial = TRUE)
}

followGrowth <- function(x, windowLength = 10){
  rollapply( lead(x, 1), 
             width = windowLength,
             FUN = mean,
             fill = NA,
             align = "left",
             na.rm = T,
             partial = TRUE)
}

peakDetection <- function(x, threshold, mindist = 20, nups = 2){
  #' @description identify the index of year when release event occur
  #' @param x a vector of absolute increase change
  #' @param threshold a minimum ai value in mm
  #' @param mindist  minimum distance between two consecutive peaks in years
  #' @param nups number of increasing steps before the peak
  
  x <- ifelse(is.na(x), -0.2, x)
  
  x <- findpeaks(x, 
                 minpeakheight = threshold,
                 minpeakdistance = mindist,
                 nups = nups) 
  
  if(is.null(x)){
    NA
  }else{
    matrix(x, ncol = 4)[,2]
  }
}

keepRelease <- function(year, type, n = 20){
  #' @description calculate the distance between gap origin and releases
  #' @param year the vector of years for event
  #' @param type type of the event (release or gap)
  #' @param n number of years to be checked
  
  keep <- rep('yes', length(year))
  
  if(any(type %in% 'gap')){
    diffyear <- year - year[type %in% 'gap']
    keep[diffyear < n & type %in% 'release'] <- 'no'
  }
  keep
}

growthCalculate <- function(data = data, windowLength = 10){
  #' @description take the list of data prepared by 'dist_get_data' function and calculate the growth change, plus age and dbh of the trees
  #' @param data a list of tree tables
  #' @param windowLength the length of the window for ai calculation
  
  # data quality check
  options(error = NULL) # not to enter debug mode
  
  # perform the checks
  if(!is.list(data)) stop('The input data is not a list of three tables')
  if(!identical(c('core',"dist_param","ring"), ls(data))) stop('The input data tables dont match with required')
  
  # calculate the age, dbh, and the growth change
  inner_join(
    data$ring,
    data$core,
    by = 'core_id'
  ) %>%
    arrange(core_id, year) %>%
    group_by(core_id) %>%
    mutate(dbh_growth = ifelse(row_number() == 1, sum(incr_mm, missing_mm, na.rm = T), incr_mm),
           dbh_growth = cumsum(dbh_growth) * 2,
           dbh_mm = ifelse(is.na(dbh_mm), max(dbh_growth), dbh_mm),
           dbh_coef = max(dbh_mm) / max(dbh_growth),
           dbh_growth = dbh_growth * dbh_coef,
           age = sum(year - min(year), missing_years, 1, na.rm = T),
           pg = priorGrowth(incr_mm, windowLength = windowLength),
           fg = followGrowth(incr_mm, windowLength = windowLength),
           ai = fg - pg) %>%
    select(dist_param, tree_id, core_id, ring_id, year, incr_mm, age, dbh_mm = dbh_growth, ai, fg, pg) ->
    data$ring
    
  return(data)
}

releaseCalculate <- function(data = NULL,  gapAge = c(5:15), nprol = 7){
  #' @description function calculate the releases for individual trees
  #' @param data a list of three dataframes, output of growthCalculate function
  #' @param nprol number of years to consider that release is sustaind
  #' @param gapAge age of the tree when it shall be tested for gap origin
  
  # data quality check
  options(error = NULL) # not to enter debug mode
  
  # perform the checks
  if(!is.list(data)) stop('The input data is not a list of three tables')
  if(!identical(c('core',"dist_param","ring"), ls(data))) stop('The input data tables dont match with required')
  
  aith <- data$dist_param  %>% select(dist_param, ai_mm) %>% deframe()
  gapth <- data$dist_param  %>% select(dist_param, gap_mm) %>% deframe()
  
  # calculate releases
  data$ring %>%
    arrange(year) %>%
    group_by(core_id) %>%
    mutate(event = ifelse(row_number() %in%  peakDetection(x = ai, threshold = aith[first(as.character(dist_param))], nups = 1,  mindist = 30), 'release', NA),
           event = ifelse(lead(fg, nprol) <= pg, NA, event),
           event = ifelse(lag(pg, nprol) >= fg, NA, event)) %>%
    filter(!is.na(event)) %>%
    select(core_id, year, event) ->
    release.event
  
  # calculate the gap origin 
  data$ring %>% 
    filter(age %in% gapAge) %>%
    arrange(year) %>%
    group_by(core_id) %>%
    summarise(dist_param = first(dist_param),
              gapGrowth = mean(incr_mm, na.rm = T),
              N = n(),
              year = min(year)) %>%
    filter(N >= 5,
           gapGrowth >= gapth[as.character(dist_param)]) %>%
    mutate(event = 'gap') %>%
    select(core_id, year, event) ->
    gap.event
  
  # add those that don't have any event
  data$ring %>%
    filter(!core_id %in% c(unique(gap.event$core_id), unique(release.event$core_id))) %>%
    group_by(core_id) %>%
    summarise(year = min(year)) %>%
    mutate(event = 'no event') ->
    no.event
    
  # add together the events
  bind_rows(release.event, gap.event, no.event) %>%
    arrange(year) %>%
    group_by(core_id) %>%
    mutate(keeprel = keepRelease(year, event, n = 30)) %>%
    ungroup() %>%
    filter(keeprel %in% 'yes') %>%
    inner_join(., data$ring, by = c('core_id', 'year')) %>%
    select(ring_id, dist_param, year, age, dbh_mm, ai, event) ->
    data$event
    
  return(data)
  
}



mdsFun <- function(ca, k = 30, bw = 5, st = 7){
  #' @description return a vector of the fited KDE function
  #' @param ca arranged vector of the canopy area values
  #' @param k a windows length, default 30
  #' @param bw a smoothing bandwidth to be used, default = 5
  #' @param st a standartization value, to scale back to canopy area
  
  rollapply( ca, 
             width = k,
             FUN = function(x){n <- length(x); density(1:n, weights = x, bw = bw, n = n)$y[round((n+1)/2)]* 100/st},
             fill = 0,
             align = "center",
             partial = TRUE)
}

predictDisYear <- function(data){
  #' @description predict the year of the disturbance for the all plot
  
  data <- data %>% select(x_m, y_m, year)
  sp::coordinates(data) <- c('x_m', 'y_m')
  
  grd <- tibble(x_m = seq(-22,22,1), y_m = seq(-22,22,1)) %>% complete(x_m, y_m)
  sp::coordinates(grd) <- c('x_m', 'y_m')
  gridded(grd) <- TRUE
  
  idw <- gstat::idw(formula = year ~ 1, locations = data, newdata = grd, idp = 8)
  
  idw.output <- as.data.frame(idw) %>%
    select(x_m, y_m, year = var1.pred)
  
  return(idw.output)
}

predictDisStand <- function(data){
  #' @description predict the year of the disturbance for the all plot
  
  data <- data %>% select(x_m = lng, y_m = lat, year)
  sp::coordinates(data) <- c('x_m', 'y_m')
  
  grd <- cbind(x_m = seq(min(data$x_m),max(data$x_m),0.0001), y_m = seq(min(data$y_m),max(data$y_m),0.0001)) %>% as.tibble() %>% complete(x_m, y_m)
  sp::coordinates(grd) <- c('x_m', 'y_m')
  gridded(grd) <- TRUE
  
  idw <- gstat::idw(formula = year ~ 1, locations = data, newdata = grd, idp = 8)
  
  idw.output <- as.data.frame(idw) %>%
    select(x_m, y_m, year = var1.pred)
  
  return(idw.output)
}


circleFun <- function( r = 12.5 ){
  #' @description Function to create the sircle.
  #' @return a data frame with x and y coordinate
  #' @param r A radius of a sircle
  
  tt <- seq(0,2*pi,length.out = 100)
  xx <- r * cos(tt)
  yy <- r * sin(tt)
  return(data.frame(X = xx, Y = yy))
}

base_size <- 7
gstyle <- list(
  theme_bw(base_size = base_size),
  theme(axis.text = element_text(size = base_size *1.2, angle = 0,  colour = "grey10"),
        axis.title=element_text(size=base_size*1.4),
        axis.ticks.length=unit(base_size*-0.15, "mm"),
        axis.ticks = element_line(size = base_size * 0.05),
        axis.text.x = element_text(margin=margin(2,0,1,0,"mm")),
        axis.text.y = element_text(margin=margin(0,2,0,1,"mm"))) ,
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black")) ,
  theme(strip.background = element_rect(colour = "white", fill = "white", size = base_size*1.1),
        strip.text.x = element_text(colour = "black", angle = 0, size = base_size*1.1,
                                    hjust = 0.5, vjust = 0.5)),
  theme(legend.key.size =  unit(base_size*1, "mm"), legend.key.height =  unit(base_size*1, "mm")),
  theme(legend.text = element_text(size=base_size * 1.2),
        legend.title = element_text(size=base_size * 1.5)),
  theme(legend.key = element_rect(colour = 'white', fill = 'white', linetype='dashed', size =0.1)),
  theme(legend.justification=c(1,1), legend.position=c(1,1)),
  scale_color_manual(values = c("Abies" = "#7CAE00",
                                'Acer' = 'orange',
                                "Fagus" = "#da2c3a",
                                "Others" = "grey50",
                                "Picea" = "#78c2ef",
                                'Pinus' = 'blue',
                                'NA' = "grey50",
                                "gap" = "#7CAE00",
                                "release" = "#da2c3a",
                                "no event" = "grey50",
                                'alive' = "#7CAE00",
                                'dead' = "#da2c3a",
                                'stump' =  "grey50",
                                "10_20_10" = "#7CAE00",
                                #"10_10_3" = "#da2c3a",
                                "10_10_5" = "#78c2ef"), drop = F),
  scale_shape_manual("Status",values = c("dead" = 21,
                                         "alive" = 19,
                                         "stump" = 4), drop = F),
  scale_alpha_manual("Used for disturbance",values = c('yes' = 1, 'no' = 0.3), guide = F))

gfill <- list(  scale_fill_manual(values = c("Abies" = "#7CAE00",
                               'Acer' = 'orange',
                               "Fagus" = "#da2c3a",
                               "Others" = "grey50",
                               "Picea" = "#78c2ef",
                               'Pinus' = 'blue',
                               'NA' = "grey50",
                               "gap" = "#7CAE00",
                               "release" = "#da2c3a",
                               "no event" = "grey50",
                               'alive' = "#7CAE00",
                               'dead' = "#da2c3a",
                               'stump' =  "grey50"), drop = F))

plotDistPos <- function(data, name){
  #data <- data.all.nest[80,] %>% unnest()
  # name <- 'test'
  
  if(!is.null(data[['dist.plot']][[1]])){
    p <- data[['position']][[1]] %>%
      ggplot( aes(x_m, y_m))+
      geom_tile(data = data[['dist.plot']][[1]], aes(x_m, y_m, fill = year), alpha = 0.4) +
      geom_point( aes(size = dbh_mm, color = Species, shape = status)) +
      geom_text_repel(aes(size = dbh_mm,label = year, color = species, alpha = dist_use), size  = 3)
    
  }else{
    p <- data[['position']][[1]] %>%
      ggplot( aes(x_m, y_m))+
      geom_point( aes(size = dbh_mm, color = Species, shape = status)) +
      geom_text_repel(aes(size = dbh_mm,label = year, color = species, alpha = dist_use), size  = 3)
  }
  
  p +
    scale_fill_gradient2(low="#78c2ef",  mid = "#7CAE00", high = "#da2c3a", guide = F, midpoint = 1850, limits = c(1650, 2010)) +
    scale_size_continuous("DBH (mm)", 
                          limits = c(0,1400),
                          breaks=c(0,200, 400, 600, 1400),
                          range = c(2,7), guide = F) +
    gstyle +
    geom_point( aes(0,0), shape = 3, color = "red",size = 3) +
    geom_path(data = circleFun(r = 12.5), aes(x = X, y = Y), color = "black", size = 0.3)+
    geom_path(data = circleFun(r = 17.84), aes(x = X, y = Y), color = "black", size = 0.3)+
    geom_path(data = circleFun(r = 21.85), aes(x = X, y = Y), color = "black", size = 0.3) +
    theme(legend.position='right') +
    ggtitle(name) ->
    position_gg
  
  # dbh
  data[['age_dbh']][[1]] %>%
    filter(status == 'alive') %>%
    ggplot() +
    geom_histogram(aes(dbh_mm, fill = Species, alpha = dist_use), binwidth = 50) +
    coord_cartesian(xlim = c(0, 800)) +
    gstyle + xlab('DBH (mm)') + ylab('Count') +
    theme(legend.position = "none") +
    gfill ->
    dbh_gg
  
  # age
  data[['age_dbh']][[1]] %>%
    ggplot() +
    geom_histogram(aes(year_min, fill = Species, alpha = dist_use), binwidth = 5) +
    coord_cartesian(xlim = c(1700, 2000)) +
    gstyle + xlab('Calendar year') + ylab('Count') +
    theme(legend.position = "none") +
    gfill ->
    age_gg
  
  # disturbacne history plot
  data[['dist']][[1]] %>%
    ggplot()+
    geom_bar(data = data[['data.dist.curr']][[1]], aes(year, ca), stat = 'identity', colour = 'orange') +
    geom_histogram(aes(year, weight = ca), binwidth = 10, fill = 'grey80', breaks = seq(1700, 2010, 10)) +
    geom_histogram(aes(year, weight = ca, fill = Species), binwidth = 1) +
    geom_line(data = data[['mds']][[1]], aes(year, value), size = 1, colour = 'grey20') +
    geom_point(data = data[['peaks']][[1]], aes(year, value, colour = method), size = 3) +
    geom_text_repel(data = data[['peaks']][[1]], aes(year, value, label = year, colour = method), size = 3) +
    geom_hline(aes(yintercept = 10), linetype = 2, colour = 'grey80') +
    geom_hline(aes(yintercept = 15), linetype = 2, colour = 'grey80') +
    facet_wrap(~plot_type, ncol = 1) +
    gstyle + 
    #theme(legend.position = "none") +
    theme(legend.justification = c(0, 1), legend.position = c(0, 1)) + guides(fill=FALSE) + 
    coord_cartesian(xlim = c(1700, 2020)) + coord_cartesian(ylim = c(0, 60)) +
    geom_point(aes(2020, 60), color = 'white') +
    ylab('Canopy area (%)') +xlab('Year') +
    gfill ->
    dist_gg
  
  # growth trend
  data[['growth.trend']][[1]] %>%
    ggplot()+
    geom_line(aes(year, incr_mm, colour = species, alpha = dist_use), size = 1) +
    geom_point(aes(1700, 3), color = 'white') +
    geom_point(aes(2020, 3), color = 'white') +
    gstyle + 
    theme(legend.position = "none") +
    coord_cartesian(xlim = c(1700, 2020)) + coord_cartesian(ylim = c(0, 3)) +
    ylab('TRW (mm)') +xlab('Year') +
    annotation_custom(grob = tableGrob(data[['data.dist.info']][[1]] %>% t(),theme = ttheme_minimal(base_size = 8)), xmax=1760, ymax=4) ->
    growth_gg
  
  # Join all figures
  age_dbh_gg <- plot_grid(dbh_gg, age_gg, nrow = 1, ncol = 2)
  tree_plot_gg <- plot_grid(position_gg, age_dbh_gg, nrow = 2, ncol = 1, rel_heights = c(4,2))
  dist_mds_gg <- plot_grid(dist_gg, growth_gg, nrow = 2, ncol = 1, rel_heights = c(2, 1.05))
  
  fin_plot_gg <- plot_grid(dist_mds_gg, tree_plot_gg, nrow = 1, ncol = 2, rel_widths = c(1.2, 1))
  
  cat(name)
  fin_plot_gg
}

plotDistStand <- function(data, name = 'test'){
  # data <- data.all.nest[5,] %>% unnest()
  
  data[['mds']][[1]] %>%
    ggplot()+
    geom_histogram(aes(year, weight = ca), fill = 'grey80', breaks = seq(1700, 2010, 10)) +
    geom_histogram(aes(year, weight = ca), fill = 'grey40', breaks = seq(1700, 2010, 1)) +
    geom_line(aes(year, kde), size = 1, colour = 'grey20') -> 
    p
    if(!is.null(data[['peaks']][[1]])){
      p + 
        geom_point(data = data[['peaks']][[1]], aes(year, kde), size = 3, colour = 'red') +
        geom_text_repel(data = data[['peaks']][[1]], aes(year, kde, label = year), size = 3) ->
        p
    }  
  p +
    geom_hline(aes(yintercept = 10), linetype = 2, colour = 'grey80') +
    gstyle + 
    coord_cartesian(xlim = c(1700, 2000), ylim = c(0, 40)) +
    ylab('Canopy area (%)') +xlab('Year') +
    ggtitle(name) ->
    dist_gg
  
  data[['mds.plot']][[1]] %>%
    ggplot() +
    geom_line(aes(year, kde, group = plot_id), size = 1, colour = 'grey70') +
    geom_line(data = data[['mds']][[1]], aes(year, kde), size = 1, colour = 'grey20') +
    geom_hline(aes(yintercept = 10), linetype = 2, colour = 'grey80') +
    gstyle + 
    coord_cartesian(xlim = c(1700, 2000), ylim = c(0, 80)) +
    ylab('Canopy area (%)') +xlab('Year') ->
    dist_plot_gg
  
  
  # peaks based on the number of plots
  data[['plot_n_mds']][[1]] %>%
    ggplot()+
    geom_histogram(aes(year, weight = ca), fill = 'grey80', breaks = seq(1700, 2010, 10)) +
    geom_histogram(aes(year, weight = ca), fill = 'grey40', breaks = seq(1700, 2010, 1)) +
    geom_line(aes(year, kde), size = 1, colour = 'grey20') -> 
    p
  if(!is.null(data[['plot_n_mds.peaks']][[1]])){
    p + 
      geom_point(data = data[['plot_n_mds.peaks']][[1]], aes(year, kde), size = 3, colour = 'red') +
      geom_text_repel(data = data[['plot_n_mds.peaks']][[1]], aes(year, kde, label = year), size = 3) ->
      p
  }  
  p +
    geom_hline(aes(yintercept = 10), linetype = 2, colour = 'grey80') +
    gstyle + 
    coord_cartesian(xlim = c(1700, 2000), ylim = c(0, 40)) +
    ylab('Canopy area (%)') +xlab('Year') +
    ggtitle("Based on number of plot") ->
    dist_n_gg
  
  
  # peaks based on the percent from plot
  data[['plot_per_mds']][[1]] %>%
    ggplot()+
    geom_histogram(aes(year, weight = ca), fill = 'grey80', breaks = seq(1700, 2010, 10)) +
    geom_histogram(aes(year, weight = ca), fill = 'grey40', breaks = seq(1700, 2010, 1)) +
    geom_line(aes(year, kde), size = 1, colour = 'grey20') -> 
    p
  if(!is.null(data[['plot_per_mds.peaks']][[1]])){
    p + 
      geom_point(data = data[['plot_per_mds.peaks']][[1]], aes(year, kde), size = 3, colour = 'red') +
      geom_text_repel(data = data[['plot_per_mds.peaks']][[1]], aes(year, kde, label = year), size = 3) ->
      p
  }  
  p +
    geom_hline(aes(yintercept = 10), linetype = 2, colour = 'grey80') +
    gstyle + 
    coord_cartesian(xlim = c(1700, 2000), ylim = c(0, 40)) +
    ylab('Canopy area (%)') +xlab('Year') +
    ggtitle("Based on percentage of plot peaks") ->
    dist_per_gg
  
  
  # data[['peaks.plot']][[1]] %>%
  #   ggplot()+
  #   geom_histogram(aes(year), fill = 'grey50' , breaks = seq(1700, 2010, 5)) +
  #   gstyle + 
  #   coord_cartesian(xlim = c(1700, 2000), ylim = c(0, 10)) +
  #   ylab('Number of plots') +xlab('Year') ->
  #   event_plot_gg
  
  # prediction
  if(!is.null(data[['pred.dist']][[1]])){
    p <- ggplot() +
      geom_tile(data = data[['pred.dist']][[1]], aes(x_m, y_m, fill = year), alpha = 0.4)
  }else{
    p <- ggplot()
  }
  
  p + 
    geom_point(data = data[['peaks.plot']][[1]], aes(lng, lat, size = kde)) +
    geom_text_repel(data = data[['peaks.plot']][[1]], aes(lng, lat, label = year, size = kde/2)) +
    scale_fill_gradient2(low="#78c2ef",  mid = "#7CAE00", high = "#da2c3a", guide = F, midpoint = 1850, limits = c(1650, 2010)) +
    gstyle +
    theme(legend.position = 'none') ->
    dist.spat
  
  up.pan <- plot_grid(dist_gg, dist_n_gg, dist_per_gg, ncol = 3, nrow = 1)
  low.pan <- plot_grid(dist_plot_gg, dist.spat, nrow = 1, ncol = 2, rel_widths = c(3,2))
  fin_plot_gg <- plot_grid(up.pan, low.pan, nrow = 2, ncol = 1)
  
  cat(name)
  fin_plot_gg
}
