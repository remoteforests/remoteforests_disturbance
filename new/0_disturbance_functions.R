distGetData <- function(ID){
  #' @description download data from the database based on the specified core id
  #' @param ID vector of core id to be used for the calcuation
  
  core.tbl <- tbl(KELuser, "core") %>%
    filter(id %in% ID) %>%
    inner_join(., tbl(KELuser, "tree") %>% select(tree_id = id, dbh_mm, species), by = "tree_id") %>%
    inner_join(., tbl(KELuser, "species_fk") %>% select(species = id, sp_group_dist), by = "species") %>%
    select(sp_group_dist, tree_id, core_id = id, missing_mm, missing_years, dbh_mm) %>%
    collect()
  
  ring.tbl <- tbl(KELuser, "ring") %>%
    filter(core_id %in% ID) %>%
    select(core_id, ring_id = id, year, incr_mm) %>%
    collect() %>%
    arrange(core_id, year) %>%
    group_by(core_id) %>%
    mutate(incr_mm = if_else(incr_mm %in% 0, NA_real_, incr_mm),
           incr_mm = na.approx(incr_mm)) %>%
    ungroup()
  
  dist_param.tbl <- tbl(KELuser, "dp") %>% collect()
  
  data <- list(dist_param = dist_param.tbl, core = core.tbl, ring = ring.tbl)
  
  return(data)
}

priorGrowth <- function(x, windowLength = 10){
  #' @description for each year calculate average growth in the specified preceeding period (including the focal year)
  #' @param x vector of incr_mm values
  #' @param windowLength length of the running window in years
  rollapply( x, 
             width = windowLength,
             FUN = mean,
             fill = NA,
             align = "right",
             na.rm = T,
             partial = TRUE)
}

followGrowth <- function(x, windowLength = 10){
  #' @description for each year calculate average growth in the specified following period (excluding the focal year)
  #' @param x vector of incr_mm values
  #' @param windowLength length of the running window in years
  rollapply( lead(x, 1), 
             width = windowLength,
             FUN = mean,
             fill = NA,
             align = "left",
             na.rm = T,
             partial = TRUE)
}

growthCalculate <- function(data.in, windowLength = 10){
  #' @description calculate the growth change, age, and dbh for individual trees
  #' @param data.in list of 3 tables ('dist_param', 'core', 'ring'), output of 'distGetData' function
  #' @param windowLength length of the running window for ai (absolute increase) calculation (in years)
  
  # data quality check
  # options(error = NULL) # not to enter debug mode
  
  # perform data quality checks
  if(!is.list(data.in)) stop("The input data is not a list of 3 data tables.")
  if(!identical(c("core","dist_param","ring"), ls(data.in))) stop("The input data tables should be 'dist_param', 'core', and 'ring'.")
  
  # calculate the dbh, age, and growth change
  growth.tbl <- inner_join(data.in$core, data.in$ring, by = "core_id") %>%
    arrange(core_id, year) %>%
    group_by(core_id) %>%
    mutate(#missing_mm = ifelse(is.na(missing_mm), 0, missing_mm),
           #missing_years = ifelse(is.na(missing_years, 0, missing_years)),
           dbh_growth = ifelse(row_number() == 1, incr_mm + missing_mm, incr_mm),
           dbh_growth = cumsum(dbh_growth) * 2,
           dbh_mm = ifelse(is.na(dbh_mm), max(dbh_growth), dbh_mm),
           dbh_coef = max(dbh_mm) / max(dbh_growth),
           dbh_growth = round(dbh_growth * dbh_coef, 0),
           age = year - min(year) + missing_years + 1,
           pg = priorGrowth(incr_mm, windowLength = windowLength),
           fg = followGrowth(incr_mm, windowLength = windowLength),
           ai = fg - pg) %>%
    ungroup() %>%
    select(sp_group_dist, tree_id, core_id, ring_id, year, incr_mm, age, dbh_mm = dbh_growth, pg, fg, ai)
  
  dist_param.tbl <- data.in$dist_param
  
  data.out <- list(dist_param = dist_param.tbl, growth = growth.tbl)
  
  return(data.out)
}

peakDetection <- function(x, threshold, nups = 1, mindist = 30){
  #' @description identify the index of year when release event occurs
  #' @param x vector of ai (absolute increase) values
  #' @param threshold minimum ai (absolute increase) value in mm
  #' @param nups number of increasing steps before the peak
  #' @param mindist minimum distance between two consecutive peaks (in years)
  
  # x <- ifelse(is.na(x), -0.2, x)
  
  x <- ifelse(is.na(x), 0, x)
  
  x <- findpeaks(x, 
                 minpeakheight = threshold,
                 minpeakdistance = mindist,
                 nups = nups) 
  
  if(is.null(x)){
    NA
  } else {
    matrix(x, ncol = 4)[,2]
  }
}

keepRelease <- function(year, type, n = 30){
  #' @description calculate the proximity of gap origin and release events
  #' @param year vector of years of individual events
  #' @param type type of the event (release or gap)
  #' @param n mimimal distance between the gap and release events (years)
  
  keep <- rep("yes", length(year))
  
  if(any(type %in% "gap")){
    diffyear <- year - year[type %in% "gap"]
    keep[diffyear < n & type %in% "release"] <- "no"
  }
  
  keep
}

eventCalculate <- function(data.in, gapAge = c(5:14), nprol = 7){
  #' @description calculate the canopy accession events for individual trees
  #' @param data.in list of 2 tables ('dist_param', 'growth'), output of 'growthCalculate' function
  #' @param nprol number of years to consider the release as sustained
  #' @param gapAge period of age when the tree is tested for gap origin
  
  # data quality check
  # options(error = NULL) # not to enter debug mode
  
  # perform data quality checks
  if(!is.list(data.in)) stop("The input data is not a list of 2 data tables.")
  if(!identical(c("dist_param","growth"), ls(data.in))) stop("The input data tables should be 'dist_param' and 'growth'.")
  
  aith <- data.in$dist_param %>% select(sp_group_dist, ai_mm) %>% deframe()
  dbhth <- data.in$dist_param %>% select(sp_group_dist, dbh_mm) %>% deframe()
  gapth <- data.in$dist_param %>% select(sp_group_dist, gap_mm) %>% deframe()
  
  # calculate releases
  release.event <- data.in$growth %>%
    arrange(core_id, year) %>%
    group_by(core_id) %>%
    mutate(event = ifelse(row_number() %in% peakDetection(x = ai, threshold = aith[first(sp_group_dist)], nups = 1, mindist = 30), "release", NA),
           event = ifelse(lead(fg, nprol) <= pg, NA, event),
           event = ifelse(lag(pg, nprol) >= fg, NA, event),
           event = ifelse(dbh_mm >= dbhth[first(sp_group_dist)], NA, event)) %>%
    ungroup() %>%
    filter(!is.na(event)) %>%
    select(core_id, year, event)
  
  # calculate gap origin 
  gap.event <- data.in$growth %>% 
    filter(age %in% gapAge) %>%
    arrange(core_id, year) %>%
    group_by(core_id, sp_group_dist) %>%
    summarise(incr_mean = mean(incr_mm, na.rm = T),
              nyears = length(incr_mm[!is.na(incr_mm)]),
              year = min(year) - (min(age) - 1)) %>%
    ungroup() %>%
    filter(nyears %in% 10,
           incr_mean >= gapth[sp_group_dist]) %>%
    mutate(event = "gap") %>%
    select(core_id, year, event)
  
  # add trees without any event
  no.event <- data.in$growth %>%
    filter(!core_id %in% c(unique(gap.event$core_id), unique(release.event$core_id))) %>%
    group_by(core_id) %>%
    summarise(year = min(year) - (min(age) - 1)) %>%
    ungroup() %>%
    mutate(event = "no event")
  
  # merge the events together -> perform 'keepRelease' check
  data.out <- bind_rows(release.event, gap.event, no.event) %>%
    arrange(core_id, year) %>%
    group_by(core_id) %>%
    mutate(keeprel = keepRelease(year, type = event, n = 30)) %>%
    ungroup() %>%
    filter(keeprel %in% "yes") %>%
    inner_join(., data.in$growth %>% distinct(., sp_group_dist, core_id, tree_id), by = "core_id") %>%
    select(sp_group_dist, tree_id, year, event) 
  
  return(data.out)
}

mdsFun <- function(ca, k = 30, bw = 5, st = 7){
  #' @description return a vector of the fitted KDE function
  #' @param ca arranged vector of the canopy area values
  #' @param k window length
  #' @param bw smoothing bandwidth
  #' @param st standartization value to scale back to canopy area
  
  rollapply( ca, 
             width = k,
             FUN = function(x){n <- length(x); density(1:n, weights = x, bw = bw, n = n)$y[round((n+1)/2)]* 100/st},
             # fill = 0,
             align = "center",
             partial = TRUE)
}
