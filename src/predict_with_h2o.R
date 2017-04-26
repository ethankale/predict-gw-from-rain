
# Functions for predicting water elevations from rainfall
# Uses random forest function from H2O.ai.

##################
# Helper functions
##################

# Calculate rolling sums on slices of data
# Regular R/zoo rollapply is too slow.
rollSum <- function(series, begin, end) {
  library(zoo)
  library(RcppRoll)
  
  width = begin - end
  precip.rollsum <- roll_sumr(series,
                              n = width)
  precip.rollsum <- precip.rollsum[1:(length(precip.rollsum)-begin)]
  precip.rollsum <- c(rep(NA, times = begin), precip.rollsum)
  return(precip.rollsum)
  # return(rollapplyr(series, 
  #                   width = list(-begin: -end), 
  #                   FUN = sum,
  #                   fill = NA))
}

# Build a precip series with cumulataive windows over various time periods
precipCumSums <- function(precip, windows) {
  # This should probably be done with apply, but I have a hell
  #   of a time getting apply to work properly, plus it's a 
  #   pretty short loop over a data.frame - not much time benefit
  
  message("Generating cumulative rainfall")
  
  i <- 0
  df <- data.frame(pr = precip);
  
  pb <- txtProgressBar(min = 0, max = nrow(windows), initial = 0)
  
  while (i < nrow(windows)) {
    i <- i+1
    df[as.character(windows[i,1])] <- rollSum(precip, windows[i,2], windows[i,3]);
    setTxtProgressBar(pb, i)
  }
  
  close(pb);
  return(df);
}


#####################
# Hydrology functions
#####################

# Fill precip data
#   precip and dates must be in the same, sequential order and of same length
#   Returns a tibble
#   Interval must work with seq()
fillPrecip <- function(precip, date.time, interval) {
  
  library(tidyverse)
  library(zoo)
  
  joined <- tibble(precip = precip, date.time = date.time)
  
  message(paste0(sum(is.na(joined)), " NA values filled with 0."))
  
  date.time.full <- seq(from = min(date.time, na.rm = TRUE), 
                     to = max(date.time, na.rm = TRUE), 
                     by = interval)
  
  
  joined.full <- tibble(date.time = date.time.full) %>%
    left_join(joined)
  
  joined.full$precip <- na.fill(joined.full$precip, 0)
  
  return(joined.full)
}

makeHourlyWindows <- function(breaks) {

  begin <- breaks[2:length(breaks)]
  end <- breaks[1:length(breaks)-1]
  title <- paste0("b", seq(from = 1, to = length(breaks)-1, by = 1))
  
  windows <- data.frame(title, begin, end)
  windows$title <- as.character(title)
  
  return(windows)
}

# Characterize the relationship between rainfall and elevations.
#   precip is a list/vector of precipitation data
#   precip.datetime is the list/vector of dates that correspond with 
#     values in precip
#   elevation is a list/vector of elevations (well, lake)
#   elevation.datetime is the list/vector of dates that correspond with
#     values in elevation
#
#   dates & times must match exactly between precip & elevation; suggest using
#     either hourly or 15 minute data
defineHydrology <- function(precip, precip.datetime, 
                            elevation, elevation.datetime,
                            windows) {
  library(h2o)
  library(zoo)
  
  # Create a series of windows
  #   Consider making the extent configurable
  # weeks <- seq(7, 358, by = 7)
  # begin <- c(1, weeks*24)
  # end <- c(weeks, 365)*24
  # title <- paste0("w", seq(1, length(begin), by = 1))
  # 
  # windows <- data.frame(title, begin, end)
  # windows$title <- as.character(title)
  
  # Create columns representing cumulative rainfall over multiple weeks
  precip.join <- cbind(tibble(dt = precip.datetime, p = precip), 
                       precipCumSums(precip, windows))
  elev <- tibble(dt = elevation.datetime, elev = elevation)
  
  e.and.p <- elev %>%
    left_join(precip.join)
  
  # Import the data into H2O and build the model
  #h2o.init()
  
  temp.path <- "./data/tmp.csv"
  write_csv(e.and.p, temp.path)
  e.and.p.hex <- h2o.importFile(temp.path, "e.and.p.hex")
  
  # Split for training, then make the model
  e.and.p.split <- h2o.splitFrame(e.and.p.hex, ratios = c(0.8))
  e.and.p.drf.level <- h2o.randomForest(x = windows$title,
                                        y = "elev",
                                        training_frame = e.and.p.split[[1]],
                                        model_id = "e.and.p.drf.elev",
                                        validation_frame = e.and.p.split[[2]])
  
  # Pull out and return the variable importance
  variable.import <- as.data.frame(h2o.varimp(e.and.p.drf.level)) %>%
    mutate(title = variable) %>%
    left_join(windows) %>%
    mutate(weeks = (((begin+end)/2) / (24*7))) %>%
    arrange(weeks)
  
  # Pull out and return predicted elevations
  e.and.p.drf.predict <- as.data.frame(h2o.predict(e.and.p.drf.level,
                                     e.and.p.hex))
  
  # Clean up
  file.remove(temp.path)
  #h2o.shutdown()
  
  return(list(variable.import, 
              h2o.r2(e.and.p.drf.level),
              e.and.p.drf.predict))
}

# Build a model by automatically collapsing neighboring variables (weeks)
#   with low importance scores together

buildModel <- function(precip, precip.datetime, 
                       elevation, elevation.datetime) {
  library(moments)  

  # Look at how well the model performs as variables are grouped
  metrics.iter <- c()
  metrics.nvar <- c()
  metrics.sd <- c()
  metrics.r2 <- c()
  metrics.skew <- c()
  
  # Track breakpoints
  b <- list()
  
  breaks <- c(1, seq(from = 1, to = 52, by = 1)*24*7)
  windows <- makeHourlyWindows(breaks)
  n <- nrow(windows)
  
  # We return the breaks with the best R2.  Somewhat naive metric.
  breaks.best <- breaks
  i <- 0
  end <- 0
  
  while(i < n) {
    
    i <- i+1
    importance <- defineHydrology(precip, precip.datetime, 
                                  elevation, elevation.datetime,
                                  windows)
    imp.percent <- importance[[1]]$percentage
        
    r2.new <- importance[[2]]

    sd.new <- sd(imp.percent)

    message("R Squared: ", r2.new)
    message("Standard Deviation: ", sd.new)
    message("Number of Variables: ", nrow(windows))
    
    # Update metrics for this iteration
    metrics.nvar <- c(metrics.nvar, nrow(windows))
    metrics.iter <- c(metrics.iter, i)
    metrics.sd <- c(metrics.sd, sd.new)
    metrics.r2 <- c(metrics.r2, r2.new)
    metrics.skew <- c(metrics.skew, skewness(imp.percent))
    
    # We want to evaluate not just the importance of one specific variable,
    #   but the importance of the variables around it.  This is important
    #   because lumping or splitting effects neighboring variables as well.
    vars.imp.mean <- rollmean(imp.percent,
                              3,
                              align = "center",
                              fill = median(imp.percent))
    
    # Switch between lumping and splitting variables, depending on skewness.
    #   Distribution of importance starts out skewed right; we bring
    #   it down by lumping variables.  If it drops too far we can split
    #   variables.  This will tend to reduce the number of variables to a
    #   stable minimum and drive up R2 (hopefully).
    if(skewness(imp.percent) >= 1) {
      message("Lumping...")
      # Select variables to remove, w/ lowest explanatory power
      vars.lowest <- which.min(vars.imp.mean)
      breaks <- breaks[-vars.lowest]
      windows <- makeHourlyWindows(breaks)

    } else {
      message("Splitting...")
      vars.highest <- which.max(vars.imp.mean)
      b2 <- breaks[c(vars.highest-1, vars.highest+1)]
      b3 <- seq(from = b2[1], to = b2[2], by = ((b2[2] - b2[1])/3))[2:3]
      
      b <- c(breaks[c(1:(vars.highest-1))],
             b3,
             breaks[c((vars.highest+1):length(breaks))])
      breaks <- round(b)
      windows <- makeHourlyWindows(breaks)

    }
    
    #message(paste0(breaks, ":"))
    metrics <- data.frame(iter = metrics.iter,
                          nvar = metrics.nvar,
                          sd = metrics.sd,
                          r2 = metrics.r2,
                          skew = metrics.skew)
    
    # Return the breaks with the best R2
    if (r2.new >= max(metrics.r2)) {
      breaks.best <- breaks
    }
    
  }
  return(list(breaks.best, metrics))
}



