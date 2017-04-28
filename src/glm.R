
# Predict well elevations with GLM


library(tidyverse)
library(relaimpo)

#############
# Import data
#############

site.name <- c("LRS-1", "LRS-7", "LRS-8", "LRS-9",
               "LRS-11", "LRS-12", "MW2", "MW3", "MW5")
site.well <- paste0("./data/",
                    c("LRS1a_17","LRS7a_17", "LRS8","LRS9_17",
                      "LRS11a_17", "LRS12_17", "MW2_17", "MW3_17", "MW5_17"),
                    "_12hr.csv")
site.precip <- paste0("./data/",
                      c("11u17","11u17","11u17","11u17","11u17",
                        "11u17","11u17","11u17", "11u17"),
                      "_1hr.csv")
# site.name <- c("Test")
# site.well <- c("./data/test.csv")
# site.precip <- c("./data/rain_test.csv")

site <- tibble(n = site.name, 
               w = site.well, 
               p = site.precip)

# site <- site[3,]

i <- 0

###############
# Create models
###############

while(i < nrow(site)) {
  
  i <- i+1;
  dir.current <- paste0("./results/", site$n[i], "/")
  dir.create(dir.current, showWarnings = FALSE)
  
  elevations <- read_csv(site$w[i]) %>%
    mutate(dt = as.POSIXct(Date, format="%m/%d/%y %H:%M:%S",
                           tz = "PST"),
           hr = format(dt, "%H")) %>%
    filter(hr == "00") %>%
    arrange(dt) %>%
    dplyr::select(dt, level = Water_Level_ft)
  
  rain <- read_csv(site$p[i]) %>%
    mutate(dt = as.POSIXct(TimeStamp, format="%m/%d/%y %H:%M:%S")) %>%
    dplyr::select(dt, rain = `Rain (inches)`)
  
  source("./src/predict_with_h2o.R")
  
  precip <- fillPrecip(rain$rain, rain$dt, "hour")
  
  # Fill the gaps in the elevation record with NAs.
  #   Converting to date and back to POSIX messes up the time, so fix that
  #   with some math.  Specifying the correct tz might work too.
  elev.filled <- fillElevation(elevations$level, 
                               as.Date(elevations$dt), "day") %>%
    mutate(dt = as.POSIXct(date.time)+(8*60*60),
           level = elev) %>%
    dplyr::select(dt, level)
  
  ##########################################
  # Create weekly breaks going back one year,
  #   and plot the variable responses
  ##########################################
  
  breaks <- c(1, seq(from = 1, to = 52, by = 1)*24*7)
  windows <- makeHourlyWindows(breaks)
  
  # Create columns representing cumulative rainfall over multiple weeks
  precip.join <- cbind(precip, precipCumSums(precip$precip, windows)) %>%
    mutate(dt = date.time)

  e.and.p.orig <- elev.filled %>%
    left_join(precip.join) %>%
    mutate(e.01 = lag(level, 1),
           e.02 = lag(level, 2),
           e.03 = lag(level, 3),
           e.04 = lag(level, 4),
           e.05 = lag(level, 5),
           e.06 = lag(level, 6),
           e.07 = lag(level, 7),
           e.14 = lag(level, 14),
           e.30 = lag(level, 30))
  
  if (sum(complete.cases(e.and.p.orig)) <= 30) {
    message(paste0("Skipping ", site$n[i]))
    next
  }
  
  predict.window <- tibble(name = c("1 Day", "2 Days", "3 Days", "4 Days", 
                                    "5 Days", "6 Days", "7 Days", "14 Days",
                                    "30 Days"),
                           var.name = c("e.01", "e.02", "e.03", "e.04", "e.05",
                                         "e.06", "e.07", "e.14", "e.30"),
                           span = c(1,2,3,4,5,6,7,14,30))
  
  write_csv(e.and.p.orig,
            paste0(dir.current, site$n[i], "_complete.csv"))
  
  # Set aside a section of data to validate results
  split.date <- as.POSIXct("2015-10-01 00:00:00")
  e.and.p <- e.and.p.orig[which(e.and.p.orig$dt <= split.date),]
  
  # Create a generalized linear model for each prediction window
  j <- 0
  m.p95 <- c()
  m.aic <- c()
  m.r2 <- c()
  while(j < nrow(predict.window)) {
    j <- j+1
    var.current <- predict.window$var.name[j]

    fm <- as.formula(paste0("level ~ ",
                           paste(windows$title, collapse = " + "),
                           " +",  var.current))

    e.and.p.lm <- lm(fm, data = e.and.p)
    
    # Predict elevations
    p <- predict(e.and.p.lm, newdata = e.and.p.orig)
    e.and.p.orig$p <- p
    
    # Calculate the overall accuracy of the predictions
    e.and.p.score <- e.and.p.orig[which(e.and.p.orig$dt > split.date),]
    diff <- e.and.p.score$level - e.and.p.score$p
    diff.p95 <- mean(abs(quantile(diff, .95, na.rm = TRUE)),
                     abs(quantile(diff, .05, na.rm = TRUE)))
    r2 <- summary(e.and.p.lm)$adj.r.squared
    
    m.p95 <- c(m.p95, diff.p95)
    m.aic <- c(m.aic, AIC(e.and.p.lm))
    m.r2 <- c(m.r2, r2)
    
    # Plot predictions
    ggplot(e.and.p.orig, aes(x = dt,
                             y = level)) +
      geom_line(aes(x = dt,
                     y = p,
                     color = "Predicted Value")) +
      geom_line() +
      labs(title = paste0(site$n[i], " Predicted and Actual Elevation"),
           subtitle = paste0(predict.window$name[j], 
                             "; 95% of predicted values within +/- ",
                             round(diff.p95, 3),
                             "; AIC: ",
                             round(AIC(e.and.p.lm), 1),
                             "; R2: ",
                             round(r2, 3)),
           x = "Date",
           y = "Elevation (feet)") +
      theme_minimal()
    
    ggsave(paste0(dir.current, site$n[i], 
                  "_", predict.window$var.name[j],
                  "_predicted.png"),
           width = 7,
           height = 4)
  }
  
  # Graph the change in metrics 
  metrics <- tibble(`95% Confidence Interval` = m.p95,
                    `AIC` = m.aic,
                    `R Squared` = m.r2,
                    span = predict.window$span)
  
  m <- metrics %>%
    gather(metric, value, -span)
  
  
  ggplot(m, aes(x = span,
                y = value)) +
    facet_wrap(~ metric,
               ncol = 1,
               scales = "free_y") +
    geom_line(color = "gray") +
    geom_point() +
    labs(paste0(site$n[i], " Fit Metrics"),
         x = "Prediction Span (days ahead)",
         y = "Value") +
    theme_minimal()
  
  ggsave(paste0(dir.current, site$n[i], 
                "_metrics.png"),
         width = 7,
         height = 4)
}



