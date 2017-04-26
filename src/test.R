

# Make models to predict well water elevations from rainfall

library(tidyverse)
library(h2o) # 3.10.4.2

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

site <- site[3,]

i <- 0

###############
# Create models
###############

h2o.init()

while(i < nrow(site)) {
  
  i <- i+1;

  elevations <- read_csv(site$w[i]) %>%
    mutate(dt = as.POSIXct(Date, format="%m/%d/%y %H:%M:%S"),
           hr = format(dt, "%H")) %>%
    filter(hr == "00") %>%
    arrange(dt) %>%
    select(dt, level = Water_Level_ft)
  
  rain <- read_csv(site$p[i]) %>%
    mutate(dt = as.POSIXct(TimeStamp, format="%m/%d/%y %H:%M:%S", tx = "PST")) %>%
    select(dt, rain = `Rain (inches)`)
  
  source("./src/predict_with_h2o.R")
  
  precip <- fillPrecip(rain$rain, rain$dt, "hour")
  
  ##########################################
  # Create weekly breaks going back one year,
  #   and plot the variable responses
  ##########################################
  
  breaks <- c(1, seq(from = 1, to = 52, by = 1)*24*7)
  windows <- makeHourlyWindows(breaks)
  
  # Changes in how the data are divided between the test and training
  #   sets can change variable importance, so run multiple times &
  #   calculate mean and range
  
  nruns <- 0
  imp <- data.frame()
  
  while (nruns <= 8) {
    nruns <- nruns+1
    
    importance <- defineHydrology(precip$precip, precip$date.time, 
                                  elevations$level, elevations$dt,
                                  windows)
    imp <- rbind(imp, importance[[1]])
    
  }

  imp.all <- imp %>%
    group_by(title) %>%
    summarize(weeks = mean(weeks),
              imp.mean = mean(percentage),
              imp.max = max(percentage),
              imp.min = min(percentage))
  
  ggplot(imp.all, aes(x = weeks,
                      y = imp.mean)) +
    geom_ribbon(aes(x = weeks,
                    ymin = imp.min,
                    ymax = imp.max),
                fill = "#cccccc") +
    geom_line() +
    geom_point() +

    labs(title = paste0(site$n[i], " Importance of Weekly Total Precipitation"),
         x = "# Weeks Preceeding Measurement",
         y = "Importance (0-1)") +
    theme_minimal()
  
  ggsave(paste0("./results/", site$n[i], "_var_imp.png"),
         width = 7,
         height = 4)
  
  ##########################################
  # Train a better model, by combining weeks
  #   And plot the variable response
  ##########################################
  

  model.refined <- buildModel(precip$precip, precip$date.time,
                              elevations$level, elevations$dt)
  breaks.new <- model.refined[[1]]
  m <- model.refined[[2]]
  b <- model.refined[[3]]

  windows <- makeHourlyWindows(breaks.new)

  importance <- defineHydrology(precip$precip, precip$date.time,
                                elevations$level, elevations$dt,
                                windows)

  # Plot final variable importance
  ggplot(importance[[1]], aes(x = weeks,
                              y = percentage)) +
    geom_line() +
    geom_point() +
    labs(title = paste0(site$n[i], " Importance of Weekly Total Precipitation"),
         subtitle = "Refined model",
         x = "# Weeks Preceeding Measurement",
         y = "Importance (0-1)") +
    theme_minimal()

  ggsave(paste0("./results/", site$n[i], "_var_imp_refined.png"),
         width = 7,
         height = 4)

  # Plot change in r2 with consolidation of breakpoints
  ggplot(m, aes(x = nvar,
                y = r2)) +
    geom_line() +
    geom_point() +
    scale_x_reverse() +
    labs(title = paste0(site$n[i], " Change in R2"),
         subtite = "Change in model fit as rainfall windows are condensed",
         x = "Number of Variables (or windows)",
         y = "R2") +
    theme_minimal()

  ggsave(paste0("./results/", site$n[i], "_r2_change.png"),
         width = 7,
         height = 4)

  # Plot change in sd with consolidation of breakpoints
  ggplot(m, aes(x = nvar,
                y = sd)) +
    geom_line() +
    geom_point() +
    scale_x_reverse() +
    labs(title = paste0(site$n[i], " Change in Standard Deviation"),
         subtite = "Change in variable uniformity as rainfall windows are condensed",
         x = "Number of Variables (or windows)",
         y = "Standard Deviation") +
    theme_minimal()

  ggsave(paste0("./results/", site$n[i], "_sd_change.png"),
         width = 7,
         height = 4)
}


