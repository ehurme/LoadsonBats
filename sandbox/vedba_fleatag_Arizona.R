library(pacman)

p_load(data.table, seewave, tuneR, ggplot2)

path25 <- "../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230625/"
folders25 <- list.dirs(path25, full.names = TRUE)
path26 <- "../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230626/"
folders26 <- list.dirs(path26, full.names = TRUE)

folders <- c(folders25[2:4], folders26[2:4])
#folder <- folders[3:5]

files <- "C:/Users/ehumre/Dropbox/MPI/Wingbeat/Arizona/flightcage/20230610/Calibration/CalibrationTag1.csv"

sampling_rate <- 54/2
window_length = sampling_rate

i = 1
# for(i in 1:length(folders))
  files <- list.files(folders[i], pattern = ".csv", full.names = TRUE)
  df <- fread(files[1])
  # resample to match sigfox data
  idx <- seq(1, nrow(df), by = 2)
  df <- df[idx,]

  plot(df$accX_mg, ylim = c(-4000,4000))
  points(df$accY_mg, col = 2)
  points(df$accZ_mg, col = 3)
  # abline(v = 20)
  # df <- df[16520:16547,] # [16501:16527]
  sx <- zoo::rollmean(df$accX_mg, k = window_length, na.pad = TRUE) # mean(df$accX_mg) #
  sy <- zoo::rollmean(df$accY_mg, k = window_length, na.pad = TRUE)
  sz <- zoo::rollmean(df$accZ_mg, k = window_length, na.pad = TRUE)

  ax <- df$accX_mg - sx
  ay <- df$accY_mg - sy
  az <- df$accZ_mg - sz
  hist(az)

  ODBA <- abs(ax)+abs(ay)+abs(az)
  VeDBA <- sqrt(ax^2 + ay^2 + az^2)
  range(VeDBA, na.rm = T)
  plot(ODBA, VeDBA)

  av <- zoo::rollmean(VeDBA, k = sampling_rate, na.pad = TRUE)

  plot(az, ylim = c(0, 10000))
  points(ODBA, col = rgb(1,0,0,.1))
  points(VeDBA, col = rgb(0,1,0,.1))
  points(av, col = rgb(0,0,1,.1))
  legend("topright", col = c(1, 2, "green", "blue"), pch = 16, legend = c("Az", "ODBA", "VeDBA", "Rolling mean Vedba"))
  # sum(VeDBA)

  sv <- zoo::rollsum(VeDBA, k = sampling_rate, na.pad = TRUE)
  plot(sv/27, col = rgb(0,0,1,.1), ylab = "rolling sum of VeDBA")
  summary(sv)

