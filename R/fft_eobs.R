# AccelerateR loop to extract fft and params from acceleration data
# Edward Hurme
# 7/4/23

library(pacman)
p_load(data.table, janitor, lubridate, tidyverse, dplyr,
       accelerateR, move, # acc and gps analysis
       signal, tuneR, seewave # spectrogram analysis
       #lutz # look up timezones
       )
source("src/fft_peak_freq.R")
source("src/get_day_or_night.R")

path <- "./../../../ownCloud/Firetail/Acerodonjubatus/tag_1521/"
path <- "./../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/"
path <- "../../../ownCloud/Firetail/Eidolonhelvum/Model_tag_2396/"

sampling_rate = 18.74
highpass = 5
lowpass = 1
files <- list.files(path, pattern = "*.csv")
bats <- sapply(strsplit(files, split = "-"), '[', 1)
if(!dir.exists(paste0(path, "accelerateR"))){
  dir.create(paste0(path, "accelerateR"))
}
i=1
for(i in 1:length(files)){
  print(bats[i])
  df <- fread(paste0(path, files[i]))
  df$tag_local_identifier <- bats[i]
  dm <- {}
  df_clean <- df[!duplicated(df$timestamp),]

  # calculate speed from GPS points
  dm <- move(x = df_clean$location_long, y = df_clean$location_lat,
             time = df_clean$timestamp, data = df_clean,
             proj = "+proj=longlat +datum=WGS84 +no_defs")
  idx <- which(!is.na(df$location_lat))
  df$speed <- NA
  df$speed[idx] <- c(NA, speed(dm))
  plot(df$location_long, df$location_lat, asp = 1)
  points(median(df$location_long, na.rm = TRUE),
         median(df$location_lat, na.rm = TRUE), col = 2, pch = 16)
  plot(df$timestamp, df$speed)

  # Add a column to the data indicating day or night
  df$day_or_night <- mapply(get_day_or_night, mean(df$location_lat, na.rm = TRUE),
                            mean(df$location_long, na.rm = TRUE),
                            df$timestamp)
  df$behavior <- "resting"
  df$behavior[which(df$`annotation-layer-commuting` != "")] <- "commuting"
  df$behavior[which(df$`annotation-layer-foraging` != "")] <- "foraging"

  ACC <- {}
  ACC <- move_trans(data =  df[df$type == "acc",], timestamp = "timestamp", acc = "eobs_accelerations_raw",
                    id = "tag_local_identifier",
                    sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
                    naxes = 3, no_cores = 20)
  # check for any missing values
  burstcount = df$eobs_accelerations_raw[1] %>% strsplit(split = " ") %>% unlist %>% length/3
  count_test(ACC, burstcount)

  # add behavior, time, and burst to ACC

  ACC$daynight <- rep(df$day_or_night[df$type == "acc"],
                      each = burstcount)
  ACC$behavior <- rep(df$behavior[df$type == "acc"],
                      each = burstcount)
  ACC$burst <- rep(1:(nrow(ACC)/burstcount),
                      each = burstcount)

  ## fast fourier transform
  if(nrow(ACC[ACC$behavior == "commuting" |
              ACC$behavior == "foraging",]) > 0){
    bursts <- unique(ACC$burst)
    freqs <- data.table()

    j = 1
    for(j in 1:length(bursts)){
      b <- ACC[ACC$burst == bursts[j],]
      pc <- prcomp(x = b[,1:3])
      b$pc <- pc$x[,1]
      plot(b$pc, type = "l")

      duration <- nrow(b) / sampling_rate
        #difftime(max(b$timestamp), min(b$timestamp), units = "secs") %>% as.numeric()

      df <- fft_peak_freq(data = b$pc,
                          time = duration,
                          highpass = highpass,
                          lowpass = lowpass,
                          sampling_rate = sampling_rate)
      df$burst <- j
      df$burst_start <- b$timestamp[1]
      df$duration <- duration
      df$behavior <- b$behavior[1]
      freqs <- rbind(freqs, df)
    }

    # add flight numbers for long breaks in recordings
    burst_diff <- as.numeric(difftime(freqs$burst_start[-1],
                                      freqs$burst_start[-nrow(freqs)],
                                      units = "hours"))
    long_breaks <- which(abs(burst_diff) > 1) + 1

    day_counter <- 1
    for (i in 1:nrow(freqs)) {
      if (i %in% long_breaks) {
        day_counter <- day_counter + 1
      }
      freqs$day[i] <- day_counter
    }

    png(file = paste0(path, "/fft/freq_rms_", bats[i], ".png"),
        width = 800, height = 600)
    layout(rbind(1:2))
      with(freqs, plot(freq, rms, cex = amp/max(amp, na.rm = TRUE), col = behavior %>% factor))
      with(freqs, plot(wfreq, rms, cex = wamp/max(wamp, na.rm = TRUE), col = behavior %>% factor))
    dev.off()

    freqs$duration_of_burst <- ACC$burst_size[1]/ACC$sample_frequency[1]

    png(file = paste0(path, "/fft/wingbeat_", bats[i], ".png"),
        width = 800, height = 600)
      with(freqs, plot(time, frequency, cex = amp/max(amp, na.rm = TRUE),
                     pch = 1, main = bats[i], col = behavior %>% factor))
    dev.off()
  }

  # low pass filter
  sum_acc <- sum_data(ACC, time = "timestamp", x="x" ,
                      y="y" , z="z" , stats = "all",
                      behaviour = "behavior")
  sum_acc$speed <- NA
  freqs$speed <- NA
  for(j in 1:nrow(df)){
    idx1 <- which.min(abs(df$timestamp[j] - sum_acc$time))
    sum_acc$speed[idx1] <- df$speed[j]
    idx2 <- which.min(abs(df$timestamp[j] - freqs$time))
    freqs$speed[idx2] <- df$speed[j]
  }

  if(nrow(ACC[ACC$behavior == "commuting"|ACC$behavior == "foraging",]) > 0){
    save(sum_acc, freqs, file = paste0(path, "accelerateR/", bats[i], ".robj"))
  }
  if(nrow(ACC[ACC$behavior == "commuting"|ACC$behavior == "foraging",]) == 0){
    save(sum_acc, file = paste0(path, "accelerateR/", bats[i], "_no_freq.robj"))
  }
}
