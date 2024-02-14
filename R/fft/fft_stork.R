# look at stork wingbeat frequency
# sampling 4s acceleration every minute

library(pacman)
p_load(data.table, tidyverse, magrittr, dplyr,
       accelerateR,move,
       stringr,
       seewave, tuneR)
op <- options(digits.secs=3)
source("src/fft_peak_freq.R")
source("src/get_day_or_night.R")

path <- "../../../ownCloud/Firetail/WhiteStor/Affenburg/"
files <- list.files(path, full.names = TRUE, pattern = ".csv")
IDs <- extracted_substrings <- gsub(".*/individual_(\\w+)\\-annotated-bursts-gps.csv", "\\1", files)


storks <- data.frame()

j = 1
for(j in 1:length(files)){
  s <- fread(files[j])
  summary(s)
  s$day_or_night <- NA
  days <- unique(date(s$timestamp))
  i = 1
  for(i in 1:length(days)){
    idx <- which(date(s$timestamp) == days[i])
    s$day_or_night[idx] <- get_day_or_night(latitude = mean(s$location_lat, na.rm = TRUE),
                                            longitude = mean(s$location_long, na.rm = TRUE),
                                            timestamp = s$timestamp[idx])
  }

  s$behavior <- NA
  s$behavior[s$`annotation-layer-walking` != ""] <- "walking"
  s$behavior[s$`annotation-layer-soaring` != ""] <- "soaring"
  s$behavior[s$`annotation-layer-flapping` != ""] <- "flapping"

  acc <- s[s$type == "acc" & s$`annotation-layer-flapping` != "",]
  ACC <- move_trans(data =  acc,
                    timestamp = "timestamp", acc = "eobs_accelerations_raw",
                    id = "tag-serial-number",
                    sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
                    naxes = 3, no_cores = 20)
  burst_size <- unique(ACC$burst_size)
  ACC$burst <- rep(which(s$type == "acc" & s$`annotation-layer-flapping` != ""),
                   each = burst_size)
  ACC$daynight <- rep(s$day_or_night[s$type == "acc" & s$`annotation-layer-flapping` != ""],
                      each = burst_size)
  ACC$behavior <- rep(s$behavior[s$type == "acc" & s$`annotation-layer-flapping` != ""],
                      each = burst_size)

  s_diff <- ACC$timestamp |> diff()
  idx <- c(1,which(s_diff > 10)+1)

  i = 2500
  stork <- data.frame()
  for(i in 1:(length(idx)-1)){
    b <- ACC[idx[i]:(idx[i+1]-1),]
    pc <- prcomp(x = b[,1:3])
    b$pc <- pc$x[,1]

    plot(1:nrow(b), b$z, type = "l", ylim = c(-4000,4e3))
    lines(1:nrow(b), b$x, type = "l", col = 2)
    lines(1:nrow(b), b$y, type = "l", col = 3)
    lines(1:nrow(b), b$pc, type = "l", col = 4)

    duration <- nrow(b)/b$sample_frequency[1]
    try({
      df <- fft_peak_freq(data = b$pc,
                          time = duration,
                          highpass = 3,
                          lowpass = 0.1,
                          sampling_rate = b$sample_frequency[1])
      df$burst <- b$burst[1]
      df$burst_start <- b$timestamp[1]
      df$duration <- duration
      df$behavior <- acc$`annotation-layer-flapping`[i]
      stork <- rbind(stork, df)
    })
  }
  ggplot(stork, aes(burst_start, freq))+
    geom_point(alpha = 0.1)

  l$bat <- bats[j]
  l$sampling_rate <- 40
  burst_diff <- l$burst_start %>% diff
  next_day <- {}
  next_day <- which(burst_diff > 100)
  l$day <- 1
  if(length(next_day) > 0){
    l$day[(next_day+1):nrow(l)] <- 2
  }
  lcb <- rbind(lcb, l)
}

# save data
save(lcb, file = "../../../Dropbox/MPI/Wingbeat/data/Frequency_Lasiurus_cinereus.robj")

summary(lcb)

ggplot(lcb[lcb$amp > 5e5,], # & lcb$peak_freq > 7,],
       aes(burst_start, freq, size = amp))+
  geom_point(alpha = 0.5)+
  scale_size(name = "Amplitude", range = c(.1, 5))+facet_wrap(~bat+day, scales = "free_x")+
  theme_bw()

ggplot(lcb, aes(freq, amp, col = bat))+geom_point(alpha = 0.2)+
  geom_hline(yintercept = 1e5)+
  geom_hline(yintercept = 3e5)+
  geom_vline(xintercept = 7)

ggplot(lcb[lcb$peak_amp > 5e5,], # & lcb$peak_freq > 7,],
       aes(hour(burst_start), peak_freq))+
  geom_point(alpha = 0.5)+
  geom_smooth()

