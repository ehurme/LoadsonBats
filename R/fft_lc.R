# look at stork wingbeat frequency

library(pacman)
p_load(data.table, tidyverse, magrittr, dplyr, accelerateR, stringr, seewave, tuneR)
op <- options(digits.secs=3)
source("src/fft_peak_freq.R")

path <- "../../../Dropbox/MPI/Wingbeat/data/Lasiurus_cinereus/DataDepository/Accelerometer Files/"
files <- list.files(path, full.names = TRUE, pattern = ".csv")
IDs <- extracted_substrings <- gsub(".*/Full_ACC_(\\w+)\\.csv", "\\1", files)

sampling_rate <- 40
lowpass <- 5
highpass <- 15

lcb <- data.frame()

j = 1
for(j in 1:length(files)){
  lc <- fread(files[j])
  lc <- lc[,1:4]
  names(lc) <- c("x","y","z","time")

  # Convert MATLAB time to POSIXct object in R
  origin <- as.POSIXct("0000-01-01", format = "%Y-%m-%d", tz = "UTC") # Define MATLAB origin
  lc$timestamp <- origin + as.difftime(lc$time, units = "days") # Convert MATLAB time to POSIXct
  lc$timestamp %>% diff -> lc_diff

  idx <- c(1,which(lc_diff > 10)+1)
  i = 10

  l <- data.frame()
  for(i in 1:(length(idx)-1)){
    b <- lc[idx[i]:(idx[i+1]-1),]
    pc <- prcomp(x = b[,1:3])
    b$pc <- pc$x[,1]

    plot(b$timestamp, b$z, type = "l", ylim = c(-40000,4e4))
    lines(b$timestamp, b$x, type = "l", col = 2)
    lines(b$timestamp, b$y, type = "l", col = 3)
    lines(b$timestamp, b$pc, type = "l", col = 4)

    duration <- difftime(max(b$timestamp), min(b$timestamp), units = "secs") %>% as.numeric()

    df <- fft_peak_freq(data = b$pc,
                        time = duration,
                        highpass = highpass,
                        lowpass = lowpass,
                        sampling_rate = 40)
    df$burst <- i
    df$burst_start <- b$timestamp[1]
    df$duration <- duration
    l <- rbind(l, df)
  }
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

