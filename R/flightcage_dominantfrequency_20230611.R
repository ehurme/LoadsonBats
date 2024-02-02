library(pacman)

p_load(tidyverse, data.table, seewave, tuneR, dplyr)

## read data sheet
datasheet <- read.csv("./../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/WeightExperiments_Az23.csv")

datasheet |> filter(species == "E. fus") -> efus
efus$bat |> table()

base_path <- "../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/"
# bat 1
files_bat1 <- list.files(path = base_path, pattern = "Efus_bat1", recursive = TRUE, full.names = TRUE)
bat1 <- data.frame(bat = 1, files = files_bat1,
                   trial = 2:3,
                    bat_weight = c(12.34, 12.12),
                    tag_weight = c(0.85, 1.46))
# bat 48
files_bat48 <- list.files(path = base_path, pattern = "bat48", recursive = TRUE, full.names = TRUE)
bat48 <- data.frame(bat = 48, files = files_bat48,
                    trial = 2:4,
                    bat_weight = c(14.7, 15, 13.49),
                    tag_weight = c(1.49,0.8,2.14))
# bat 51
files_bat51 <- list.files(path = base_path, pattern = "bat51", recursive = TRUE, full.names = TRUE)
bat51 <- data.frame(bat = 51, files = files_bat51,
                    trial = 2:4,
                    bat_weight = c(16.72, 17.01, 16.96),
                    tag_weight = c(0.84, 2.13, 3.2))
# bat 53
files_bat53 <- list.files(path = base_path, pattern = "bat53", recursive = TRUE, full.names = TRUE)
bat53 <- data.frame(bat = 53, files = files_bat53,
                    trial = 2,
                    bat_weight = c(12.54),
                    tag_weight = c(1.03))
# bat 57
files_bat57 <- list.files(path = base_path, pattern = "bat57", recursive = TRUE, full.names = TRUE)
bat57 <- data.frame(bat = 57, files = files_bat57,
                    trial = 2,
                    bat_weight = c(15),
                    tag_weight = c(1.22))
# bat 58
files_bat58 <- list.files(path = base_path, pattern = "Epfu_bat58", recursive = TRUE, full.names = TRUE)
bat58 <- data.frame(bat = 58, files = files_bat58,
                    trial = 2:4,
                    bat_weight = c(12.45, 12.45, 12.04),
                    tag_weight = c(1.25,  2.16,  0.83))
# bat 62
files_bat62 <- list.files(path = base_path, pattern = "Epfu_bat62", recursive = TRUE, full.names = TRUE)
bat62 <- data.frame(bat = 62, files = files_bat62,
                    trial = 2:4,
                    bat_weight = c(12.75, 11.41, 11.87),
                    tag_weight = c(2.11,  0.83,  1.27))
# bat 63
files_bat63 <- list.files(path = base_path, pattern = "Epfu_bat63", recursive = TRUE, full.names = TRUE)
bat63 <- data.frame(bat = 63, files = files_bat63,
                    trial = 2:3,
                    bat_weight = c(13.18, 12.88),
                    tag_weight = c(1.71, 0.83))
# bat 64
files_bat64 <- list.files(path = base_path, pattern = "Epfu_bat64", recursive = TRUE, full.names = TRUE)
bat64 <- data.frame(bat = 64, files = files_bat64,
                    trial = 2:3,
                    bat_weight = c(14.23,13.76),
                    tag_weight = c(1.27, 1.72))

sampling_rate <- 54
window_length = sampling_rate*8

fc_freq <- rbind(bat1, bat48, bat51, bat53, bat57, bat58, bat62, bat63, bat64)

fc_freq <- fc_freq |> mutate(total_weight = bat_weight+tag_weight,
                          frequency = NA, frequency2 = NA,
                          peak_amp_diff = NA)

i = 1
for(i in 1:nrow(fc_freq)){
  df <- fread(fc_freq$files[i])

  # remove blank rows
  df <- df[!is.na(df$burstCount),]

  # add tag time
  df$time <- (1:nrow(df))/sampling_rate
  acc_duration <- max(df$time)

  # plot data
  layout(1)
  # with(df,#[1:1000,],#[9500:9700,],
  #      plot(time/60, accX_mg/1000, type = "l",
  #           xlab = "time (min)", ylab = "Acceleration (Gs)"))
  # lines(df$time/60, df$accY_mg/1000, col = rgb(1,0,0,.3))
  # lines(df$time/60, df$accZ_mg/1000, col = rgb(0,1,0,.3))
  # legend("bottomright", col = 1:3, lty = 1, legend = c("x", "y", "z"))


  # get dominant frequency from spectrogram
  df$z0 <- mean(df$accZ_mg, na.rm = TRUE) - df$accZ_mg
  w <- tuneR::Wave(left = na.omit(df$accZ_mg), samp.rate = sampling_rate, bit = 16)
  # plot(w)
  # spec <- meanspec(w, wl = window_length)
  # fpeaks(spec, nmax = 10)
  # spectro(w, wl = window_length)
  #dev.off()
  flapping <- ffilter(w, f= sampling_rate, from = 5, to = 20, bandpass = TRUE, wl = window_length)
  # spectro(flapping, f = sampling_rate)
  layout(1)
  spec <- meanspec(flapping, f=sampling_rate, wl = window_length*2, plot = FALSE)
  peak <- fpeaks(spec, nmax = 20, plot = FALSE)
  peak_order <- peak[order(peak[,2], decreasing = TRUE),]

  fc_freq$frequency[i] <- as.numeric(peak_order[1,1]) * 1000
  fc_freq$frequency2[i] <- as.numeric(peak_order[2,1]) * 1000
  fc_freq$peak_amp_diff[i] <- as.numeric(peak_order[1,2] - peak_order[2,2])
}
  # look at how wing beat frequency varies
  # step <- sampling_rate*15
  # intervals <- round(nrow(df)/(step))
  # j = 1
  #
  # for(j in 1:intervals){
  #   freq_int <- data.frame(bat_weight = bat_weight[i], tag_weight = tag_weight[i],
  #                          total_weight = bat_weight[i]+tag_weight[i],
  #                          interval = NA, time = NA, frequency = NA, frequency2 = NA,
  #                          peak_amp_diff = NA)
  #
  #   freq_int$interval <- j
  #   # freq_int$time
  #
  #   idx <-(((j-1)*step)+1):(j*step)
  #   df_int <- df[idx,]
  #   # get dominant frequency from spectrogram
  #   df_int$z0 <- mean(df_int$accZ_mg) - df_int$accZ_mg
  #   w <- tuneR::Wave(left = df_int$accZ_mg, samp.rate = sampling_rate, bit = 16)
  #   # plot(w)
  #   # spec <- meanspec(w, wl = window_length)
  #   # fpeaks(spec, nmax = 10)
  #   #spectro(w, wl = window_length)
  #   #dev.off()
  #   flapping <- ffilter(w, f= sampling_rate, from = 10, to = 20, bandpass = TRUE, wl = window_length)
  #   # spectro(flapping, f = sampling_rate, wl = window_length)
  #   layout(1)
  #   spec <- meanspec(flapping, f=sampling_rate, wl = step, plot = FALSE)
  #   peak <- {}
  #   try({peak <- fpeaks(spec, nmax = 20, plot = TRUE)
  #     if(length(peak) > 0){
  #       peak_order <- peak[order(peak[,2], decreasing = TRUE),]
  #
  #
  #       freq_int$frequency <- as.numeric(peak_order[1,1]) * 1000
  #       freq_int$frequency2 <- as.numeric(peak_order[2,1]) * 1000
  #       freq_int$peak_amp_diff <- as.numeric(peak_order[1,2] - peak_order[2,2])
  #
  #       fc_freq_int <- rbind(fc_freq_int, freq_int)
  #     }
  #   })
  # }

  # # manueverability
  # try({
  #   wf <- ffilter(w, f= sampling_rate, from = 0, to = 1, bandpass = TRUE, wl = length(w)/5)
  #   fc_freq$rms[i] <-  rms(df$z0)
  #   fc_freq$rms_filter[i] <- rms(wf*100)
  # })
# }
layout(1)
fc_freq <- fc_freq[order(fc_freq$total_weight),]
# plot(fc_freq$total_weight, fc_freq$frequency, col = bat,
#      cex = 2, pch = 16, type = "o", lwd = 1.5,
#      ylab = "wingbeat frequency (Hz)", xlab = "total bat weight (g)")
library(ggplot2)
ggplot(fc_freq, aes(x = total_weight, y = frequency,
                    group = factor(bat), col = factor(bat)))+
  geom_point(aes(size = peak_amp_diff))+
  geom_path(size = 1)+ theme_minimal()+
  xlab("Total mass (g)")+
  ylab("Wingbeat frequency (Hz)")+
  scale_color_viridis_d(name = "bat")

fc_freq$mean_freq <- (fc_freq$frequency + fc_freq$frequency2)/2
ggplot(fc_freq, aes(x = total_weight, y = mean_freq,
                    group = factor(bat), col = factor(bat)))+
  geom_point(aes(size = peak_amp_diff))+
  geom_path(size = 1)+ theme_minimal()+
  xlab("Total mass (g)")+
  ylab("Wingbeat frequency (Hz)")+
  scale_color_viridis_d(name = "bat")


ggplot(fc_freq, aes(x = total_weight, y = frequency,
                    group = factor(bat), col = factor(bat)))+
  geom_point(aes(size = peak_amp_diff))+
  geom_path(size = 1)+ theme_minimal()+
  xlab("Total mass (g)")+
  ylab("Wingbeat frequency (Hz)")+
  scale_color_viridis_d(name = "bat")

fc_freq$percent <- fc_freq$tag_weight/fc_freq$bat_weight
fc_freq <- fc_freq[order(fc_freq$percent),]
ggplot(fc_freq, aes(x = percent, y = mean_freq,
                    group = factor(bat), col = factor(bat)))+
  geom_point(aes(size = peak_amp_diff))+
  geom_path(size = 1)+ theme_minimal()+
  xlab("% body mass added")+
  ylab("Wingbeat frequency (Hz)")+
  scale_color_viridis_d(name = "bat")


fc_freq[which.min(fc_freq$frequency),]

fc_freq_int_clean <- fc_freq_int[fc_freq_int$frequency > 12 & fc_freq_int$frequency < 15,]

plot(fc_freq_int_clean$interval, fc_freq_int_clean$frequency, type = "o", ylim = c(8,15))


