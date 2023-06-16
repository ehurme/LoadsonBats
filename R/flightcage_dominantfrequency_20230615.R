library(pacman)

p_load(data.table, seewave, tuneR)

path <- "../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230615/"
folders <- list.dirs(path, full.names = TRUE)
folders <- folders[2:4]

sampling_rate <- 54
window_length = sampling_rate*8

i = 3
tag_weight <- c(0.85, 1.46, 0.82)
bat_weight <- c(12.34, 12.12, 9.81)

fc_freq <- data.frame(folders, bat_weight, tag_weight, total_weight = bat_weight+tag_weight,
                      frequency = NA, frequency2 = NA, peak_amp_diff = NA)
                      # frequency_pc = NA, frequency2_pc = NA, peak_amp_diff_pc = NA)
# fc_freq_int <- data.frame()
#rms = NA, rms_filter = NA)
i = 1
for(i in 1:length(folders)){
  files <- list.files(folders[i], pattern = ".csv", full.names = TRUE)
  df <- fread(files[1])

  # remove blank rows
  df <- df[!is.na(df$burstCount),]

  # add tag time
  df$time <- (1:nrow(df))/sampling_rate
  acc_duration <- max(df$time)

  # plot data
  # layout(1)
  # with(df,#[1:1000,],#[9500:9700,],
  #      plot(time/60, accX_mg/1000, type = "l",
  #           xlab = "time (min)", ylab = "Acceleration (Gs)"))
  # lines(df$time/60, df$accY_mg/1000, col = 2)
  # lines(df$time/60, df$accZ_mg/1000, col = 3)


  # get dominant frequency from spectrogram
  # pl.pc <- prcomp(with(df, cbind(accX_mg, accY_mg, accZ_mg)), scale. = FALSE)
  # df$pc1 <- pl.pc$x[,1]
  # plot(df$time, df$pc1, type = "l")

  df$z0 <- mean(df$accZ_mg) - df$accZ_mg

  w <- tuneR::Wave(left = df$accZ_mg, samp.rate = sampling_rate, bit = 16)
  # w_pc <- tuneR::Wave(left = df$pc1, samp.rate = sampling_rate, bit = 16)
  # plot(w)
  # spec <- meanspec(w, wl = window_length)
  # fpeaks(spec, nmax = 10)
  # spectro(w, wl = window_length)
  #dev.off()
  flapping <- ffilter(w, f= sampling_rate, from = 8, to = 15, bandpass = TRUE, wl = window_length)
  # flapping_pc <- ffilter(w_pc, f= sampling_rate, from = 8, to = 15, bandpass = TRUE, wl = window_length)
  spectro(flapping, f = sampling_rate)
  # spectro(flapping_pc, f = sampling_rate)
  layout(1)
  spec <- meanspec(flapping, f=sampling_rate, wl = window_length, plot = FALSE)
  peak <- fpeaks(spec, nmax = 20, plot = TRUE)
  peak_order <- peak[order(peak[,2], decreasing = TRUE),]

  # spec_pc <- meanspec(flapping_pc, f=sampling_rate, wl = window_length, plot = FALSE)
  # peak_pc <- fpeaks(spec_pc, nmax = 20, plot = TRUE)
  # peak_pc_order <- peak_pc[order(peak_pc[,2], decreasing = TRUE),]

  fc_freq$frequency[i] <- as.numeric(peak_order[1,1]) * 1000
  fc_freq$frequency2[i] <- as.numeric(peak_order[2,1]) * 1000
  fc_freq$peak_amp_diff[i] <- as.numeric(peak_order[1,2] - peak_order[2,2])

  # fc_freq$frequency_pc[i] <- as.numeric(peak_order[1,1]) * 1000
  # fc_freq$frequency2_pc[i] <- as.numeric(peak_order[2,1]) * 1000
  # fc_freq$peak_amp_diff_pc[i] <- as.numeric(peak_order[1,2] - peak_order[2,2])
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
plot(fc_freq$total_weight, fc_freq$frequency,
     cex = 2, pch = 16, type = "o", lwd = 1.5,
     ylab = "wingbeat frequency (Hz)", xlab = "total bat weight (g)")

fc_freq_int_clean <- fc_freq_int[fc_freq_int$frequency > 12 & fc_freq_int$frequency < 15,]

plot(fc_freq_int_clean$interval, fc_freq_int_clean$frequency, type = "o", ylim = c(8,15))
ggplot(fc_freq_int, aes(x = interval, y = frequency, group = total_weight, col = total_weight))+
  geom_point()+
  geom_path()


