library(pacman)

p_load(data.table, seewave, tuneR, stringr, lubridate, ggplot2)

path <- "../../../Dropbox/test weight/test weight/27072023/L1_m_trial1"
capture_sheet <- read.csv(paste0(path, "WeightExperiments_Az23.csv"))

files1 <- dir(path, pattern = ".csv", recursive = TRUE)
files2 <- dir(path, pattern = "Epfu", recursive = TRUE)
files <- c(files1, files2)

dates <- ymd(substr(files, 1, 8))
bats <- str_extract(pattern = "(?<=bat)[0-9]+", string = files)
trials <- str_extract(pattern = "(?<=trial)[0-9]+", string = files)
sampling_rate <- 105
window_length = sampling_rate*8
burst_length = 156

fc_freq <- data.frame(file = files, date = dates, bat = bats, trial = trials,
                      bat_weight = NA, tag_weight = NA, total_weight = NA,
                      sd = NA, median = NA,
                      frequency = NA, frequency2 = NA, peak_amp_diff = NA)

i = 1
for(i in 1:length(files)){
  try({
    fc_freq$bat_weight[i] <- capture_sheet$bat.weight[which(capture_sheet$bat == fc_freq$bat[i] &
                                                              capture_sheet$species == "E. fus" &
                                                              capture_sheet$trial == fc_freq$trial[i] &
                                                              capture_sheet$date == fc_freq$date[i])] %>% as.numeric()
    fc_freq$tag_weight[i] <- capture_sheet$tag.weight[which(capture_sheet$bat == fc_freq$bat[i] &
                                                              capture_sheet$species == "E. fus" &
                                                              capture_sheet$trial == fc_freq$trial[i] &
                                                              capture_sheet$date == fc_freq$date[i])] %>% as.numeric()
    fc_freq$total_weight[i] <- fc_freq$bat_weight[i] + fc_freq$tag_weight[i]
  })

  df <- fread(paste0(path, files[i]))

  df <- fread("C:/Users/edwar/Dropbox/MPI/Wingbeat/Belgium23/Nyctalus_leisleri_fleatag/27072023/L1_m_trial1/ACC_L1_m_trial1_tag1.25_bat11.6_8G_105Hz.csv")
  df <- fread("C:/Users/edwar/Dropbox/MPI/Wingbeat/Belgium23/Nyctalus_leisleri_fleatag/27072023/L1_m_trial2/ACC_L1_m_trial2_tag0.76_bat12.47_8G_105Hz.csv")
  df <- fread("C:/Users/edwar/Dropbox/MPI/Wingbeat/Belgium23/Nyctalus_leisleri_fleatag/27072023/L1_m_trial3/ACC_L1_m_trial3_tag1.02_bat12.5_8G_105Hz.csv")
  df <- fread("C:/Users/edwar/Dropbox/MPI/Wingbeat/Belgium23/Nyctalus_leisleri_fleatag/28072023/bat1_f_trial1/ACC_bat1_f_trial1_tag1.00_bat13.44_8g_105Hz.csv")
  df <- fread("C:/Users/edwar/Dropbox/MPI/Wingbeat/Belgium23/Nyctalus_leisleri_fleatag/28072023/bat1_f_trial2/ACC_bat1_f_trial2_tag1.25_bat13.74_8g_105Hz.csv")

  # add data by burst
  max_burst <- as.numeric(df$burstCount[nrow(df)])
  d <- data.frame(time = 1:((max_burst+1) * burst_length)/sampling_rate,
                  burst = rep(0:max_burst, each = burst_length),
                  x = NA, y = NA, z = NA)

  bursts <- unique(df$burstCount) %>% as.numeric() %>% na.omit()
  bursts <- bursts[bursts >= 0 & bursts =< max_burst]
  for(i in 1:length(bursts)){
    idx <- which(df$burstCount == bursts[i])
    didx <- which(d$burst == bursts[i])
    if(length(idx) != length(didx)){
      d$x[didx[1:length(idx)]] <- df$accX_mg[idx] %>% as.numeric
      d$y[didx[1:length(idx)]] <- df$accY_mg[idx] %>% as.numeric
      d$z[didx[1:length(idx)]] <- df$accZ_mg[idx] %>% as.numeric
    }
    if(length(idx) == length(didx)){
      d$x[didx] <- df$accX_mg[idx] %>% as.numeric
      d$y[didx] <- df$accY_mg[idx] %>% as.numeric
      d$z[didx] <- df$accZ_mg[idx] %>% as.numeric
    }
  }
  summary(d)

  d$x[which(abs(d$x) > 8000)] <- 8000
  d$y[which(abs(d$y) > 8000)] <- 8000
  d$z[which(abs(d$z) > 8000)] <- 8000

  # plot data
  layout(1)
  with(d, #[9500:9700,],
       plot(time/60, x/1000, type = "l",
            xlab = "time (min)", ylab = "Acceleration (Gs)",
            ylim = c(-8,8)))
  lines(d$time/60, d$y/1000, col = rgb(0,1,0,.4))
  lines(d$time/60, d$z/1000, col = rgb(1,0,0,.4))

  # get dominant frequency from spectrogram
  # pl.pc <- prcomp(with(df, cbind(accX_mg, accY_mg, accZ_mg)), scale. = FALSE)
  # df$pc1 <- pl.pc$x[,1]
  # plot(df$time, df$pc1, type = "l")

  d$z0 <- mean(d$z, na.rm = TRUE) - d$z
  d$z0[which(is.na(d$z0))] <- 0
  # summary(df$z0)
  # plot(df$z0)
  w <- tuneR::Wave(left = d$z0, samp.rate = sampling_rate, bit = 16)
  # w_pc <- tuneR::Wave(left = df$pc1, samp.rate = sampling_rate, bit = 16)
  # plot(w)
  # spec <- meanspec(w, wl = window_length)
  # fpeaks(spec, nmax = 10)
  # spectro(w, wl = window_length)
  #dev.off()
  flapping <- ffilter(w, f= sampling_rate, from = 10, to = 25, bandpass = TRUE,
                      wl = window_length)
  # flapping_pc <- ffilter(w_pc, f= sampling_rate, from = 8, to = 15, bandpass = TRUE, wl = window_length)
  spectro(flapping, f = sampling_rate)
  # spectro(flapping_pc, f = sampling_rate)
  layout(1)

  # csh(flapping, f = sampling_rate)
  try({
    dfrq <- dfreq(flapping, sampling_rate, clip = .5)
    fc_freq$sd[i] <- (dfrq[,2]*1000) %>% na.omit %>% sd
    fc_freq$median[i] <- (dfrq[,2]*1000) %>% na.omit %>% median
  })

  spec <- meanspec(flapping, f=sampling_rate, wl = window_length, plot = FALSE)
  peak <- fpeaks(spec, plot = TRUE, threshold = 0.5, freq = 1)
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
fc_freq <- fc_freq[order(fc_freq$bat, fc_freq$trial),]
plot(fc_freq$total_weight, fc_freq$frequency,
     cex = 2, pch = 16, col = fc_freq$bat %>% factor,
     ylab = "wingbeat frequency (Hz)", xlab = "total bat weight (g)")

ggplot(fc_freq, aes(total_weight, frequency, col = bat))+
  geom_point()+
  geom_path()

ggplot(fc_freq, aes(tag_weight, frequency, col = bat))+
  geom_point()+
  geom_path()

ggplot(fc_freq, aes(total_weight, median, col = bat))+
  geom_point()+
  geom_path()
# geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.2,
#               position=position_dodge(0.05))


plot(fc_freq$frequency, fc_freq$median)
plot(fc_freq$frequency, fc_freq$frequency2)


fit <- glm(frequency ~ total_weight + tag_weight + trial, data = fc_freq)
summary(fit)
