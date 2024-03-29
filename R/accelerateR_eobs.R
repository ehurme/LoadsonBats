# AccelerateR loop to extract fft and params from acceleration data
# Edward Hurme
# 7/4/23


# devtools::install_github("wanjarast/accelerateR")

library(pacman)
p_load(data.table, janitor, accelerateR, move, signal, tuneR, seewave)
path <- "./../../../ownCloud/Firetail/Acerodonjubatus/tag_1521/"
path <- "./../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/"
path <- "../../../ownCloud/Firetail/Eidolonhelvum/Model_tag_2396/"

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
  plot(df$timestamp, df$speed)

  ACC <- {}
  ACC <- move_trans(data =  df[df$type == "acc",], timestamp = "timestamp", acc = "eobs_accelerations_raw",
                    id = "tag_local_identifier",
                    sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
                    naxes = 3, no_cores = 20)
  # check for any missing values
  burstcount = df$eobs_accelerations_raw[1] %>% strsplit(split = " ") %>% unlist %>% length/3
  count_test(ACC, burstcount)

  df$behavior <- "resting"
  df$behavior[which(df$`annotation-layer-commuting` != "")] <- "commuting"
  df$behavior[which(df$`annotation-layer-foraging` != "")] <- "foraging"

  ACC$behavior <- rep(df$behavior[df$type == "acc"],
                      each = burstcount)
  table(ACC$behavior)
  ACC$burst <- rep(1:(nrow(ACC)/burstcount),
                      each = burstcount)
  table(ACC$burst) %>% unique

  ## fast fourier transform
  if(nrow(ACC[ACC$behavior == "commuting" |
              ACC$behavior == "foraging",]) > 0){
    fft_acc <- {}
    ACC$time <- ACC$timestamp
    fft_acc <- sum_data(ACC, time = "time",
                        burstcount = burstcount,
                        #x="x" , y="y" ,
                        z="z" ,
                        stats = "FFT")

    png(file = paste0(path, "/accelerateR/fft_", bats[i], ".png"),
        width = 800, height = 600)
      image(fft_acc[,3:((burstcount/2)+1)] %>% as.matrix)
    dev.off()

    freqs <- data.frame(time = ACC$timestamp %>% unique,
                        freq = NA, amp = NA,
                        wfreq = NA, wamp = NA,
                        rms = NA, rms_filter = NA,
                        behavior = NA)
    j <- 5
    for(j in 1:nrow(freqs)){
      # get frequency
      idx <- which.max(fft_acc[j,3:((burstcount/2) + 1)])
      # abline(v = idx)
      freqs$freq[j] <- names(idx) %>% substr(3,nchar(names(idx)[1])) %>% as.numeric
      freqs$amp[j] <- max(fft_acc[j, 3:((burstcount/2)+1)])

      # measure peak frequency and maneuverability
      b1 <- ACC[ACC$burst == j,]
      b1$z0 <- b1$z - mean(b1$z)
      w <- tuneR::Wave(left = b1$z0, samp.rate = ACC$sample_frequency[1], bit = 16)
      # plot(w)
      try({
        spec <- meanspec(w, f=ACC$sample_frequency[1], wl = nrow(b1), plot = FALSE)
        peak <- fpeaks(spec, nmax = 4, plot = FALSE)
        pidx <- which(peak[,1] > 0.0005)
        midx <- which.max(peak[,2])
        if(length(midx) > 0){
          freqs$wfreq[j] <- peak[pidx[midx],1] * 1000
          freqs$wamp[j] <- peak[pidx[midx],2]
        }
      })

      wf <- ffilter(w, f= ACC$sample_frequency[1], from = 0, to = 1, bandpass = TRUE, wl = length(w)/5)
      freqs$rms[j] <-  rms(b1$z0)
      freqs$rms_filter[j] <- rms(wf*100)

      freqs$behavior[j] <- df$behavior[which(freqs$time[j] == df$timestamp)]
    }

    freqs$duration_of_burst <- ACC$burst_size[1]/ACC$sample_frequency[1]
    freqs$frequency <- freqs$freq/freqs$duration_of_burst

    plot(freqs$frequency, freqs$wfreq)
    plot(freqs$amp, freqs$wamp)

    png(file = paste0(path, "/accelerateR/freq_rms_", bats[i], ".png"),
        width = 800, height = 600)
    layout(rbind(1:2))
      with(freqs, plot(frequency, rms, cex = amp/max(amp, na.rm = TRUE), col = behavior %>% factor))
      with(freqs, plot(wfreq, rms, cex = wamp/max(wamp, na.rm = TRUE), col = behavior %>% factor))
    dev.off()

    freqs$duration_of_burst <- ACC$burst_size[1]/ACC$sample_frequency[1]

    png(file = paste0(path, "/accelerateR/wingbeat_", bats[i], ".png"),
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
