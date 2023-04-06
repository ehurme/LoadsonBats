## function to get the dominant frequency from ACC data during commutes

# df <- fread("./../../../ownCloud/Firetail/Phyllostomushastatus/Model_tag_7CE02AF_main/individual_7CDA644-annotated-bursts-gps.csv")
# df <- fread("./../../../ownCloud/Firetail/Phyllostomushastatus/Model_tag_7CE02AF_main/individual_7CDA644-annotated-bursts-gps.csv")
df <- fread("./../../../ownCloud/Firetail/Myotisvivesi/Mviv19_18_model/individual_Mviv19_18-annotated-bursts-gps.csv")
# calculate sunset and sunrise
# add wavelet
<<<<<<< HEAD
save_path = "./../../../ownCloud/Firetail/Myotisvivesi/Mviv19_18_model/Wingbeats/PCA/"
tag_id <- "Mviv19_18"
=======
save_path = "./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/Wingbeats/PCA/"
df <- fread("./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/tag_GPA-10_8147_S1-annotated-bursts-gps.csv") |> clean_names()
tag_id <- "8147"
>>>>>>> faae0a6d6a2e4195cbbac44fc86a21a5809e3aaa
Burst = TRUE
PCA = TRUE
sampling_rate = 100
  # M. lasiopterus 100
  # M. vivesi 50 - 2018,2019
  # M. vivesi 40 - 2017
  # P. hastatus 25
  # P. lylei 18.74
  # R. aegyptiacus 50
min_freq = 2
max_freq = 16 #8  # 4
# wavelet = FALSE
saved_cores = 20
gps = TRUE
min_seg_duration = 10
dfreq_threshold = 20
use_FFT = TRUE
# use_wavelet = FALSE
tag_type = "eObs"
  # "eObs" # "wildfi"
firetail = TRUE
firetail_filter = TRUE
solar_time = TRUE
save_files = TRUE
roll_var = TRUE
location = "Mexico"
var_seg = TRUE # segment the acc by variation
sd_adjust = 1 # how many standard deviations above the mean should the variation threshold be?

# Burst
# PCA
# sampling_rate
# min_freq
# max_freq
# min_seg_duration - minimum segment duration for determining the frequency
# wavelet
# gps - is there gps data and should we plot it?
# saved_cores
time <- {}
dominant_freq <- function(df,
                          GPS = NULL,
                          save_path = NULL,
                          tag_id = NULL,
                          Burst = TRUE,
                          PCA = TRUE,
                          sampling_rate = 25, # P. hastatus
                          min_freq = 2,
                          max_freq = 8,
                          # wavelet = FALSE,
                          saved_cores = 4,
                          gps = TRUE,
                          min_seg_duration = 10,
                          dfreq_threshold = 20,
                          use_FFT = TRUE,
                          # use_wavelet = FALSE,
                          tag_type = "eObs", # "wildfi"
                          firetail = TRUE,
                          firetail_filter = FALSE,
                          solar_time = TRUE,
                          save_files = TRUE,
                          location = "Mexico",
                          var_seg = TRUE,
                          sd_adjust = 1){
  # load libraries
  ## utilities
  require(pacman)
  p_load(data.table, magrittr, dplyr,
         janitor,
         foreach, doParallel, R.utils, roll,
         lubridate, suncalc,
         seewave, Rwave, tuneR,
         mclust)

  # source("../NoctuleMigration/scr/scan.track.map.r")

  # get timestamp with milliseconds
  op <- options(digits.secs=3)
  # set up parallelization
  cores=detectCores()
  cl <- makeCluster(cores[1]-saved_cores) #not to overload your computer
  registerDoParallel(cl)

  if(nrow(df) > 0){
    bats <- {}
    try(bats <- unique(tag_id))
    if(length(bats) == 0){
      bats <- "bat"
      df$id <- "bat"
    }
    try({
      if(length(df$id) == 0){
        df$id <- tag_id
      }})

    if(length(df$timestamp) == 0){
      # if(length(time) > 0){
      #   df$timestamp <- time
      # }
      if(length(df$`sample-timestamp`) > 0){
        df$timestamp <- df$`sample-timestamp`
      }
      # else{
      #   df$timestamp <- NULL
      # }
    }
    dates <- unique(date(df$timestamp))

    if(solar_time == TRUE){
      if(location == "Mexico"){
        data <- data.frame(date = dates, lat = 29,
                           lon = -113)
      }
      if(location == "Portugal"){
        data <- data.frame(date = dates, lat = 36.990,
                           lon = 6.440)
      }
      try(data <- data.frame(date = dates, lat = median(df$location_lat, na.rm = TRUE),
                             lon = median(df$location_long, na.rm = TRUE)))

      try({
        if(nrow(GPS) > 0){
          try(data <- data.frame(date = dates, lat = median(GPS$location_lat, na.rm = TRUE),
                                 lon = median(GPS$location_long, na.rm = TRUE)))
          }
        })

      sun <- suncalc::getSunlightTimes(data = data)
      moon <- suncalc::getMoonTimes(data = data)
      df$solar_time <- df$timestamp - hour(sun$solarNoon[1])*60*60
    }

    i = 1
    for(i in 1:length(bats)){
      full_path <- paste0(save_path, tag_id, "/")
      if(!dir.exists(full_path)){
        dir.create(full_path)
      }

      bdf <- df[df$id == bats[i],]
      dates <- unique(date(bdf$solar_time))
      # date(ymd_hms(bdf$timestamp)) |> table
      j = 2
      for(j in 1:length(dates)){
        d <- bdf[which(date(ymd_hms(bdf$solar_time)) == dates[j]),]

        if(gps == TRUE){
          # plot tracks
          try({
            lon <- d$location_long
            lat <- d$location_lat
            gps_time <- d$timestamp
            })
          try({
            if(nrow(GPS) > 0){
            try({
              lon <- GPS$location_long
              lat <- GPS$location_lat
              gps_time <- GPS$timestamp
              })
            }
          })
          # get sunrise and sunset times
          if(any(!is.na(lon))){
            try({
              png(file = paste0(full_path, "GPS_", bats[i],"_", dates[j], ".png"))
              layout(rbind(c(1, 2), c(1, 3)))
              par(mar = c(0, 4, 0, 0), oma = c(4, 0, 4, 4), xpd = NA)
              plot(lon, lat, asp = 1, type = "o", pch = 19, cex = 0.5)
              if(firetail == TRUE){
                commuting_idx <- which(d$annotation_layer_commuting != "" & d$type == "gps")
                foraging_idx <- which(d$annotation_layer_foraging != "" & d$type == "gps")
                resting_idx <- which(d$annotation_layer_resting != "" & d$type == "gps")
                with(d[commuting_idx,], points(location_long, location_lat, col = "green"))
                with(d[foraging_idx,], points(location_long, location_lat, col = "red"))
                with(d[resting_idx,], points(location_long, location_lat, col = "blue"))
                legend("topleft", legend = c("commuting", "foraging", "resting"),
                       col = c("green", "red", "blue"), pch = 16)
              }
              plot(gps_time, lon, type = "o", pch = 19, cex = 0.5,
                   xlab = "", xaxt = "none",
                   xlim = c(gps_time[1]-3*3600,
                            gps_time[length(gps_time)]+3*3600))

              lines(gps_time, lon)
              if(firetail == TRUE){
                with(d[commuting_idx,], points(timestamp, location_long, col = "green"))
                with(d[foraging_idx,], points(timestamp, location_long, col = "red"))
                with(d[resting_idx,], points(timestamp, location_long, col = "blue"))
              }
              abline(v = sun$sunset, col = "blue", xpd = FALSE)
              abline(v = sun$sunrise, col = "orange", xpd = FALSE)
              plot(gps_time, lat, type = "o", pch = 19, cex = 0.5,
                   xlim = c(gps_time[1]-3*3600,
                            gps_time[length(gps_time)]+3*3600))
              lines(gps_time, lat)
              if(firetail == TRUE){
                with(d[commuting_idx,], points(timestamp, location_lat, col = "green"))
                with(d[foraging_idx,], points(timestamp, location_lat, col = "red"))
                with(d[resting_idx,], points(timestamp, location_lat, col = "blue"))
              }
              abline(v = sun$sunset, col = "blue", xpd = FALSE)
              abline(v = sun$sunrise, col = "orange", xpd = FALSE)
              legend("bottomright", legend = c("sunset", "sunrise"), col = c("blue", "orange"), lty = 1)
              dev.off()
              # get speed

              # label commutes
              ## first passage time?
            })
          }
        }
        x <- {}
        y <- {}
        z <- {}
        time <- {}

        if(Burst == TRUE){
          burst <- {}

          if(nrow(d) > 100){
            ii = 1
            temp_df <- data.frame()
            db <- foreach(ii = 1:nrow(d), .combine = rbind) %dopar% {
              temp <- {}
              x <- NA
              y <- NA
              z <- NA
              time <- NA
              burst <- NA
              if(tag_type == "Wildfi"){
                temp <- d$accInGBurst[ii] |> strsplit(" ") |>
                  unlist() |>
                  as.numeric()
              }
              if(tag_type == "eObs"){
                temp <- d$eobs_accelerations_raw[ii] |> strsplit(" ") |>
                  unlist() |>
                  as.numeric()
              }

              if(length(temp) > 0){
                x <- temp[seq(1, length(temp), 3)]
                y <- temp[seq(2, length(temp), 3)]
                z <- temp[seq(3, length(temp), 3)]
                temp_time <- seq.POSIXt(from = d$timestamp[ii],
                                        to = d$timestamp[ii]+length(temp)/sampling_rate/3,
                                        length.out = length(temp)/3)
                time <- format(temp_time, "%Y-%m-%d %H:%M:%OS3") |> as.character()
                burst <- rep(ii, length(time))
              }
              temp_df <- data.frame(x,y,z,time,burst)
            }
            db <- na.omit(db)
            db$time <- ymd_hms(db$time)
            if(solar_time == TRUE){
              db$time <- db$time - hour(sun$solarNoon[1])*3600
            }
            x <- db$x
            y <- db$y
            z <- db$z
            time <- db$time
            burst <- db$burst
            if(firetail == TRUE){
              commuting_idx <- which(d$annotation_layer_commuting != "") # & d$type == "acc")
              foraging_idx <- which(d$annotation_layer_foraging != "") # & d$type == "acc")
              resting_idx <- which(d$annotation_layer_resting != "") # & d$type == "acc")

              db$firetail <- NA
              db$firetail[db$burst %in% commuting_idx] <- "c"
              db$firetail[db$burst %in% foraging_idx] <- "f"
              db$firetail[db$burst %in% resting_idx] <- "r"
              db$firetail |> table()
              behav <- db$firetail
              if(firetail_filter == TRUE){
                fire_idx <- which(behav != "r")
                x <- x[fire_idx]
                y <- y[fire_idx]
                z <- z[fire_idx]
                time <- time[fire_idx]
                db <- db[fire_idx,]
                table(db$firetail)
              }
            }
          }
        }
        if(Burst == FALSE){
          x <- d$ACCX
          y <- d$ACCY
          z <- d$ACCZ
          time <- df$time
          if(length(x) == 0){
            x <- d$`acceleration-x`
            y <- d$`acceleration-y`
            z <- d$`acceleration-z`
            time <- d$timestamp
          }
          if(length(x) == 0){
            x <- d$acceleration_x
            y <- d$acceleration_y
            z <- d$acceleration_z
            time <- d$timestamp
          }

          db <- data.frame(x,y,z,time)
          db$time <- ymd_hms(db$time)
          if(solar_time == TRUE){
            db$time <- db$time - hour(sun$solarNoon[1])*3600
          }
        }

        # check for breaks in recording
        check_break <- TRUE
        if(check_break == TRUE){
          acc_break_threshold <- 1 # in seconds
          acc_breaks <- table(cut(1:length(time),labels = FALSE,
                                  breaks = c(1, which(diff(time) > acc_break_threshold),
                                             length(time))))
          acc_break_max_idx <- which.max(acc_breaks)

          acc_start <- sum(acc_breaks[1:(acc_break_max_idx)-1])
          acc_end <- sum(acc_breaks[1:(acc_break_max_idx)])
          db <- db[acc_start:acc_end,]
          x <- x[acc_start:acc_end]
          y <- y[acc_start:acc_end]
          z <- z[acc_start:acc_end]
        }
        # !! add a feature to cycle through breaks longer than a certain amount
        # probably useful for Phast


        if(length(x) > 1){
          # pca
          png(file = paste0(full_path, "FullSpectro_", bats[i],"_", dates[j], ".png"))
          if(PCA == TRUE){
            try({
              pc <- prcomp(cbind(x, y, z), scale. = FALSE)
              w <- Wave(left = pc$x[,1], samp.rate = sampling_rate, bit = 16)
              if(roll_var == TRUE){
                #plot(time, pc$x[,1], type = "l")
                rz <- roll_var(pc$x[,1], width = 100000, complete_obs = TRUE, na_restore = TRUE)
                # lines(rz, col = 2)
                c <- Mclust((rz |> na.omit()), G = 2)
                # plot(c, what = "density")
                thresh <- c$parameters$mean[2] - sqrt(c$parameters$variance$sigmasq[2])
                roll_idx <- which(rz > thresh)
                roll_time <- time[roll_idx]
                #points(roll_time, rz[roll_idx], col = 3)
                w <- Wave(left = pc$x[roll_idx,1], samp.rate = sampling_rate, bit = 16)
              }
              wf <- ffilter(w, f= sampling_rate, from = min_freq, to = max_freq, bandpass = TRUE)
              spectro(wf, f = sampling_rate, wl = 1024*4, ovlp = 50, fastdisp = TRUE)
              # add sunrise and sunset?
            })
          }
          if(PCA == FALSE){
            # create wave form
            w <- Wave(left = z, samp.rate = sampling_rate, bit = 16)
            wf <- ffilter(w, f= sampling_rate, from = min_freq, to = max_freq, bandpass = TRUE)
            spectro(wf, f = sampling_rate, wl = 1024*4, ovlp = 50, fastdisp = TRUE)
          }
          par(new=TRUE)
          pf <- dfreq(wf, f=sampling_rate, wl = 1024*4, ovlp=75, threshold=dfreq_threshold,
                      type="l", col="red", lwd=0.5, xlab = "", ylab = "")
          dev.off()

          # get segments where wing beat can be estimated
          idx <- which(!is.na(pf[,2])) -1
          stable_var <- seqToIntervals(idx) |> as.data.frame()
          stable_var$diff <- NA

          # if track is all one segment, divide track by variation in wingbeat
          ## or if you manually want to divide by track by variation
          if(nrow(stable_var) < 2 | var_seg == TRUE){
            # divide track by variation in freq
            pf_var <- roll_var(pf[,2], width = 10)
            if(any(!is.na(pf_var))){
              png(file = paste0(full_path, "VarSeg_", bats[i],"_", dates[j], ".png"),
                  height = 600, width = 1000)
              layout(cbind(1,2))
              hist(pf_var, breaks = length(pf_var)/10, main = "")

              big <- 1000000
              c <- Mclust((pf_var |> na.omit())*big, G = 2)
              # plot(c, what = "density")
              thresh <- c$parameters$mean[1]/big + sd_adjust*sqrt(c$parameters$variance$sigmasq[1])/big
              abline(v = thresh, col = 2)
              #layout(1)
              # plot(pf_var, col = (pf_var > thresh)+1)
              plot(time, pf[,2], xlab = "Time (sec)",
                   col = (pf_var < thresh & pf[,2] > min_freq/1000 & pf[,2] < max_freq/1000)+1,
                   pch = 16)
              idx <- which(pf_var < thresh & pf[,2] > min_freq/1000 & pf[,2] < max_freq/1000)
              stable_var <- seqToIntervals(idx) |> as.data.frame()
              abline(v = pf[stable_var$from])
              abline(v = pf[stable_var$to], col = "orange")
              dev.off()
            }
          }

          reg_segs <- 10 # regular segmentation (min)
          reg_seg <- FALSE
          if(nrow(stable_var) < 2 | reg_seg == TRUE){
            seqs <- c(seq(1, length(time), by = 10*sampling_rate), )

          }

          ## join periods that are less than # samples apart
          min_diff = 5 # what are the time units on this?
          kk = 1
          for(kk in 1:nrow(stable_var)){
            stable_var$diff[kk] <- stable_var$from[kk+1] - stable_var$to[kk]
          }
          while(any(stable_var$diff < min_diff, na.rm = TRUE)){
            rm_idx <- {}
            for(kk in 1:(nrow(stable_var)-1)){
              if(stable_var$diff[kk] < min_diff){
                stable_var$from[kk+1] <- stable_var$from[kk]
                rm_idx <- c(rm_idx, kk)
              }
            }
            stable_var <- stable_var[-rm_idx,]
            for(kk in 1:nrow(stable_var)){
              stable_var$diff[kk] <- stable_var$from[kk+1] - stable_var$to[kk]
            }
          }

          stable_var$length <- stable_var$to - stable_var$from + 1 # adding one includes all parts of the segment

          ## create
          # multiply by sampling interval of pf
          w_stable <- stable_var * (pf[2,1])
          df_stable <- stable_var * (pf[2,1]) * sampling_rate

          if(roll_var){
            # fix to work for multiple intervals
            df_stable$from <- df_stable$from + roll_idx[1]
            df_stable$to <- df_stable$to + roll_idx[1]
          }
          ## fix last value??
          # df_stable$to[nrow(df_stable)] <- nrow(db)
          # nrow(db)-6233 * (pf[2,1]) * sampling_rate

          # fix first value
          if(w_stable[1,1] == 0) w_stable[1,1] <- 1
          if(df_stable[1,1] == 0) df_stable[1,1] <- 1

          # get FFT peak freq
          if(use_FFT == TRUE){
            try({
              w_stable$peak_freq <- NA
              w_stable$peak_amp_diff <- NA # difference from 1st peak amplitude to 8th peak
              w_stable$amplitude <- NA
              w_stable$firetail <- NA
              # jj = which.max(w_stable$length)
              jj = 1
              for(jj in 1:nrow(w_stable)){
                if(firetail == TRUE){
                  tf <- table(db$firetail[df_stable$from[jj]:df_stable$to[jj]])
                  w_stable$firetail[jj] <- names(tf[which.max(tf)])
                }
                if(stable_var$length[jj] > 1){
                  # filter short periods
                  if(difftime(db$time[df_stable$to[jj]], db$time[df_stable$from[jj]], units = "mins") > min_seg_duration){
                    if(PCA == TRUE){
                      w_cut <- Wave(left = pc$x[,1][df_stable$from[jj]:df_stable$to[jj]] |> na.omit(),
                                    samp.rate = sampling_rate, bit = 16)
                    }
                    if(PCA == FALSE){
                      w_cut <- Wave(left = db$z[df_stable$from[jj]:df_stable$to[jj]] |> na.omit(),
                                    samp.rate = sampling_rate, bit = 16)
                    }

                    png(file = paste0(full_path, "Wav_", bats[i],"_", dates[j],"_segment_", jj, ".png"))
                      layout(1)
                      plot(w, main = paste0(bats[i]," ", dates[j]," segment ", jj))
                      abline(v = w_stable$to[jj], col = 2, lwd = 2)
                      abline(v = w_stable$from[jj], col = 3, lwd = 2)
                    dev.off()
                    wf_cut <- ffilter(w_cut, f= sampling_rate, from = min_freq, to = max_freq, bandpass = TRUE)

                    png(file = paste0(full_path, "Spectro_", bats[i],"_", dates[j],"_segment_", jj, ".png"))
                      spectro(wf_cut, f = sampling_rate, xpd = FALSE)
                    dev.off()
                    spec <- meanspec(wf_cut, f=sampling_rate, plot = FALSE)
                    layout(1)
                    png(file = paste0(full_path, "fpeaks_", bats[i],"_", dates[j],"_segment_", jj, ".png"))
                      peak <- fpeaks(spec, nmax = 8)
                    dev.off()
                    pidx <- which(peak[,1] > 0.0025)
                    midx <- which.max(peak[pidx,2])
                    w_stable$peak_freq[jj] <- peak[pidx[midx],1]
                    w_stable$peak_amp_diff[jj] <- peak[,2] |> range() |> diff() # doesn't seem like a useful value
                    # df$ACCZ[df_stable$from[i]:df_stable$to[i]] |> hist()
                    maxz <- roll_max(db$z[df_stable$from[jj]:df_stable$to[jj]], width = 6) |>
                      median(na.rm = TRUE)
                    minz <- roll_min(db$z[df_stable$from[jj]:df_stable$to[jj]], width = 6) |>
                      median(na.rm = TRUE)
                    # abline(v = c(minz, maxz), col = 2, lty = 2)
                    w_stable$amplitude[jj] <- maxz - minz
                  }
                }
              }
            })
          }

          if(any(!is.na(w_stable$peak_freq))){
            try({
              png(filename = paste0(full_path, "Domfreq_segments_",  bats[i], "_", dates[j],".png"))
              layout(1)
              xlab = "time (UTC)"
              if(solar_time == TRUE){
                xlab = "Time since solar noon"
              }
              if(firetail == TRUE){
                w_stable$firetail <- factor(w_stable$firetail,
                                            levels = c("c", "f", "r"))
                plot(x = ymd_hms(db$time[df_stable$from]), y = w_stable$peak_freq*1000,
                     cex = 2*scales::rescale(w_stable$amplitude)+1,
                     pch = 16,
                     col = w_stable$firetail,
                     xlab = xlab, ylab = "dominant wingbeat frequency (Hz)")
              }
              if(firetail == FALSE){
                plot(x = ymd_hms(db$time[df_stable$from]), y = w_stable$peak_freq*1000,
                     cex = 2*scales::rescale(w_stable$amplitude),
                     pch = 16,
                     col = 2,
                     xlab = xlab, ylab = "dominant wingbeat frequency (Hz)")
              }
              segments(x0 = ymd_hms(db$time[df_stable$from]),
                       y0 = w_stable$peak_freq*1000,
                       x1 = ymd_hms(db$time[df_stable$to]),
                       y1 = w_stable$peak_freq*1000)
              if(solar_time == TRUE){
                abline(v = sun$sunset - hour(sun$solarNoon[1])*3600, col = "blue")
                abline(v = sun$sunrise - hour(sun$solarNoon[1])*3600, col = "orange")
              }
              if(solar_time == FALSE){
                abline(v = sun$sunset, col = "blue")
                abline(v = sun$sunrise, col = "orange")
              }
              if(firetail == TRUE){
                legend("topright", legend = c("commuting", "foraging", "resting"),
                       col = c(1:3), pch = 16)
              }
              dev.off()
            })
          }

          if(save_files == TRUE){
            save(db, w, wf, pf, sampling_rate, stable_var, w_stable, df_stable, bats, dates,
                 save_path,
                 tag_id,
                 Burst,
                 PCA,
                 sampling_rate,
                 min_freq,
                 max_freq,
                 #wavelet,
                 saved_cores,
                 gps,
                 min_seg_duration,
                 dfreq_threshold,
                 use_FFT,
                 #use_wavelet,
                 tag_type,
                 firetail,
                 firetail_filter,
                 solar_time,
                 save_files,
                 location,
                 var_seg,
                 sd_adjust,
                 file = paste0(full_path, bats[i], "_", dates[j],"_wingbeatfreq.robj"))
          }
        }
      }
    }
  }
}
