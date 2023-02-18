## function to get the dominant frequency from ACC data during commutes

df <- fread("./../../../ownCloud/Firetail/Phyllostomushastatus/Model_tag_7CE02AF_main/individual_7CDA644-annotated-bursts-gps.csv")
# calculate sunset and sunrise
# add wavelet
Burst = TRUE
PCA = TRUE
sampling_rate = 50 # P. hastatus
  # P. lylei 18.74
min_freq = 2
max_freq = 8  # 4
wavelet = TRUE
saved_cores = 20
gps = FALSE
min_seg_duration = 30
dfreq_threshold = 20
use_FFT = TRUE

# Burst
# PCA
# sampling_rate
# min_freq
# max_freq
# min_seg_duration - minimum segment duration for determining the frequency
# wavelet
# gps - is there gps data and should we plot it?
# saved_cores
time <- dmy_hms(df$utcDate)
dominant_freq <- function(df, tag_id = NULL, time = NULL, Burst = FALSE,
                          tag_type == "Wildfi",
                          PCA = FALSE, sampling_rate = 25,
                          min_freq = 2, max_freq = 8,
                          min_seg_duration = 10,
                          wavelet = FALSE, gps = FALSE, saved_cores = 4){
  # load libraries
  ## utilities
  require(data.table)
  require(magrittr)
  require(foreach)
  require(doParallel)
  require(R.utils)
  require(roll)
  # time
  require(lubridate)
  require(suncalc)
  # movement analysis
  require(marcher)
  # sound analysis
  require(seewave)
  require(Rwave)
  require(tuneR)
  require(dplR)

  source("../NoctuleMigration/scr/scan.track.map.r")

  # get timestamp with milliseconds
  op <- options(digits.secs=3)
  # set up parallelization
  cores=detectCores()
  cl <- makeCluster(cores[1]-saved_cores) #not to overload your computer
  registerDoParallel(cl)

  bats <- unique(tag_id)
  if(length(bats) == 0){
    bats <- "bat"
    df$id <- "bat"
  }

  if(length(df$timestamp) == 0){
    if(length(time) > 0){
      df$timestamp <- time
    }
    else{
      df$timestamp <- NULL
    }
  }
  d <- unique(date(df$timestamp))
  data <- data.frame(date = d, lat = median(df$location_lat, na.rm = TRUE),
                     lon = median(df$location_long, na.rm = TRUE))
  sun <- suncalc::getSunlightTimes(data = data)
  moon <- suncalc::getMoonTimes(data = data)
  i = 1
  for(i in 1:length(bats)){
    bdf <- df[df$id == bats[i],]
    dates <- unique(date(ymd_hms(bdf$timestamp)))
    # date(ymd_hms(bdf$timestamp)) %>% table
    j = 3
    for(j in 1:length(dates)){
      d <- bdf[which(date(ymd_hms(bdf$timestamp)) == dates[j]),]

      if(gps == TRUE){
        # plot tracks
        lon <- d$location_long
        lat <- d$location_lat
        gps_time <- d$timestamp
        # get sunrise and sunset times
        layout(rbind(c(1, 2), c(1, 3)))
        par(mar = c(0, 4, 0, 0), oma = c(4, 0, 4, 4), xpd = NA)
        plot(lon, lat, asp = 1, type = "o", pch = 19)
        plot(gps_time, lon, type = "o", pch = 19,
             xlab = "", xaxt = "none",
             xlim = c(gps_time[1]-3*3600,
                      gps_time[length(gps_time)]+3*3600))
        abline(v = sun$sunset, col = "blue", xpd = FALSE)
        abline(v = sun$sunrise, col = "orange", xpd = FALSE)
        plot(gps_time, lat, type = "o", pch = 19,
             xlim = c(gps_time[1]-3*3600,
                      gps_time[length(gps_time)]+3*3600))
        abline(v = sun$sunset, col = "blue", xpd = FALSE)
        abline(v = sun$sunrise, col = "orange", xpd = FALSE)
        legend("bottomright", legend = c("sunset", "sunrise"), col = c("blue", "orange"), lty = 1)
        # get speed

        # label commutes
        ## first passage time?
      }

      x <- {}
      y <- {}
      z <- {}
      time <- {}

      if(Burst == TRUE){
        burst <- {}

        if(nrow(d) > 100){
          db <- foreach(ii = 1:nrow(d), .combine = rbind) %dopar% {
            require(magrittr)
            temp <- {}
            x <- NA
            y <- NA
            z <- NA
            time <- NA
            burst <- NA
            if(tag_type == "Wildfi"){
              temp <- d$accInGBurst[ii] %>% strsplit(" ") %>%
                unlist %>%
                as.numeric
            }
            if(tag_type == "eObs"){
              temp <- d$eobs_accelerations_raw[ii] %>% strsplit(" ") %>%
                unlist %>%
                as.numeric
            }

            if(length(temp) > 0){
              x <- temp[seq(1, length(temp), 3)]
              y <- temp[seq(2, length(temp), 3)]
              z <- temp[seq(3, length(temp), 3)]
              temp_time <- seq.POSIXt(from = d$timestamp[ii],
                                      to = d$timestamp[ii]+length(temp)/sampling_rate/3,
                                      length.out = length(temp)/3)
              time <- format(temp_time, "%Y-%m-%d %H:%M:%OS3") %>% as.character
              burst <- rep(ii, length(time))
            }
            temp <- data.frame(x,y,z,time,burst)
          }
          db <- na.omit(db)
          db$time <- ymd_hms(db$time)
          x <- db$x
          y <- db$y
          z <- db$z
          time <- db$time
          burst <- db$burst
        }
      }
      if(Burst == FALSE){
        x <- df$ACCX
        y <- df$ACCY
        z <- df$ACCZ
        time <- df$time
      }

      # pca
      if(PCA == TRUE){
        pc <- prcomp(cbind(x, y, z), scale. = FALSE)
        w <- Wave(left = pc$x[,1], samp.rate = sampling_rate, bit = 16)
        wf <- ffilter(w, f= sampling_rate, from = min_freq, to = max_freq, bandpass = TRUE)
        spectro(wf, f = sampling_rate, wl = 1024*4, ovlp = 50, fastdisp = TRUE)
        # add sunrise and sunset?
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

      # get periods where wing beat can be estimated
      idx <- which(!is.na(pf[,2])) -1
      stable_var <- seqToIntervals(idx) %>% as.data.frame
      stable_var$diff <- NA

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

      stable_var$length <- stable_var$to - stable_var$from# + 1

      ## create
      # multiply by sampling interval of pf
      w_stable <- stable_var * (pf[2,1])

      df_stable <- stable_var * (pf[2,1]) * sampling_rate
      ## fix last value
      # df_stable$to[nrow(df_stable)] <- nrow(db)
      # nrow(db)-6233 * (pf[2,1]) * sampling_rate

      # fix first value
      if(w_stable[1,1] == 0) w_stable[1,1] <- 1
      if(df_stable[1,1] == 0) df_stable[1,1] <- 1

      # get FFT peak freq
      if(use_FFT == TRUE){
        try({
          w_stable$peak_freq <- NA
          # w_stable$peak_amp <- NA
          w_stable$amplitude <- NA
          # jj = which.max(w_stable$length)
          for(jj in 1:nrow(w_stable)){
            if(stable_var$length[jj] > 1){
              # filter short periods
              if(difftime(db$time[df_stable$to[jj]],db$time[df_stable$from[jj]], units = "mins") > min_seg_duration){
                w_cut <- Wave(left = db$z[df_stable$from[jj]:df_stable$to[jj]] %>% na.omit,
                              samp.rate = sampling_rate, bit = 16)

                plot(w)
                abline(v = w_stable$to[jj], col = 2)
                abline(v = w_stable$from[jj], col = 3)
                wf_cut <- ffilter(w_cut, f= sampling_rate, from = min_freq, to = max_freq, bandpass = TRUE)
                spectro(wf_cut, f = sampling_rate)
                spec <- meanspec(wf_cut, f=sampling_rate, plot = FALSE)
                layout(1)
                peak <- fpeaks(spec, nmax = 8)
                pidx <- which(peak[,1] > 0.0025)
                midx <- which.max(peak[pidx,2])
                w_stable$peak_freq[jj] <- peak[pidx[midx],1]
                # w_stable$peak_amp[jj] <- peak[pidx[midx],2] # doesn't seem like a useful value
                # df$ACCZ[df_stable$from[i]:df_stable$to[i]] %>% hist
                maxz <- roll_max(db$z[df_stable$from[jj]:df_stable$to[jj]], width = 6) %>%
                  median(na.rm = TRUE)
                minz <- roll_min(db$z[df_stable$from[jj]:df_stable$to[jj]], width = 6) %>%
                  median(na.rm = TRUE)
                # abline(v = c(minz, maxz), col = 2, lty = 2)
                w_stable$amplitude[jj] <- maxz - minz
              }
            }
          }
        })
      }

      # wavelet give frequencies lower than expected
      if(use_wavelet == TRUE){
        try({
          w_stable$peak_freq <- NA
          # w_stable$peak_amp <- NA
          w_stable$amplitude <- NA
          jj = which.max(w_stable$length)
          for(jj in 1:nrow(w_stable)){
            if(stable_var$length[jj] > 1){
              # filter short periods
              if(difftime(db$time[df_stable$to[jj]],db$time[df_stable$from[jj]], units = "mins") > min_seg_duration){
                w_cut <- Wave(left = db$z[df_stable$from[jj]:df_stable$to[jj]] %>% na.omit,
                              samp.rate = sampling_rate, bit = 16)
                spectro(w_cut)

                wvl <- dplR::morlet(y1 = db$z[df_stable$from[jj]:df_stable$to[jj]][20001:40000],
                                    x1 = db$time[df_stable$from[jj]:df_stable$to[jj]][1:20000] %>% as.numeric,
                                    p2 = 3, dj = 0.01, siglvl = 0.99)
                layout(1)
                # wvl$Power %>% image
                plot(wvl$x, wvl$Scale[apply(wvl$Power,1,which.max)],
                     cex = apply(wvl$Power,1,max)/max(apply(wvl$Power,1,max)),
                     ylim = c(0,8), col = rgb(0,0,0,.1))
                lines(wvl$x, roll_mean(wvl$Scale[apply(wvl$Power,1,which.max)], width = 100), col = 2)
                abline(h = 6, col = 2)
                (wvl$Scale[apply(wvl$Power,1,which.max)]) %>% hist(breaks = 1000)
                dplR::wavelet.plot(wvl, useRaster = NA)
              }
            }
          }
        })
      }


      try({
        png(filename = paste0(path, "wingbeat/", studies[k], bats[j], "_", dates[i],"_domfreq_segments.png"))
        plot(x = ymd_hms(db$time[df_stable$from])+time_offset, y = w_stable$peak_freq*1000,
             cex = 2*scales::rescale(w_stable$amplitude),
             xlab = "local time", ylab = "dominant wingbeat frequency (Hz)")
        segments(x0 = ymd_hms(db$time[df_stable$from])+time_offset,
                 y0 = w_stable$peak_freq*1000,
                 x1 = ymd_hms(db$time[df_stable$to])+time_offset,
                 y1 = w_stable$peak_freq*1000)
        dev.off()
      })

      save(db,w,wf,pf, sampling_rate, stable_var, w_stable, df_stable,
           file = paste0(path, "wingbeat/", studies[k], bats[j], "_", dates[i],"_wingbeatfreq.robj"))
  }
  }
}
