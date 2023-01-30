## function to get the dominant frequency from ACC data during commutes

# calculate sunset and sunrise
# add wavelet
Burst = TRUE
PCA = TRUE
sampling_rate = 18.74
max_freq = 4
dominant_freq <- function(df, Burst = FALSE, PCA = FALSE, sampling_rate = 25,
                          min_freq = 2, max_freq = 8,
                          wavelet = FALSE, gps = FALSE, saved_cores = 4){
  # load libraries
  require(lubridate)
  require(foreach)
  require(doParallel)
  require(suncalc)
  require(marcher)

  source("../NoctuleMigration/scr/scan.track.map.r")

  # get timestamp with milliseconds
  op <- options(digits.secs=3)
  # set up parallelization
  cores=detectCores()
  cl <- makeCluster(cores[1]-saved_cores) #not to overload your computer
  registerDoParallel(cl)

  bats <- unique(df$tag_local_identifier)
  d <- unique(date(df$timestamp))
  data <- data.frame(date = d, lat = median(df$location_lat, na.rm = TRUE),
                     lon = median(df$location_long, na.rm = TRUE))
  sun <- suncalc::getSunlightTimes(data = data)
  moon <- suncalc::getMoonTimes(data = data)
  i = 1
  for(i in 1:length(bats)){
    bdf <- df[df$tag_local_identifier == bats[i],]
    dates <- unique(date(ymd_hms(bdf$timestamp)))
    j = 1
    for(j in 1:length(dates)){
      d <- bdf[which(date(ymd_hms(bdf$timestamp)) == dates[j]),]

      if(gps == TRUE){
        # plot tracks
        x <- d$location_long
        y <- d$location_lat
        time <- d$timestamp
        # get sunrise and sunset times
        layout(rbind(c(1, 2), c(1, 3)))
        par(mar = c(0, 4, 0, 0), oma = c(4, 0, 4, 4), xpd = NA)
        plot(x, y, asp = 1, type = "o", pch = 19)
        plot(time, x, type = "o", pch = 19,
             xlim = c(time[1]-3*3600,
                      time[length(time)]+3*3600))
        abline(v = sun$sunset, col = "blue", xpd = FALSE)
        abline(v = sun$sunrise, col = "orange", xpd = FALSE)
        plot(time, y, type = "o", pch = 19,
             xlim = c(time[1]-3*3600,
                      time[length(time)]+3*3600))
        abline(v = sun$sunset, col = "blue", xpd = FALSE)
        abline(v = sun$sunrise, col = "orange", xpd = FALSE)
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
            temp <- d$eobs_accelerations_raw[ii] %>% strsplit(" ") %>%
              unlist %>%
              as.numeric
            if(length(temp) > 0){
              x <- temp[seq(1, length(temp), 3)]
              y <- temp[seq(2, length(temp), 3)]
              z <- temp[seq(3, length(temp), 3)]
              temp_time <- seq.POSIXt(from = d$timestamp[ii],
                                      to = d$timestamp[ii]+length(temp)/sampling_rate/3,
                                      length.out = length(temp)/3)
              time <- format(temp_time, "%Y-%m-%d %H:%M:%OS3") %>% as.character
              burst <- rep(ii, length(temp))
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
      }
      if(PCA == FALSE){
        # create wave form
        w <- Wave(left = z, samp.rate = sampling_rate, bit = 16)
        wf <- ffilter(w, f= sampling_rate, from = min_freq, to = max_freq, bandpass = TRUE)
        spectro(wf, f = sampling_rate, wl = 1024*4, ovlp = 50, fastdisp = TRUE)
      }
      par(new=TRUE)
      pf <- dfreq(wf, f=sampling_rate, wl = 1024*4, ovlp=50, threshold=30,
                  type="l", col="red", lwd=0.5, xlab = "", ylab = "")

      # get periods where wing beat can be estimated
      idx <- which(!is.na(pf[,2]))
      stable_var <- seqToIntervals(idx) %>% as.data.frame
      stable_var$length <- stable_var$to - stable_var$from + 1

      ## create
      # multiply by sampling interval of pf
      w_stable <- stable_var * (pf[2,1])
      df_stable <- stable_var * (pf[2,1]) * sampling_rate

      try({
        w_stable$peak_freq <- NA
        w_stable$peak_amp <- NA
        w_stable$amplitude <- NA
        jj = 1
        for(jj in 1:nrow(w_stable)){
          if(stable_var$length[jj] > 1){
            w_cut <- Wave(left = db$z[df_stable$from[jj]:df_stable$to[jj]] %>% na.omit,
                          samp.rate = sampling_rate, bit = 16)
            plot(w_cut)
            wf_cut <- ffilter(w_cut, f= sampling_rate, from = 2, to = 4, bandpass = TRUE)
            spectro(wf_cut, f = sampling_rate)
            spec <- meanspec(wf_cut, f=sampling_rate, plot = TRUE)
            peak <- fpeaks(spec, nmax = 8)
            pidx <- which(peak[,1] > 0.0025)
            midx <- which.max(peak[pidx,2])
            w_stable$peak_freq[jj] <- peak[pidx[midx],1]
            w_stable$peak_amp[jj] <- peak[pidx[midx],2]
            # df$ACCZ[df_stable$from[i]:df_stable$to[i]] %>% hist
            maxz <- roll_max(db$z[df_stable$from[jj]:df_stable$to[jj]], width = 6) %>%
              median(na.rm = TRUE)
            minz <- roll_min(db$z[df_stable$from[jj]:df_stable$to[jj]], width = 6) %>%
              median(na.rm = TRUE)
            # abline(v = c(minz, maxz), col = 2, lty = 2)
            w_stable$amplitude[jj] <- maxz - minz
          }
        }
      })
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