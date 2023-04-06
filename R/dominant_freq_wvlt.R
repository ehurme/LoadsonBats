p_load(janitor, rnaturalearth)

countries <- ne_download(scale = 10, type = "countries")

df <- fread(files[3]) %>% clean_names()

save_path = "./../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/Wingbeats/PCA/"
tag_id <- "2268"
PCA = TRUE
sampling_rate = 18.74
min_freq = 2
max_freq = 8
wavelet = TRUE
saved_cores = 8
gps = TRUE
FFT = TRUE
tag_type = "eObs" # "wildfi"
firetail = TRUE
firetail_filter = TRUE
FFT_amp_thresh = -0.2
solar_time = TRUE
save_files = TRUE
levels = 1000


dominant_freq_wvlt <- function(df,
                          save_path = NULL,
                          tag_id = NULL,
                          #Burst = TRUE,
                          PCA = TRUE,
                          sampling_rate = 25, # P. hastatus
                          min_freq = 2,
                          max_freq = 8,
                          wavelet = TRUE,
                          levels = 1000,
                          saved_cores = 4,
                          gps = TRUE,
                          FFT = TRUE,
                          FFT_amp_thresh = -0.2,
                          tag_type = "eObs", # "wildfi"
                          firetail = TRUE,
                          firetail_filter = FALSE,
                          solar_time = TRUE,
                          save_files = TRUE){
  # load libraries
  ## utilities
  require(pacman)
  p_load(data.table, magrittr, dplyr,
         foreach, doParallel, R.utils, roll,
         lubridate, suncalc,
         seewave, Rwave, tuneR,
         mclust)

  # get timestamp with milliseconds
  op <- options(digits.secs=3)
  # set up parallelization
  cores=detectCores()
  cl <- makeCluster(cores[1]-saved_cores) #not to overload your computer
  registerDoParallel(cl)

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
    try(data <- data.frame(date = dates, lat = median(df$location_lat, na.rm = TRUE),
                           lon = median(df$location_long, na.rm = TRUE)))
    sun <- suncalc::getSunlightTimes(data = data)
    moon <- suncalc::getMoonTimes(data = data)
    df$solar_time <- df$timestamp - hour(sun$solarNoon[1])*60*60
  }

  i = 1
  bat_freq <- data.frame()
  for(i in 1:length(bats)){
    full_path <- paste0(save_path, tag_id, "/")
    if(!dir.exists(full_path)){
      dir.create(full_path)
    }

    bdf <- df[df$id == bats[i],]
    dates <- unique(date(bdf$solar_time))

    date_freq <- data.frame()

    j = 1
    for(j in 3:length(dates)){
      d <- bdf[which(date(ymd_hms(bdf$solar_time)) == dates[j]),]

      if(gps == TRUE){
        # plot tracks
        lon <- d$location_long
        lat <- d$location_lat
        gps_time <- d$timestamp
        idx <- which(!is.na(lon))
        # get sunrise and sunset times
        if(any(!is.na(lon))){
          try({
            png(file = paste0(full_path, "GPS_", bats[i],"_", dates[j], ".png"), width = 1000, height = 800)
            layout(rbind(c(1, 2), c(1, 3)))
            par(mar = c(0, 4, 0, 0), oma = c(4, 0, 4, 4), xpd = NA)
            plot(lon[idx], lat[idx], asp = 1, type = "o", pch = 19, cex = 0.5,
                 xlab = "Longitude", ylab = "Latitude")
            #lines(lon[idx], lat[idx])
            # lines(countries)
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

            lines(gps_time[idx], lon[idx])
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
            lines(gps_time[idx], lat[idx])
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
      burst <- {}

      if(nrow(d) > 100){
        ii = 1
        temp_df <- data.frame()
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
          commuting_idx <- which(d$annotation_layer_commuting != "" & d$type == "acc")
          foraging_idx <- which(d$annotation_layer_foraging != "" & d$type == "acc")
          resting_idx <- which(d$annotation_layer_resting != "" & d$type == "acc")

          db$firetail <- NA
          db$firetail[db$burst %in% commuting_idx] <- "c"
          db$firetail[db$burst %in% foraging_idx] <- "f"
          db$firetail[db$burst %in% resting_idx] <- "r"
          behav <- db$firetail
          if(firetail_filter == TRUE){
            fire_idx <- which(behav != "r")
            x <- x[fire_idx]
            y <- y[fire_idx]
            z <- z[fire_idx]
            time <- time[fire_idx]
          }
        }

        # which bursts to use?
        bursts <- {}
        # table(db$firetail)
        if(length(unique(db$firetail)) <= 1){Freqs <- data.frame()}
        if(length(unique(db$firetail)) > 1){
          bursts <- unique(db$burst)
          if(firetail_filter == "c"){
            bursts <- unique(db$burst[which(db$firetail == "c")])
          }
          if(firetail_filter == "f"){
            bursts <- unique(db$burst[which(db$firetail == "f")])
          }
          if(firetail_filter == "cf"){
            bursts <- unique(db$burst[which(db$firetail == "c" |
                                              db$firetail == "f")])
          }
          B = bursts[1]
          Freqs <- data.frame()
          if(length(bursts) > 100) bursts <- bursts[1:10]
          for(B in bursts){
            b <- db[db$burst == B,]
            C <- wavelet_domfreq(b, levels,
                                 sampling_rate, min_freq, max_freq,
                                 PCA, FFT)
            peaks <- C[[2]]
            while(peaks[2,3] > FFT_amp_thresh){
              peak <- matrix(c(freq = mean(peaks[1:2,1]),
                               amp = mean(peaks[1:2,2]),
                               diff = NA), nrow = 1)
              peaks <- peaks[2:nrow(peaks),1:3]
              peaks[1,] <- peak
              peaks[,3] <- c(NA, diff(peaks[,2]))
              if(nrow(peaks) < 3) break
            }

            c <- C[[1]]
            # hist(c$amp[c$amp > 0.1])
            # with(c[c$amp > 0.1,], plot(time, freqs, cex = amp/max(amp)))
            c0 <- na.omit(c)

            # abline(h = c0$freqs[which.max(c0$amp)])
            # abline(h = peaks[1,1]*1000, lty = 2)
            freqs <- {}
            if(nrow(c0) == 0){
              freqs <- data.frame(burst = B, time = b$time[1],
                                  wvlt = NA,
                                  wvlt_amp = NA,
                                  fft = peaks[1,1]*1000,
                                  fft_amp = peaks[1,2],
                                  fft_amp_diff = peaks[2,3])
            }
            if(nrow(c0) > 0){
              freqs <- data.frame(burst = B, time = b$time[1],
                                  wvlt = c0$freqs[which.max(c0$amp)],
                                  wvlt_amp = c0$amp[which.max(c0$amp)],
                                  fft = peaks[1,1]*1000,
                                  fft_amp = peaks[1,2],
                                  fft_amp_diff = peaks[2,3])
            }
            Freqs <- rbind(Freqs , freqs)
          }

          png(file = paste0(full_path, "Freqs_", bats[i],"_", dates[j], ".png"), width = 1000, height = 800)
          plot(Freqs$time, Freqs$fft, cex = Freqs$fft_amp/max(Freqs$fft_amp),
               col = 5, pch = 16, type = "o",
               # ylim = c(min(c(Freqs$fft, Freqs$wvlt)),
               #          max(c(Freqs$fft, Freqs$wvlt))),
               ylim = c(min_freq, max_freq),
               ylab = "Wingbeat Frequency (Hz)",
               xlab = "Time (UTC)",
               main = paste0("Bat: ", bats[i],", Date: ", dates[j]))
          points(Freqs$time, Freqs$wvlt, cex = Freqs$wvlt_amp/max(Freqs$wvlt_amp),
                 col = 4, pch = 16, type = "o")
          abline(v = sun$sunset, col = "blue", xpd = FALSE)
          abline(v = sun$sunrise, col = "orange", xpd = FALSE)
          legend("topleft", legend = c("Commuting","Foraging","Resting"),
                 col = 3:1, pch = 16, pt.cex = (3:1)/6, cex = 0.8)
          legend("topright", legend = c("Wavelet", "FFT", "Sunset", "Sunrise"),
                 pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 1),
                 col = c(5:4, "blue", "orange"), cex = 0.8)
          fire <- db$firetail |> factor(levels = c("r", "f", "c"))
          points(db$time, rep(min_freq, nrow(db)), cex = as.numeric(fire)/10, col = fire)
          dev.off()
        }

      }
      date_freq <- rbind(date_freq, Freqs)
    }
    date_freq$bat <- bats[i]
    bat_freq <- rbind(bat_freq, date_freq)
  }
  if(save_files == TRUE){
    save(bat_freq, save_path,
         tag_id,
         PCA,
         sampling_rate,
         min_freq,
         max_freq,
         wavelet,
         saved_cores,
         gps,
         FFT,
         tag_type,
         firetail,
         firetail_filter,
         FFT_amp_thresh,
         solar_time,
         save_files,
         levels,
         file = paste0(full_path, "Freqs_", bats[i], ".robj"))
  }
  return(bat_freq)
}
