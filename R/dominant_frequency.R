## function to get the dominant frequency from ACC data during commutes

dominant_freq <- function(df, Burst = FALSE, PCA = TRUE){
  require(lubridate)
  
  dates <- unique()
  for(j in 1:length(dates)){
    d <- df[which(date(ymd_hms(df$timestamp)) == dates[j]),]
    x <- {}
    y <- {}
    z <- {}
    time <- {}
    
    if(Burst == TRUE){
      burst <- {}
      
      if(nrow(d) > 1000){
        b <- 1
        for(i in 1:nrow(d)){
          temp <- {}
          temp <- d$eobs_accelerations_raw[i] %>% strsplit(" ") %>% 
            unlist %>% 
            as.numeric
          # d$resting %>% nchar
          # if(length(temp) > 0 & nchar(d$resting[i]) == 0 & nchar(d$foraging[i]) == 0){
          x <- c(x, temp[seq(1, length(temp), 3)])
          y <- c(y, temp[seq(2, length(temp), 3)])
          z <- c(z, temp[seq(3, length(temp), 3)])
          temp_time <- seq.POSIXt(from = d$timestamp[i], 
                                  to = d$timestamp[i]+length(temp)/25/3, 
                                  length.out = length(temp)/3)
          time <- c(time, format(temp_time, "%Y-%m-%d %H:%M:%OS2") %>% as.character)
          burst <- c(burst, rep(b, length(temp)))
        }
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
      wpc <- Wave(left = pc$x[,1], samp.rate = 25, bit = 16)
      wfpc <- ffilter(wpc, f= 25, from = 5, to = 7, bandpass = TRUE)
      spectro(wfpc, f = 25, wl = 1024*2)
    }
    
  # create wave form
    w <- Wave(left = z, samp.rate = 25, bit = 16)
    wf <- ffilter(w, f= 25, from = 5, to = 7, bandpass = TRUE)
    spectro(wf, f = 25, wl = 1024*4, ovlp = 50, fastdisp = TRUE)
  }
  
  
  
  
  
  
}