fft_peak_freq <- function(data, time, lowpass = 5, highpass = 15, sampling_rate = 40){
  require(tuneR)
  require(seewave)
  # calculate fft of data
  test <- fft(data)

  # extract magnitudes and phases
  magn <- Mod(test) # sqrt(Re(test)*Re(test)+Im(test)*Im(test))
  # phase <- Arg(test) # atan(Im(test)/Re(test))
  # plot(phase, type = "l")

  # select only first half of vectors
  magn.1 <- magn[1:(length(magn)/2)]
  #phase.1 <- Arg(test)[1:(length(test)/2)]

  # plot(magn.1,type="l")

  # generate x-axis with frequencies
  x.axis <- 1:length(magn.1)/time
  res <- x.axis[1]
  idx <- which(x.axis > highpass | x.axis < lowpass)
  magn.1[idx] <- 0

  plot(x=x.axis,y=magn.1,type="l")
  idx <- which.max(magn.1)
  freq <- x.axis[idx]
  amp <- max(magn.1)
  abline(v = freq, col = 2)

  w <- Wave(left = data, samp.rate = sampling_rate, bit = 16)
  wfreq <- NA
  wamp <- NA
  # plot(w)
  try({
    spec <- meanspec(w, f=sampling_rate, wl = round(length(w)/5,-1), plot = FALSE)
    peak <- fpeaks(spec, nmax = 4, plot = FALSE)
    pidx <- which(peak[,1] > 0.0005)
    midx <- which.max(peak[,2])
    if(length(midx) > 0){
      wfreq <- peak[pidx[midx],1] * 1000
      wamp <- peak[pidx[midx],2]
    }
  })
  rms <- NA
  rms <- rms(w@left*100)
  wf <- ffilter(w, f= sampling_rate, from = 0, to = 1,
                bandpass = TRUE,
                wl = round(length(w)/5,-1))
  # plot(wf)
  rms_filter <- NA
  rms_filter <- rms(wf*100)

  return(data.frame(freq, amp, res, wfreq, wamp, rms, rms_filter))
}
