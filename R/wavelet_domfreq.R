db$burst[which(db$behavior == "commuting")] |> unique()
# db$burst[which(db$behavior == "foraging")] |> unique()
temp <- db[db$burst == 586,]
#
# min_freq <- 2
# max_freq <- 4
# levels <- 1000

C <- wavelet_domfreq(temp, levels = 1000, sampling_rate = 18.74, min_freq = 2, max_freq = 8, PCA = FALSE, FFT = TRUE)
peaks <- C[[2]]
c <- C[[1]]
hist(c$amp[c$amp > 0.1])
with(c[c$amp > 0.1,], plot(time, freqs, cex = amp/max(amp)))
c0 <- na.omit(c)
abline(h = c0$freqs[which.max(c0$amp)])
abline(h = peaks[1,1]*1000, lty = 2)

wavelet_domfreq <- function(temp, levels, sampling_rate, min_freq, max_freq, PCA = FALSE, FFT = FALSE){
  require(WaveletComp)
  pl.data <- data.frame(x = temp$z, date = temp$time) |> na.omit()
  if(PCA){
    pl.pc <- prcomp(with(temp, cbind(x, y, z)), scale. = FALSE)
    pl.data <- data.frame(x = pl.pc$x[,1], date = temp$time)
  }

  if(FFT){
    pl.w <- Wave(left = rep(pl.data$x,3), samp.rate = sampling_rate, bit = 16)
    plot(pl.w)
    # try(pl.w <- ffilter(pl.w, f = sampling_rate, from = 2, to = 4, bandpass = TRUE))
    spectro(pl.w, f = sampling_rate, wl = 16, ovlp = 75, fastdisp = TRUE)
    dev.off()
    spec <- meanspec(pl.w, f=sampling_rate, plot = TRUE)
    peaks <- fpeaks(spec)
    peaks <- unlist(peaks[order(peaks[,2], decreasing = TRUE),])
    #peaks <- cbind(peaks, diff = c(NA, diff(peaks[,2])))
    peaks <- peaks[peaks[,1] > 0.002 & peaks[,1] < 0.004,]
    peaks <- cbind(peaks, diff = c(NA, diff(peaks[,2])))
    peaks
  }

  pl.wl <- analyze.wavelet(pl.data, "x",
                           loess.span = 0,
                           dt = 1/sampling_rate,
                           dj = 1/levels,
                           lowerPeriod = min_freq,
                           upperPeriod = max_freq,
                           make.pval = TRUE,
                           n.sim = 10)
  wt.image(pl.wl, color.key = "quantile", n.levels = levels, col.ridge = "purple",
           legend.params = list(lab = "wavelet power levels", mar = 4.7))
  my_freq <- pl.wl$Period[apply(pl.wl$Power,2,which.max)]
  my_amp <- pl.wl$Ampl[apply(pl.wl$Power,2,which.max)]

  my_freq[which(my_freq == min_freq | my_freq == max_freq)] <- NA
  #plot(my_freq)

  # remove points outside the coi
  my_coi <- 2^pl.wl$coi.2[3:(length(pl.wl$coi.2)-2)]
  idx <- which(my_freq > my_coi)
  my_freq[idx] <- NA

  # ridge <- pl.wl$Ridge
  # rownames(ridge) <- pl.wl$axis.2
  # colnames(ridge) <- pl.wl$axis.1
  # ridge[which(ridge == 0)] <- NA
  # table(ridge)

  # image(t(ridge), col = 1)
  # plot(pl.wl$axis.1,
  #      my_freq,
  #      ylim = c(min_freq, max_freq),
  #      #pl.wl$Period[apply(pl.wl$Ridge,2,which.max)],
  #      #type = "l")
  # pch = 16, col = rgb(0,0,0,.1))
  # lines(pl.wl$coi.1, 2^pl.wl$coi.2, col = 2)

# hist(my.wC$Period[apply(my.wC$Ridge,2,which.max)][idx])
# summary(my.wC$Period[apply(my.wC$Ridge,2,which.max)][idx])
  summary(my_freq)
  c <- data.frame(freqs = my_freq, amp = my_amp, time = pl.wl$axis.1)# |> na.omit()
  # plot(c$time, c$freqs, cex = c$amp)
  ifelse(FFT, return(list(c, peaks)), return(list(c)))
}

