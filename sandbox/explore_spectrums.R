lowpass <- 5
highpass <- 15

x <- b$z0
wingbeat <- function(x, sampling_rate, highpass, lowpass){
  s <- spectrum(x, method = "pgram", plot = T)
  freq_bins <- s$freq*sampling_rate
  res <- freq_bins[2]
  idx <- which(freq_bins < lowpass | freq_bins > highpass)
  s$spec[idx] <- 0
  freq <- s$freq[which.max(s$spec)]*(sampling_rate)
  spec <- max(s$spec)
}



p@starts
p@spec[[2]] %>% plot
p@energy
p@variance
p@samp.rate

require(lomb)
t <- ts(x)
per <- lsp(t)
res <- per$scanned[1]*sampling_rate
freq_bins <- per$scanned[1:(length(per$scanned)-1)]*sampling_rate
idx <- which(freq_bins < lowpass | freq_bins > highpass)
per$power[idx] <- 0
freq <- per$scanned[which.max(per$power)]*(sampling_rate)

per$power
p <- getpeaks(per, plotit = FALSE, npeaks = 5)
freq <- p$data$time[1] * sampling_rate
npower <- p$data$peaks[1]



fft_peak_freq <- function(data, time){

  # calculate fft of data
  test <- fft(data)

  # extract magnitudes and phases
  magn <- Mod(test) # sqrt(Re(test)*Re(test)+Im(test)*Im(test))
  # phase <- Arg(test) # atan(Im(test)/Re(test))
  # plot(phase, type = "l")

  # select only first half of vectors
  magn.1 <- magn[1:(length(magn)/2)]
  #phase.1 <- Arg(test)[1:(length(test)/2)]

  plot(magn.1,type="l")

  # generate x-axis with frequencies
  x.axis <- 1:length(magn.1)/time

  plot(x=x.axis,y=magn.1,type="l")
  idx <- which.max(magn.1)
  peak_freq <- x.axis[idx]
  peak_amp <- max(magn.1)
  abline(v = peak_freq, col = 2)

  return(data.frame(peak_freq, peak_amp))
}


## fft

plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))

  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2]

  plot(plot.data, t="h", lwd=2, main="",
       xlab="Frequency (Hz)", ylab="Strength",
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

trajectory <- sapply(ts, function(t) f(t,w))
X.k <- fft(b$pc)                   # find all harmonics with fft()
plot.frequency.spectrum(X.k, xlimits=c(0,nrow(b)))


library(Rmpfr) # to avoid an integer overflow problem in `WHICH.MOD`
d <- w
s.rate <- 40
microbenchmark(
  WHICH.MOD = which((mfft<-Mod(fft(d@left)))[1:(length(d@left)/2)] == max(abs(mfft[1:(length(d@left)/2)]))) * mpfr(s.rate,100) / length(d@left),
  MEANSPEC = (ms<-meanspec(d,f=s.rate,wl=sampling_rate*8,plot=F))[which.max(ms[,2]),1]*1000,
  #DFREQ.HIST = (h<-hist(dfreq(d,f=s.rate,wl=sampling_rate*8,plot=F)[,2],200,plot=F))$mids[which.max(h$density)]*1000,
  #DFREQ.DENS = (dens <- density(dfreq(d,f=s.rate,wl=sampling_rate*8,plot=F)[,2],na.rm=T))$x[which.max(dens$y)]*1000,
  FPEAKS.MSPEC = fpeaks(meanspec(d,f=s.rate,wl=sampling_rate*8,plot=F),nmax=1,plot=F)[,1]*1000 ,
  times=25)

library(TSA)


which(Mod(fft(b$pc)) == max(abs(Mod(fft(b$pc))))) * sampling_rate / length(b$pc)
# Frequency_of_Peak = Data_Sample_Rate * Bin_number_of_Peak / Length_of_FFT
fft_result <- Mod(fft(b$pc))
fft_result <- fft_result[1:(nrow(b)%/%2)]
plot(fft_result)

n <- length(fft_result)
freq_bins <- (1:n) * (1 / (n * sampling_rate))

# Apply FFT to the data
fft_result <- fft(b$pc)/nrow(b)

# Calculate frequency bins
n <- length(fft_result)
freq_bins <- (1:n) * (1 / (n * sampling_rate))
freq_bins <- seq(0, sampling_rate/2, by = sampling_rate/2/nrow(b))

# Filter frequencies below 0.002
fft_result[freq_bins < lowpass_frequency | freq_bins > highpass_frequency] <- 0

# Calculate the magnitudes of amplitudes
index_of_amplitudes <- Mod(fft_result)

# # Plot the frequency spectrum
plot(freq_bins,
     index_of_amplitudes,
     type = "l", xlab = "Frequency (Hz)", ylab = "Amplitude")

# Identify dominant frequency
index_of_max_amplitude <- which.max(index_of_amplitudes)
dominant_frequency <- freq_bins[index_of_max_amplitude]

df <- data.frame(burst = i, burst_start = idx[i], dominant_frequency,
                 max_amplitude = max(index_of_amplitudes[2:(n/2)]))
cat("Dominant frequency:", round(dominant_frequency * 1000, 2), "Hz\n")
lcb <- rbind(lcb, df)
}

# lcb$id <- "29507A30"
lcb$sampling_rate <- sampling_rate
ggplot(lcb, aes(x = burst, y = dominant_frequency*1000, size = max_amplitude))+
  geom_point()

b$z0 <- mean(b$z) - b$z
w <- Wave(left = b$z0, samp.rate = sampling_rate, bit = 16)
p <- periodogram(w, width = sampling_rate*4)
plot(p, which = 1)
plot(p, which = 2)
summary(p)
# image(p)
p@spec
abline(v = FF(p, peakheight = 0.4), col = 2)
FF(p, peakheight = 0.4)

p1 <- TSA::periodogram(y = b$pc)
p1$freq[which.max(p1$spec)]*(sampling_rate)

s1 <- spectrum(rep(b$z0,1))
freq_bins <- s1$freq*sampling_rate
idx <- which(freq_bins < lowpass_frequency | freq_bins > highpass_frequency)
s1$spec[idx] <- 0
s1$freq[which.max(s1$spec)]*(sampling_rate)
# s1$freq[1]*sampling_rate

spec <- meanspec(w, wl = sampling_rate*4)
flapping <- ffilter(w, f= sampling_rate, from = 5, to = 15,
                    bandpass = TRUE, wl = sampling_rate*4)
# spectro(flapping, f = sampling_rate*4)
spec <- meanspec(flapping, f=sampling_rate, wl = sampling_rate*8,
                 plot = FALSE)
spec[2,1]*1000
peak <- fpeaks(spec, nmax = 20, plot = TRUE)
peak_order <- peak[order(peak[,2], decreasing = TRUE),]