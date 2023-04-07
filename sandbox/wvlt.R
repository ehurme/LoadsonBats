# Load necessary libraries
library(Rwave)
library(seewave)
library(tuneR)

# Perform wavelet transform
wt <- cwt(my.data, dj=0.125, dj.t=0.125, wavelet="morlet", parallel=TRUE)

# Extract dominant frequency and amplitude
freq <- rowMeans(wt$freq[apply(abs(wt$CWT), 2, which.max)])
amp <- rowMeans(abs(wt$CWT))

# Plot the results
plot(freq, type="l", ylab="Frequency (Hz)", xlab="Time (s)")
plot(amp, type="l", ylab="Amplitude", xlab="Time (s)")

x <- 1:512
chirp <- sin(2*pi * (x + 0.002 * (x-256)^2 ) / 16)
plot.ts(chirp)
retChirp <- cwt(chirp, noctave=5, nvoice=12)

x <- 1:512
chirp <- sin(2*pi * (x + 0.002 * (x-256)^2 ) / 16)
chirp <- chirp + 1i * sin(2*pi * (x + 0.004 * (x-256)^2 ) / 16)
retChirp <- cwtp(chirp, noctave=5, nvoice=12)

x <- 1:512
chirp <- sin(2*pi * (x + 0.002 * (x-256)^2 ) / 16)
retChirp <- cwtTh(chirp, noctave=5, nvoice=12, moments=20)

m1 = gabor(1024, 512,  2 * pi, 20 )
plot.ts(m1)
w1 <- Wave(left = Re(m1), samp.rate = 10, bit = 16)
layout(1)
plot(w1)
spectro(w1)
retm <- cwt(m1, noctave=5, nvoice=12)
