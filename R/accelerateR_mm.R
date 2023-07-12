library(pacman)
p_load(data.table, tidyverse, magrittr, dplyr, accelerateR, stringr)

mm <- fread("../../../Dropbox/MPI/Wingbeat/data/Myotis_myotis/mm18_1807_L1.csv")
names(mm) <- c("datetime","num","sampling_rate","raw_acc")

mm$time <- as.POSIXct(sub(":(?=[^:]*$)", ".", mm$datetime, perl = TRUE), format = "%Y-%m-%d %H:%M:%OS")
mm$raw_acc <- str_squish(mm$raw_acc)
mm$raw_acc[4]

# Function to count values separated by spaces
count_values <- function(x) {
  lengths(strsplit(x, " "))
}

# Apply the function to the column
result <- unname(sapply(mm$raw_acc, count_values))
(result/3) %>% summary
which(result != 153)

amm <- accelerateR::move_trans(data = mm[4:(nrow(mm)-1),], timestamp = "time", acc = "raw_acc",
                               id = "num", sample_frequency = "sampling_rate", no_cores = 10)

summ <- accelerateR::sum_data(amm, time = "timestamp", x = "x", y = "y", z = "z", stats = "all")
fftm <- accelerateR::sum_data(amm[1:5100,], time = "timestamp", z = "z", burstcount = 51, stats = "FFT")

save(fftm, summ, amm, mm, file = "../../../Dropbox/MPI/Wingbeat/data/Myotis_myotis/mm18_1807_L1.robj")
summ


plot(summ$Pitch)
plot(summ$ODBA)
plot(amm$z)


freqs <- data.frame(time = ACC$timestamp %>% unique,
                    freq = NA, amp = NA,
                    wfreq = NA, wamp = NA,
                    rms = NA, rms_filter = NA,
                    behavior = NA)
j <- 5
for(j in 1:nrow(summ)){
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


if ("FFT" %in% stats) {
  fastf = function(axis) {
    output <- Mod(fft(axis)/length(axis))
    output <- output[1:(length(axis)%/%2)]
    return(output)
  }
  if (!is.null(x)) {
    xnames <- paste("x", c(1:(burstcount%/%2)), sep = ".")
    fft_datax <- data[, transpose(list(fastf(x))), by = time]
    old_namesx <- names(fft_datax)[-1]
    setnames(fft_datax, old_namesx, xnames)
    results <- merge(results, fft_datax, by = time,
                     all = TRUE)
  }
  if (!is.null(y)) {
    ynames <- paste("y", c(1:(burstcount%/%2)), sep = ".")
    fft_datay <- data[, transpose(list(fastf(y))), by = time]
    old_namesy <- names(fft_datay)[-1]
    setnames(fft_datay, old_namesy, ynames)
    results <- merge(results, fft_datay, by = time,
                     all = TRUE)
  }
  if (!is.null(z)) {
    znames <- paste("z", c(1:(burstcount%/%2)), sep = ".")
    fft_dataz <- data[, transpose(list(fastf(z))), by = time]
    old_namesz <- names(fft_dataz)[-1]
    setnames(fft_dataz, old_namesz, znames)
    results <- merge(results, fft_dataz, by = time,
                     all = TRUE)
  }
}