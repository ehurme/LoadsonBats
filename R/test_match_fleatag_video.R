# flea tag DLT calibration
library(pacman)
p_load(data.table, tuneR, seewave, tidyverse)

# video sync
df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230610/CalibrationTag.csv", header = T, fill = T)
dt <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230610/testfleatagxypts.csv", header = T, fill = T)

df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230610/CalibrationTag2.csv", header = T, fill = T)
dt <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230610/testfleatag2xypts.csv", header = T, fill = T)

df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230610/CalibrationTag3.csv", header = T, fill = T)
dt <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230610/testfleatag3xypts.csv", header = T, fill = T)

# start time
df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230611/Calibration_starttime.csv", header = T, fill = T)
dt <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230611/testfleatag3xypts.csv", header = T, fill = T)
videostart_time <- hms("16:35:08.137")
videoend_time <- hms("16:37:12.353")

videostart <- 21
videoend <- 2*60+42
video_duration <- videoend - videostart
video_duration

acc_duration <- df$time[nrow(df)]
delay <- video_duration - acc_duration

# duration
df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230611/CalibrationDuration/CalibrationDuration_tag3016.csv", header = T, fill = T)
dt <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230611/testfleatag3xypts.csv", header = T, fill = T)

videostart <- 26
videoend <-14*60+15
video_duration <- videoend - videostart
video_duration # 829, 13.82
829/60

acc_duration <- df$time[nrow(df)]
delay <- video_duration - acc_duration

# remove blank rows
df <- df[!is.na(df$burstCount),]



frame_rate <- 29.924
sampling_rate <- 54
layout(1)
plot(df$accZ_mg, type = "l")
points(1:1400, df$accZ_mg[1:1400])
#df <- df[1:1000,]

df$time <- 1:nrow(df)/sampling_rate
dt$time <- 1:nrow(dt)/frame_rate

layout(cbind(1:3))
plot(df$time, df$accX_mg, type = "l")
plot(df$time, df$accY_mg, type = "l")
plot(df$time, df$accZ_mg, type = "l")
plot(dt$time, dt$pt1_cam1_X, type = "l")
plot(dt$time, dt$pt1_cam1_Y, type = "l")

acc_filter <- function(z, min, max){
  w <- tuneR::Wave(left = z, samp.rate = sampling_rate, bit = 16)
  # duration(w)
  # plot(w)
  # meanspec(w, wl = round(length(w),-1))
  # spectro(w, wl = length(w)/5)
  #dev.off()
  wf <- seewave::ffilter(w, f= sampling_rate, from = min, to = max, bandpass = TRUE, wl = length(w)/5)
  #plot(w)
  z0 <- scale(wf[,1])[,1]
  # plot(df_sub$time, df_sub$z0, type = "l")
  return(z0)
}

wy <- tuneR::Wave(left = na.omit(scale(dt$pt1_cam1_Y)), samp.rate = frame_rate, bit = 16)
layout(1)
plot(wy)
spec <- meanspec(wy)
peaks <- fpeaks(spec, nmax = 2)

min = min(peaks[,1])*1000
max = max(peaks[,1])*1000
min = 0.1
max = 4

df$x0 <- acc_filter(z = df$accX_mg, min, max)
df$y0 <- acc_filter(z = df$accY_mg, min, max)
df$z0 <- acc_filter(z = df$accZ_mg, min, max)
df$sum <- acc_filter(z = df$z+df$y+df$x, min, max)
#df$sum0 <- df_sub$x0 + df$y0 + df$z0

plot(df$time, df$x0, type = "l")
lines(df$time, scale(df$accX_mg), col = 4)

layout(cbind(1:4))
plot(df$time, df$x0, type = "l")
plot(df$time, df$y0, type = "l")
plot(df$time, df$z0, type = "l")
plot(df$time, df$sum, type = "l")

resample <- function(x,t) {
  #
  # Resample time series `x`, assumed to have unit time intervals, at time `t`.
  # Uses quadratic interpolation.
  #
  n <- length(x)
  if (n < 3) stop("First argument to resample is too short; need 3 elements.")
  i <- median(c(2, floor(t+1/2), n-1)) # Clamp `i` to the range 2..n-1
  u <- t-i
  x[i-1]*u*(u-1)/2 - x[i]*(u+1)*(u-1) + x[i+1]*u*(u+1)/2
}


e.frequency <- 54   # Herz
v.frequency <- 29.924 # Herz
x <- sapply(1:(nrow(df)*v.frequency/e.frequency),
            function(t) resample(df$x0, e.frequency*t/v.frequency))
y <- sapply(1:(nrow(df)*v.frequency/e.frequency),
            function(t) resample(df$y0, e.frequency*t/v.frequency))
z <- sapply(1:(nrow(df)*v.frequency/e.frequency),
            function(t) resample(df$z0, e.frequency*t/v.frequency))
sum_xyz <- sapply(1:(nrow(df)*v.frequency/e.frequency),
            function(t) resample(df$sum, e.frequency*t/v.frequency))
# df_sub <- data.table(x = x - mean(x, na.rm = TRUE),
#                      y = y-mean(y, na.rm = TRUE),
#                      z = z-mean(z, na.rm = TRUE),
#                      time = 1:length(e)/v.frequency)
df_sub <- data.frame(x = scale(x),
                     y = scale(y),
                     z = scale(z),
                     sum = scale(sum_xyz),
                     time = 1:length(x)/v.frequency)
layout(cbind(1:4))
with(df_sub, plot(time, x, type = "l"))
with(df_sub, plot(time, y, type = "l"))
with(df_sub, plot(time, z, type = "l"))
with(df_sub, plot(time, sum, type = "l"))

summary(df_sub)
# low pass filter

layout(cbind(1:4))
# plot(df_sub$time, df_sub$x0, type = "l")
# plot(df_sub$time, df_sub$y0, type = "l")
# plot(df_sub$time, df_sub$z0, type = "l")
plot(dt$time, dt$pt1_cam1_X, type = "l")
plot(dt$time, dt$pt1_cam1_Y, type = "l")

dt$speed <- NA
i = 2538

for(i in 1:nrow(dt)){
  #if(!is.na(dt$pt1_cam1_X[i])){
    try(dt$speed[i] <- sqrt(#(dt$pt1_cam1_X[i-1]-dt$pt1_cam1_X[i])^2 +
          (dt$pt1_cam1_Y[i-1]-dt$pt1_cam1_Y[i])^2))
  #}
}
plot(dt$speed, type = "l")

dt$acc <- c(NA, diff(dt$speed))
dt$acc0 <- scale(dt$acc)
plot(dt$acc0, type = "l")

# complete dt
c <- rle(is.na(dt$acc))

idx <- which(c$lengths == max(c$lengths[c$values == FALSE]))
c$values[idx]
start <- sum(c$lengths[1:(idx-1)])
end <- sum(c$lengths[1:(idx)])

dtc <- dt[start:end,] |> na.omit()
summary(dtc)

layout(cbind(1:2))
plot(df_sub$x, type = "l")
plot(dtc$acc, type = "l")
#plot(dtc$pt1_cam1_Y, type = "l")

# points(1:350,dtc$pt1_cam1_Y[1:350])

layout(1)

lag = 400
cor <- ccf(dtc$pt1_cam1_Y, df_sub$sum, lag.max = lag)
shift <- cor$acf |> abs() |> which.max() - lag

layout(1)
plot(df_sub$x/2+2, type = "l", ylim = c(-2, 10))
lines(1:nrow(df_sub), df_sub$y/2+4, col = 2)
lines(1:nrow(df_sub), df_sub$z/2+6, col = 3)
lines(1:nrow(df_sub), df_sub$sum/2+8, col = 4)
#lines(1:nrow(dtc)-shift,scale(-dtc$pt1_cam1_X), col = 4)
lines(1:nrow(dtc)-shift,scale(-dtc$pt1_cam1_Y), col = 6)
# lines(1:nrow(dtc)-shift,dtc$acc0/2, col = 2)
# lines(1:nrow(dtc)-shift,dtc$speed, col = 3)

layout(1)
plot((1:nrow(dtc)-shift)/sampling_rate,scale(-dtc$pt1_cam1_Y), col = 6,
     ylim = c(-2,8), type = "l")
lines(df_sub$time, df_sub$x+4)

