library(pacman)
p_load(tidyverse, data.table, lubridate)
op <- options(digits.secs=3)
matlab2POS = function(x, timez = "UTC") {
  days = x - 719529 	# 719529 = days from 1-1-0000 to 1-1-1970
  secs = days * 86400 # 86400 seconds in a day
  # This next string of functions is a complete disaster, but it works.
  # It tries to outsmart R by converting the secs value to a POSIXct value
  # in the UTC time zone, then converts that to a time/date string that
  # should lose the time zone, and then it performs a second as.POSIXct()
  # conversion on the time/date string to get a POSIXct value in the user's
  # specified timezone. Time zones are a goddamned nightmare.
  return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1',
                                        tz = 'UTC'), format = '%Y-%m-%d %H:%M:%OS',
                             tz = 'UTC', usetz = FALSE), tz = timez))
}


acc <- fread("../../../Dropbox/MPI/Wingbeat/data/Lasiurus_cinereus/DataDepository/Accelerometer Files/Full_ACC_29507A30.csv")
acc <- acc[,1:4]
acc$time <- matlab2POS(acc$ACCtime)
dt <- diff(acc$time) %>% as.numeric
summary(dt)
sampling_rate <- 1/.025
which(dt > 1) %>% diff
burst_size <- 393*.025
acc$ACCZ %>% plot(type = "l")

plot(acc$ACCtime, acc$ACCZ, type = "l")

# identify bursts
diff(acc$ACCtime) %>% plot

