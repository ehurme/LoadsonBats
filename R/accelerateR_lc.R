library(pacman)
p_load(data.table, tidyverse, magrittr, dplyr, accelerateR, stringr)
op <- options(digits.secs=3)

lc <- fread("../../../Dropbox/MPI/Wingbeat/data/Lasiurus_cinereus/DataDepository/Accelerometer Files/Full_ACC_29507A30.csv")
lc <- lc[,1:4]
names(lc) <- c("x","y","z","time")

# Convert MATLAB time to POSIXct object in R
origin <- as.POSIXct("0000-01-01", format = "%Y-%m-%d", tz = "UTC") # Define MATLAB origin
lc$timestamp <- origin + as.difftime(lc$time, units = "days") # Convert MATLAB time to POSIXct
lc$timestamp %>% diff -> lc_diff

idx <- c(1,which(lc_diff > 10)+1)
i = 1
lcb <- data.frame()
for(i in 1:(length(idx)-1)){
  b <- lc[idx[i]:(idx[i+1]-1),]
  temp_acc <- {}
  for(j in 1:nrow(b)){
    temp_acc <- paste(temp_acc, round(b$x[j]),
                      round(b$y[j]), round(b$z[j]))
  }
  df <- data.frame(time = b$timestamp[1], raw_acc = temp_acc,
                     burst_size = nrow(b),
                   sampling_rate = 1/diff(b$timestamp) %>% median %>% as.numeric)
  lcb <- rbind(lcb, df)
}

lcb$id <- "29507A30"
lcb$sampling_rate <- 40
alc <- accelerateR::move_trans(data = lcb, timestamp = "time", acc = "raw_acc",
                               id = "id", sample_frequency = "sampling_rate", no_cores = 10)

summ <- accelerateR::sum_data(lc[1:34880,], time = "timestamp", x = "x", y = "y", z = "z",
                                  burstcount = 393, stats = "ODBA")

summ$meany %>% plot(type = "l")

