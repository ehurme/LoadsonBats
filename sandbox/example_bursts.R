library(pacman)
p_load(data.table, tidyverse, lubridate)
op <- options(digits.secs=3)

#
species <- "Nyctalus_lasiopterus"
bat <- "tag_GPA-18_66BA-S1"
df <- fread("../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/tag_GPA-18_66BA_S1-annotated-bursts.csv")
idx <- which(df$`annotation-layer-commuting` != "")

data.frame(idx, c(NA, diff(idx)))
sampling_rate <- 100
bursts <- c(36306, 36002, 36128, 36140, 36155)
for(j in 1:length(bursts)){
  # 15 second bursts
  burst15 <- c(bursts[j]:(bursts[j]+8))
  i = 1
  db <- data.frame()
  for(i in 1:length(burst15)){
    b <- {}
    temp <- {}
    x <- NA
    y <- NA
    z <- NA
    time <- NA
    burst <- NA
    temp <- df$eobs_accelerations_raw[burst15[i]] |> strsplit(" ") |>
      unlist() |>
      as.numeric()
    if(length(temp) > 0){
      x <- temp[seq(1, length(temp), 3)]
      y <- temp[seq(2, length(temp), 3)]
      z <- temp[seq(3, length(temp), 3)]
      temp_time <- seq.POSIXt(from = df$`burst-start-timestamp`[burst15[i]],
                              to = df$`burst-start-timestamp`[burst15[i]]+length(temp)/sampling_rate/3,
                              length.out = length(temp)/3)
      time <- format(temp_time, "%Y-%m-%d %H:%M:%OS3") |> as.character()
      #burst <- rep(ii, length(temp_time))
    }
    b <- data.frame(x,y,z,time)
    db <- rbind(db, b)
  }

  fwrite(db, file = paste0("../../../Dropbox/MPI/Wingbeat/data/sample_bursts/", species, "_", bat, "_", burst15[1], "to", burst15[length(burst15)], ".csv"),
         row.names = FALSE)
}


plot(db$z, type= "l")
