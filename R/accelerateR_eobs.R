# AccelerateR loop to extract fft and params from acceleration data
# Edward Hurme
# 7/4/23


# devtools::install_github("wanjarast/accelerateR")

library(pacman)
p_load(data.table, janitor, accelerateR, move)
# df <- fread("./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/tag_GPA-10_8147_S1-annotated-bursts-gps.csv")
# df$individual_local_identifier <- "GPA-10_8147_S1"
# acc <- move_trans(data =  df, timestamp = "timestamp", acc = "eobs_accelerations_raw",
#                   id = "individual_local_identifier",
#                   sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
#                   naxes = 3, no_cores = 4)
path <- "./../../../ownCloud/Firetail/Eidolonhelvum/Model_tag_2396/"
files <- list.files(path, pattern = "*.csv")
bats <- sapply(strsplit(files, split = "-"), '[', 1)
if(!dir.exists(paste0(path, "accelerateR"))){
  dir.create(paste0(path, "accelerateR"))
}
i=2
for(i in 6:length(files)){
  print(bats[i])
  df <- fread(paste0(path, files[i]))
  df$tag_local_identifier <- bats[i]
  ACC <- {}
  ACC <- move_trans(data =  df[df$type == "acc",], timestamp = "timestamp", acc = "eobs_accelerations_raw",
                    id = "tag_local_identifier",
                    sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
                    naxes = 3, no_cores = 4)
  burstcount = df$eobs_accelerations_raw[1] %>% strsplit(split = " ") %>% unlist %>% length/3
  count_test(ACC, burstcount)

  df$behavior <- "resting"
  df$behavior[which(df$`annotation-layer-commuting` != "")] <- "commuting"
  df$behavior[which(df$`annotation-layer-foraging` != "")] <- "foraging"

  ACC$behavior <- rep(df$behavior[df$type == "acc"],
                      each = burstcount)
  table(ACC$behavior)
  ## fast fourier transform
  if(nrow(ACC[ACC$behavior == "commuting",]) > 0){
    fft_acc <- sum_data(ACC[ACC$behavior == "commuting",], time = "timestamp",
                        burstcount = 792/3,
                        #x="x" , y="y" ,
                        z="z" ,
                        stats = "FFT")

    image(fft_acc[,3:133] %>% as.matrix)
    freqs <- data.frame(time = df$timestamp[df$behavior == "commuting"], freq = NA, amp = NA)
    j <- 5
    for(j in 1:nrow(fft_acc)){
      # fft_acc[i,3:133] %>% as.numeric() %>% plot
      idx <- which.max(fft_acc[j,3:133])
      # abline(v = idx)
      freqs$freq[j] <- names(idx) %>% substr(3,nchar(names(idx)[1])) %>% as.numeric
      freqs$amp[j] <- max(fft_acc[j, 3:133])
    }
    freqs$duration_of_burst <- ACC$burst_size[1]/ACC$sample_frequency[1]
    freqs$frequency <- freqs$freq/freqs$duration_of_burst
    png(file = paste0(path, "/accelerateR/wingbeat_", bats[i], ".png"),
        width = 800, height = 600)
    with(freqs, plot(time, frequency, cex = 2*amp/max(amp, na.rm = TRUE),
                     pch = 1, main = bats[i]))
    dev.off()
  }

  sum_acc <- sum_data(ACC, time = "timestamp", x="x" ,
                      y="y" , z="z" , stats = "all",
                      behaviour = "behavior")
  if(nrow(ACC[ACC$behavior == "commuting",]) > 0){
    save(sum_acc, freqs, file = paste0(path, "accelerateR/", bats[i], ".robj"))
  }
  if(nrow(ACC[ACC$behavior == "commuting",]) == 0){
    save(sum_acc, file = paste0(path, "accelerateR/", bats[i], "_no_freq.robj"))
  }
}
