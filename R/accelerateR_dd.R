# AccelerateR loop to extract fft and params from acceleration data
# Edward Hurme
# 11/4/23

library(pacman)
p_load(data.table, janitor, accelerateR, move)
path <- "./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/"
df <- fread("./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/tag_GPA-10_8147_S1-annotated-bursts-gps.csv")

files <- list.files(path, pattern = "*.csv")
bats <- sapply(strsplit(files, split = "_S1"), '[', 1)
if(!dir.exists(paste0(path, "accelerateR"))){
  dir.create(paste0(path, "accelerateR"))
}
i=2
for(i in 1:length(files)){
  print(bats[i])
  df <- fread(paste0(path, files[i]))
  df$tag_local_identifier <- bats[i]
  if(any(names(df) ==  "burst-start-timestamp")){
    df$timestamp <- df$`burst-start-timestamp`
  }
  ACC <- {}
  ACC <- move_trans(data =  df[df$type == "acc",], timestamp = "timestamp", acc = "eobs_accelerations_raw",
                    id = "tag_local_identifier",
                    sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
                    naxes = 3, no_cores = 20)
  burstcount = df$eobs_accelerations_raw[1] %>% strsplit(split = " ") %>% unlist %>% length/3
  count_test(ACC, burstcount)

  df$behavior <- "resting"
  df$behavior[which(df$`annotation-layer-commuting` != "")] <- "commuting"
  df$behavior[which(df$`annotation-layer-foraging` != "")] <- "foraging"

  ACC$behavior <- rep(df$behavior[df$type == "acc"],
                      each = burstcount)
  if(length(unique(ACC$behavior)) == 3){
    ACC$behavior <- factor(ACC$behavior, levels = c("commuting", "foraging", "resting"))
  }
  if(length(unique(ACC$behavior)) == 4){
    ACC$behavior <- factor(ACC$behavior, levels = c("commuting", "foraging", "resting", "unknown"))
  }

  # change burst size
  new_burst <- TRUE
  if(new_burst == TRUE){
    duration <- 14 # seconds
    burstcount <- duration*ACC$sample_frequency[1]
    ACC$burst_size <- burstcount

    # step burstcount factors at a time and estimate the median factor
    behavior <- sapply(seq(1, length(ACC$behavior)/burstcount), function(k) {
      start <- (k-1)*burstcount+1
      end <- min((k-1)*burstcount+burstcount, length(ACC$behavior))
      median(as.numeric(ACC$behavior[start:end]))
    })

    # create a new vector of the same length as the original vector
    new_behavior <- rep(behavior, each=burstcount)

    # add the remaining median factor to the new vector
    if(length(new_behavior) < length(ACC$behavior)){
      remaining_median <- median(as.numeric(as.character(ACC$behavior[(length(new_behavior)+1):length(ACC$behavior)])))
      new_behavior <- c(new_behavior, rep(remaining_median, length(ACC$behavior) - length(new_behavior)))
    } else if (length(new_behavior) > length(ACC$behavior)) {
      new_behavior <- new_behavior[1:length(ACC$behavior)]
    }

    # ACC <- ACC[1:length(new_behavior),]

    # get timestamp every burstcount
    ACC$time <- NA
    time <- rep(ACC$timestamp[seq(1, length(ACC$behavior),
                                  by=burstcount)], each = burstcount)
    ACC$time <- time[1:length(new_behavior)]

    # check the length of the new vector
    if(length(new_behavior) == nrow(ACC)){
      if(length(unique(ACC$behavior)) == 3){
        ACC$behavior[which(new_behavior == 1)] <- "commuting"
        ACC$behavior[which(new_behavior == 2)] <- "foraging"
        ACC$behavior[which(new_behavior == 3)] <- "resting"
      }
    }
  }
  ACC <- count_test(ACC, burstcount)

  fft_acc <- sum_data(ACC[1:14000,], time = "time",
           burstcount = burstcount,
           #x="x" , y="y" ,
           z="z" ,
           stats = "FFT")


  table(ACC$time[ACC$behavior == "commuting"])
  ## fast fourier transform
  if(nrow(ACC[ACC$behavior == "commuting" |
              ACC$behavior == "foraging",]) > 0){
    fft_acc <- sum_data(ACC[which(ACC$behavior == "commuting" |
                                    ACC$behavior == "foraging"),], time = "time",
                        burstcount = burstcount,
                        #x="x" , y="y" ,
                        z="z" ,
                        stats = "FFT")

    image(fft_acc[,3:((burstcount/2)+1)] %>% as.matrix)
    freqs <- data.frame(time = df$timestamp[df$behavior == "commuting"|
                                              df$behavior == "foraging"],
                        freq = NA, amp = NA,
                        behavior = df$behavior[df$behavior == "commuting"|
                                                 df$behavior == "foraging"])
    j <- 5
    for(j in 1:nrow(fft_acc)){
      # fft_acc[i,3:133] %>% as.numeric() %>% plot
      idx <- which.max(fft_acc[j,3:((burstcount/2) + 1)])
      # abline(v = idx)
      freqs$freq[j] <- names(idx) %>% substr(3,nchar(names(idx)[1])) %>% as.numeric
      freqs$amp[j] <- max(fft_acc[j, 3:((burstcount/2)+1)])
    }
    freqs$duration_of_burst <- ACC$burst_size[1]/ACC$sample_frequency[1]
    freqs$frequency <- freqs$freq/freqs$duration_of_burst
    png(file = paste0(path, "/accelerateR/wingbeat_", bats[i], ".png"),
        width = 800, height = 600)
      with(freqs, plot(time, frequency, cex = 2*amp/max(amp, na.rm = TRUE),
                     pch = 1, main = bats[i],
                     col = behavior %>% factor(levels = c("commuting", "foraging"))))
    dev.off()
  }

  sum_acc <- sum_data(ACC, time = "timestamp", x="x" ,
                      y="y" , z="z" , stats = "all",
                      behaviour = "behavior")
  if(nrow(ACC[ACC$behavior == "commuting"|ACC$behavior == "foraging",]) > 0){
    save(sum_acc, freqs, file = paste0(path, "accelerateR/", bats[i], ".robj"))
  }
  if(nrow(ACC[ACC$behavior == "commuting"|ACC$behavior == "foraging",]) == 0){
    save(sum_acc, file = paste0(path, "accelerateR/", bats[i], "_no_freq.robj"))
  }
}
