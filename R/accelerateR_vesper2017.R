# AccelerateR loop to extract fft and params from acceleration data
# Edward Hurme
# 11/4/23

library(pacman)
p_load(data.table, janitor, accelerateR, move)
path <- "./../../../ownCloud/Firetail/Myotisvivesi/Mviv17_60_model/"
path <- "./../../../ownCloud/Firetail/Myotisvivesi/Mviv18_07_model/"
path <- "./../../../ownCloud/Firetail/Myotisvivesi/Mviv19_10_model/"
path <- "./../../../ownCloud/Firetail/Phyllostomushastatus/Model_tag_7CE02AF_main/"
path <- "./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/"
# df <- fread("./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/tag_GPA-10_8147_S1-annotated-bursts-gps.csv")

files <- list.files(path, pattern = "*.csv")
bats <- sapply(strsplit(files, split = "-a"), '[', 1)

if(!dir.exists(paste0(path, "accelerateR"))){
  dir.create(paste0(path, "accelerateR"))
}

i=1
for(i in 2:length(files)){
  print(bats[i])
  df <- fread(paste0(path, files[i]))
  df$tag_local_identifier <- bats[i]
  idx <- grep(x = names(df), pattern = "time")[1]
  df$timestamp <- as.data.frame(df)[,idx]

  if(year(df$timestamp[1]) > 2017){
    sampling_rate <- 50
  }
  if(grepl(x = bats[i], pattern = "GPA")){
    sampling_rate <- 100
  }
  if(year(df$timestamp[1]) == 2017){
    sampling_rate <- 40
  }
  if(year(df$timestamp[1]) > 2021){
    sampling_rate <- 25
  }

  df$behavior <- "resting"
  if(any(names(df) == "annotation_layer_commuting")){
    df$behavior[which(df$annotation_layer_commuting != "")] <- "commuting"
    df$behavior[which(df$annotation_layer_foraging != "")] <- "foraging"
  }
  if(!any(names(df) == "annotation_layer_commuting")){
    df$behavior[which(df$`annotation-layer-commuting` != "")] <- "commuting"
    df$behavior[which(df$`annotation-layer-foraging` != "")] <- "foraging"
  }

  df$behavior %>% table
  # set burst size
  new_burst <- TRUE
  if(new_burst == TRUE){
    duration <- 14 # seconds
    burstcount <- duration*sampling_rate

    if(any(names(df) == "eobs_accelerations_raw")){
      ACC <- {}
      ACC <- move_trans(data =  df[df$type == "acc",], timestamp = "timestamp", acc = "eobs_accelerations_raw",
                        id = "tag_local_identifier",
                        sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
                        naxes = 3, no_cores = 20)
      ACC$behavior <- rep(df$behavior[df$type == "acc"],
                          each = burstcount)[1:nrow(ACC)]
      table(ACC$behavior)/burstcount
    }
    if(!any(names(df) == "eobs_accelerations_raw")){
      ACC <- data.table(x = df$`acceleration-x`, y = df$`acceleration-y`, z = df$`acceleration-z`,
                        timestamp = df$timestamp, sample_frequency = sampling_rate,
                        behavior = factor(df$behavior, levels = c("commuting", "foraging", "resting")),
                        burst_size = burstcount)
    }
    count_test(ACC, burstcount = ACC$burst_size[1])

    # floor(nrow(ACC)/burstcount)
    ACC <- ACC[1:(burstcount*floor(nrow(ACC)/burstcount)),]
    table(ACC$behavior)/burstcount
    ACC$burst_size <- burstcount
    # step burstcount factors at a time and estimate the median factor
    behavior <- sapply(seq(1, length(ACC$behavior)/burstcount), function(k) {
      start <- (k-1)*burstcount+1
      end <- min((k-1)*burstcount+burstcount, length(ACC$behavior))
      median(as.numeric(ACC$behavior[start:end]), na.rm = TRUE)
    })
    behavior <- round(behavior, 0)

    # create a new vector of the same length as the original vector
    new_behavior <- rep(behavior, each=burstcount)

    # get timestamp every burstcount
    ACC$time <- NA
    time <- rep(ACC$timestamp[seq(1, length(ACC$behavior),
                                  by=burstcount)], each = burstcount)

    nrow(ACC)/burstcount
    length(new_behavior)/burstcount
    length(time)/burstcount

    ACC$time <- time[1:nrow(ACC)]

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
  ACC[which(ACC$keep == FALSE),]
  # ACC <- ACC[-which(ACC$keep == FALSE),]
  nrow(ACC)/burstcount
  table(ACC$behavior)/burstcount
  summary(ACC)

  ## fast fourier transform
  # foraging_idx <- which(ACC$behavior == "commuting" |
  #                         ACC$behavior == "foraging")
  foraging_idx <- 1:nrow(ACC)
  if(nrow(ACC[foraging_idx,]) > 0){
    fft_acc <- {}
    try({
      fft_acc <- sum_data(data = ACC[foraging_idx,] %>% na.omit(),
                        time = "time",
                        burstcount = burstcount,
                        #x="x" , y="y" ,
                        z="z",
                        stats = "FFT")
      })
    png(file = paste0(path, "/accelerateR/fft_", bats[i], ".png"),
        width = 800, height = 600)
      image(fft_acc[,3:((burstcount/2)+1)] %>% as.matrix)
    dev.off()

    freqs <- data.frame(time = ACC$time[foraging_idx] %>% unique(),
                        freq = NA, amp = NA,
                        behavior = ACC$behavior[foraging_idx[seq(1, length(foraging_idx), by=burstcount)]])
    j <- 5
    for(j in 1:nrow(fft_acc)){
      # fft_acc[i,3:133] %>% as.numeric() %>% plot
      idx <- which.max(fft_acc[j,3:((burstcount/2) + 1)])
      # abline(v = idx)
      freqs$freq[j] <- names(idx) %>% substr(3,nchar(names(idx)[1])) %>% as.numeric
      freqs$amp[j] <- max(fft_acc[j, 3:((burstcount/2)+1)])
    }
    freqs$duration_of_burst <- ACC$burst_size[1]/median(ACC$sample_frequency)
    # freqs$fft_resolution <- median(ACC$sample_frequency)/ACC$burst_size[1]
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
  #if(nrow(ACC[ACC$behavior == "commuting"|ACC$behavior == "foraging",]) > 0){
    save(sum_acc, freqs, file = paste0(path, "accelerateR/", bats[i], ".robj"))
  #}
  # if(nrow(ACC[ACC$behavior == "commuting"|ACC$behavior == "foraging",]) == 0){
  #   save(sum_acc, file = paste0(path, "accelerateR/", bats[i], "_no_freq.robj"))
  # }
}
