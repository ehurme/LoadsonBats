# devtools::install_github("wanjarast/accelerateR")

library(pacman)
p_load(data.table, janitor, accelerateR, move)
df <- fread("./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/tag_GPA-10_8147_S1-annotated-bursts-gps.csv")
df$individual_local_identifier <- "GPA-10_8147_S1"
acc <- move_trans(data =  df, timestamp = "timestamp", acc = "eobs_accelerations_raw",
                  id = "individual_local_identifier",
                  sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
                  naxes = 3, no_cores = 1)

df <- fread("./../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/tag_1940-annotated-bursts-gps.csv")
df$tag_local_identifier <- "1940"
ACC <- {}
ACC <- move_trans(data =  df[df$type == "acc",], timestamp = "timestamp", acc = "eobs_accelerations_raw",
                  id = "tag_local_identifier",
                  sample_frequency = "eobs_acceleration_sampling_frequency_per_axis",
                  naxes = 3, no_cores = 10)
burstcount = df$eobs_accelerations_raw[1] %>% strsplit(split = " ") %>% unlist %>% length/3
count_test(ACC, burstcount)

df$behavior <- "resting"
df$behavior[which(df$`annotation-layer-commuting` != "")] <- "commuting"
df$behavior[which(df$`annotation-layer-foraging` != "")] <- "foraging"

ACC$behavior <- rep(df$behavior[df$type == "acc"],
                    each = burstcount)
# acto(data = ACC, time = "timestamp", behaviour = "behavior" , target_bev = c("commuting"),
#      sun = T, long = 100, lat = 14, daily = TRUE)

###############################################################################################
## fast fourier transform
fft_acc <- sum_data(ACC[ACC$behavior == "commuting",], time = "timestamp",
                    burstcount = 792/3,
                    #x="x" , y="y" ,
                    z="z" ,
                    stats = "FFT")

image(fft_acc[,3:133] %>% as.matrix)
freqs <- data.frame(time = df$timestamp[df$behavior == "commuting"], freq = NA, amp = NA)
i <- 5
for(i in 1:nrow(fft_acc)){
  fft_acc[i,3:133] %>% as.numeric() %>% plot
  idx <- which.max(fft_acc[i,3:133])
  abline(v = idx)
  freqs$freq[i] <- names(idx) %>% substr(3,nchar(names(idx)[1])) %>% as.numeric
  freqs$amp[i] <- max(fft_acc[i, 3:133])
}
freqs$duration_of_burst <- ACC$burst_size[1]/ACC$sample_frequency[1]
freqs$frequency <- freqs$freq/freqs$duration_of_burst

with(freqs, plot(time, frequency, cex = 2*amp/max(amp)))

sum_acc <- sum_data(ACC, time = "timestamp", x="x" ,
                    y="y" , z="z" , stats = "all",
                    behaviour = "behavior")

# psych::pairs.panels(sum_acc[,c("behavior", "ODBA", "meanz", "sdz", "maxz", "minz", "rangez", "covyz", "coryz")])
# psych::pairs.panels(sum_acc[,c("behavior", "ODBA", "sddyz", "mdocpz", "sdocpz", "varz", "q")])
# psych::pairs.panels(sum_acc[,c("behavior", "ODBA","Pitch", "ICVz", "CVz", "Kurtosisz", "Skewnessz")])

table(sum_acc$behavior)



# join more bats together to increase sample size
acc.rf <- randomForest::randomForest(behavior~., data = sum_acc[sum_acc$behavior != "foraging",],
                                     importance = TRUE, proximity = TRUE)

# 75% IQR index
iqr75_idx <- which(quantile(sum_acc$ODBA)[4] < sum_acc$ODBA |
             quantile(sum_acc$meanz)[4] < sum_acc$meanz |
             quantile(sum_acc$varz)[4] < sum_acc$varz )
median_sd_idx <- which((quantile(sum_acc$ODBA)[3]+sd(sum_acc$ODBA)) < sum_acc$ODBA |
                         (quantile(sum_acc$meanz)[3]+sd(sum_acc$meanz)) < sum_acc$meanz |
                         (quantile(sum_acc$varz)[3]+sd(sum_acc$varz)) < sum_acc$varz )
sum_acc$ODBA %>% plot()
sum_acc$ODBA[median_sd_idx] %>% points(median_sd_idx, ., col = 2)

sum_acc$ODBA %>% plot(col = (quantile(sum_acc$ODBA)[4] < sum_acc$ODBA)+1)
sum_acc$meanz %>% plot
sum_acc$varz %>% plot(col = (quantile(sum_acc$varz)[4] < sum_acc$varz)+1)


sum_acc$Pitch %>% plot
sum_acc$Roll %>% plot
sum_acc$Yaw %>% plot

names(sum_acc)


sum_acc$varz %>% hist(breaks = 1000, xlim = c(0,10000))


acceleration <- data.frame(time = rep(seq(5),each=20) , x = runif(n = 100,min = 1900,max=2100) ,
                           y = runif(n = 100,min = 2100,max=2300) , z = runif(n = 100,min = 1800,max=2000))

sumstats <- sum_data(data=acceleration , time="time" , x="x" ,
                     y="y" , z="z" , stats="all")
