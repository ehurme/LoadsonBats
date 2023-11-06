# read accelerateR data

library(pacman)
<<<<<<< HEAD
p_load(data.table, janitor, accelerateR, signal)
paths <- c("../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Acerodonjubatus/tag_1521/",
           "../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Pteropuslylei/Model_tag_2268/",
           "../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Eidolonhelvum/Model_tag_2396/",
           "../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/",
           "../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Phyllostomushastatus/Model_tag_7CE02AF_main/",
           "../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Myotisvivesi/Mviv17_60_model/",
           "../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Myotisvivesi/Mviv18_07_model/",
           "../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Myotisvivesi/Mviv19_10_model/")

locations = c("Philippines","Thailand","Ghana", "Spain", "Panama", "Mexico", "Mexico", "Mexico")
Frequency <- data.frame()
i = 4
for(i in 1:length(locations)){
  files <- list.files(paste0(paths[i], "accelerateR/"),
                      pattern = "*.robj", full.names = TRUE)
  files <- files[!grepl(x = files, pattern = "no_freq")]
  files <- files[!grepl(x = files, pattern = "S1.robj")]
  bats <- sapply(strsplit(list.files(paste0(paths[i], "accelerateR/"),pattern = "*.robj", recursive = TRUE), split = ".robj"), "[", 1)
  Freq <- data.frame()
  j = 1
  for(j in 1:length(files)){
    freqs <- data.frame()
    load(files[j])
    print(files[j])
    freqs <- na.omit(freqs)
    plot(freqs$time, freqs$frequency, main = paste0(bats[j], " ", freqs$time[1]),
         col = factor(freqs$behavior), pch = 16)
    if(nrow(freqs) > 0){
      freqs$bat <- bats[j]
      idx <- which(sum_acc$time %in% freqs$time)
      freqs$amplitude <- NA
      if(length(idx) == nrow(freqs)){
        freqs$amplitude <- sum_acc$rangez[idx]
      }
      if(length(idx) != nrow(freqs)){
        k = 1
        for(k in 1:nrow(freqs)){
          idx <- which.min(abs(freqs$time[k] - sum_acc$time))
          freqs$amplitude <- sum_acc$rangez[idx]
        }
      }

      if(length(names(freqs)) != 8){
        freqs$bat <- NA
        freqs$amplitude <- NA
      }
      Freq <- rbind(Freq, freqs)
    }
=======
p_load(data.table, janitor, accelerateR, move)
path <- "../../../ownCloud/Firetail/Acerodonjubatus/tag_1521/"
path <- "../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/"
path <- "../../../ownCloud/Firetail/Eidolonhelvum/Model_tag_2396/"
files <- list.files(paste0(path, "accelerateR/"),
                 pattern = "*.robj", full.names = TRUE)
bats <- sapply(strsplit(list.files(paste0(path, "accelerateR/"),pattern = "*.robj"), split = ".robj"), "[", 1)
location = "Philippines"
location = "Thailand"
location = "Ghana"
location = "Spain"

Freq <- data.frame()
i = 1
for(i in 1:length(files)){
  freqs <- data.frame()
  load(files[i])
  if(nrow(freqs) > 0){
    freqs$bat <- bats[i]
    Freq <- rbind(Freq, freqs)
>>>>>>> eb72c51c7d12df6004711475bf3be2b3540a89cd
  }
}


if(location == "Philippines"){
  long <- 120.4
  lat <- 14.78
}
if(location == "Thailand"){
  long <- 101.1
  lat <- 13.52
}
if(location == "Ghana"){
  long <- -0.19
  lat <- 5.55
}
if(location == "Spain"){

}

Freq
Freq$sunrise <- maptools::sunriset(crds = matrix(c(long,
                                              lat), nrow = 1), dateTime = as_datetime(Freq$time),
                              POSIXct.out = T, direction = "sunrise")$time
Freq$sunset <- maptools::sunriset(crds = matrix(c(long, lat),
                                           nrow = 1), dateTime = as_datetime(Freq$time),
                             POSIXct.out = T, direction = "sunset")$time
Freq[which(Freq$time > Freq$sunset | Freq$time < Freq$sunrise),]

with(Freq[which(Freq$frequency > 2 &
                   Freq$frequency < 5 &
                  hour(Freq$time) > hour(Freq$sunset[1]) &
                  hour(Freq$time) < hour(Freq$sunrise[1]) &
                   Freq$behavior == "commuting"),],
     plot(hour(time), frequency, col = rgb(0,0,0,.1), pch = 16))
abline(v = hour(Freq$sunset[1]), col = "blue")
abline(v = hour(Freq$sunrise[1]), col = "orange")

ggplot(Freq[which(Freq$frequency > 2 &
                     Freq$frequency < 5 &
                     hour(Freq$time) > hour(Freq$sunset[1]) &
                     hour(Freq$time) < hour(Freq$sunrise[1]) &
                     Freq$behavior == "commuting"),],
       aes(x = hour(time), y = frequency))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm")
Freq$hour <- hour(Freq$time)
Freq$hours_from_sunrise # ask chat gpt

library(glmmTMB)

fit <- glmmTMB(frequency~amp+hour+(1|bat),
           data = Freq[which(Freq$frequency > 2 &
                               Freq$frequency < 6 &
                               # hour(Freq$time) > hour(Freq$sunset[1]) &
                               # hour(Freq$time) < hour(Freq$sunrise[1]) &
                               Freq$behavior == "commuting"),])
summary(fit)
# wingbeat mass
f1 <- predict(fit, newdata = data.frame(hour = hour(Freq$sunset[1])))
f2 <- predict(fit, newdata = data.frame(hour = hour(Freq$sunrise[1])))

m1 <- 1450#
m2 <- (f2/f1)^2 * m1
m2

# for frugivores, how does time between commutes influence weight gain?

