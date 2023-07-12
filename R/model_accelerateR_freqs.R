# read accelerateR data

library(pacman)
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

