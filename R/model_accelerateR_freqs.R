# read accelerateR data

library(pacman)
p_load(data.table, janitor, accelerateR, move)
paths <- c("../../../ownCloud/Firetail/Acerodonjubatus/tag_1521/",
           "../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/",
           "../../../ownCloud/Firetail/Eidolonhelvum/Model_tag_2396/",
           "../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/")

locations = c("Philippines","Thailand","Ghana", "Spain")
Frequency <- data.frame()
i = 4
for(i in 1:length(locations)){
  files <- list.files(paste0(paths[i], "accelerateR/"),
                      pattern = "*.robj", full.names = TRUE)
  bats <- sapply(strsplit(list.files(paste0(paths[i], "accelerateR/"),pattern = "*.robj"), split = ".robj"), "[", 1)
  Freq <- data.frame()
  for(j in 1:length(files)){
    freqs <- data.frame()
    load(files[j])
    if(nrow(freqs) > 0){
      freqs$bat <- bats[j]
      Freq <- rbind(Freq, freqs)
    }
  }

  if(locations[i] == "Philippines"){
    Freq$long <- 120.273
    Freq$lat <- 14.788
    Freq$species <- "Ajubatus"
  }
  if(locations[i] == "Thailand"){
    Freq$long <- 101.1
    Freq$lat <- 13.52
    Freq$species <- "Plylei"
  }
  if(locations[i] == "Ghana"){
    Freq$long <- -0.19
    Freq$lat <- 5.55
    Freq$species <- "Ehelvum"
  }
  if(locations[i] == "Spain"){
    Freq$long <- 6.44
    Freq$lat <- 36.99
    Freq$species <- "Nlasiopterus"
  }

  # Get the sunset time for the given location and time
  sunset_times <- suncalc::getSunlightTimes(date = as.Date(Freq$time) %>% unique(), lon = long, lat = lat, keep = c("sunset", "sunrise"))
  Freq$hours_since_sunset <- NA

  for(j in 1:nrow(Freq)){
    sunset_idx <- which.min(abs(Freq$time[j] - sunset_times$sunset))
    Freq$sunset[j] <- sunset_times$sunset[sunset_idx] %>% as.character()
    Freq$hours_since_sunset[j] <- difftime(Freq$time[j], sunset_times$sunset[sunset_idx], units = "hours") %>% as.numeric
  }
  Frequency <- rbind(Frequency, Freq)
}
table(Frequency$species)

ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 6 &
                   Frequency$amp > 10 &
                   Frequency$behavior == "commuting" &
                   Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = frequency))+geom_histogram()+facet_wrap(~species)

ggplot(Frequency[Frequency$frequency > 2 &
                 Frequency$frequency < 6 &
                 Frequency$amp > 10 &
                 Frequency$behavior != "commuting" &
                 Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = amp))+geom_histogram()+facet_wrap(~species)

Frequency$power <- Frequency$frequency^3 * Frequency$amp^3
ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 6 &
                   Frequency$amp > 10 &
                   Frequency$behavior != "commuting" &
                   Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = power))+geom_histogram(bins = 50)+facet_wrap(~species, scales = "free")

ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 6 &
                   Frequency$amp > 10 &
                   Frequency$behavior != "resting" &
                   Frequency$hours_since_sunset > -2,] %>% na.omit,
       aes(x = amp, y = frequency))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm")+
  facet_wrap(~species)

ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 6 &
                   Frequency$amp > 10 &
                   Frequency$behavior == "commuting" &
                   Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = hours_since_sunset, y = frequency, group = bat, col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~species, scales = "free")+theme(legend.position = "none")

ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 6 &
                   Frequency$amp > 10 &
                   Frequency$behavior == "commuting" &
                   Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = hours_since_sunset, y = amp, group = bat, col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~species, scales = "free")+theme(legend.position = "none")

ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 6 &
                   Frequency$amp > 10 &
                   Frequency$behavior == "commuting" &
                   Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = hours_since_sunset, y = frequency^3*amp^3, group = bat, col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~species, scales = "free")+
  theme(legend.position = "none")



library(glmmTMB)
library(DHARMa)

Frequency$date <- date(Frequency$sunset)

fit_power <- glmmTMB(power~round(hours_since_sunset,0)+(1|bat+date),
                     family = Gamma(link = "log"),
                  data = Frequency[Frequency$frequency > 2 &
                                     Frequency$frequency < 6 &
                                     Frequency$amp > 10 &
                                     Frequency$behavior == "commuting" &
                                     Frequency$hours_since_sunset > -2 &
                                   Frequency$species == "Ehelvum",] %>% na.omit)
summary(fit_power)
simres <- simulateResiduals(fit_power)
plot(simres)

fit_log_power <- glmmTMB(log(power)~round(hours_since_sunset,0)+(date|bat),
                     family = Gamma(link = "log"),
                     data = Frequency[Frequency$frequency > 2 &
                                        Frequency$frequency < 6 &
                                        Frequency$amp > 10 &
                                        Frequency$behavior == "commuting" &
                                        Frequency$hours_since_sunset > -2 &
                                        Frequency$species == "Ajubatus",] %>% na.omit)
summary(fit_log_power)
simres <- simulateResiduals(fit_log_power)
plot(simres)


fit_eh <- glmmTMB(frequency~amp+round(hours_since_sunset,0)+(1|bat+date),
           data = Frequency[Frequency$frequency > 2 &
                               Frequency$frequency < 6 &
                               Frequency$amp > 10 &
                               Frequency$behavior == "commuting" &
                               Frequency$hours_since_sunset > -2 &
                               Frequency$species == "Ehelvum",] %>% na.omit)
summary(fit_eh)




# wingbeat mass
f1 <- predict(fit, newdata = data.frame(hour = hour(Freq$sunset[1])))
f2 <- predict(fit, newdata = data.frame(hour = hour(Freq$sunrise[1])))

m1 <- 1450#
m2 <- (f2/f1)^2 * m1
m2

# for frugivores, how does time between commutes influence weight gain?

