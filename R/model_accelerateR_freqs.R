# read accelerateR data

library(pacman)
p_load(data.table, janitor, accelerateR)
paths <- c("../../../ownCloud/Firetail/Acerodonjubatus/tag_1521/",
           "../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/",
           "../../../ownCloud/Firetail/Eidolonhelvum/Model_tag_2396/",
           "../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/",
           "../../../ownCloud/Firetail/Phyllostomushastatus/Model_tag_7CE02AF_main/",
           "../../../ownCloud/Firetail/Myotisvivesi/Mviv17_60_model/",
           "../../../ownCloud/Firetail/Myotisvivesi/Mviv18_07_model/",
           "../../../ownCloud/Firetail/Myotisvivesi/Mviv19_10_model/")

locations = c("Philippines","Thailand","Ghana", "Spain", "Panama", "Mexico", "Mexico", "Mexico")
Frequency <- data.frame()
i = 1
for(i in 4:length(locations)){
  files <- list.files(paste0(paths[i], "accelerateR/"),
                      pattern = "*.robj", full.names = TRUE)
  files <- files[!grepl(x = files, pattern = "no_freq")]
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
      Freq <- rbind(Freq, freqs)
    }
  }
  # Freq[Freq$bat == 'individual_Mviv18_27',]
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
  if(locations[i] == "Panama"){
    Freq$long <- -82.3
    Freq$lat <- 9.4
    Freq$species <- "Phastatus"
  }
  if(locations[i] == "Mexico"){
    Freq$long <- -113
    Freq$lat <- 28.9
    Freq$species <- "Mvivesi"
    if(year(Freq$time[i]) == 2017){
      #Freq$local_time <- Freq$time
      Freq$time <- Freq$time + 7*3600
    }
    b <- unique(Freq$bat)
    k = 7
    for(k in 1:length(b)){
      f_idx <- which(Freq$bat == b[k])
      if(Freq$time[f_idx[1]] %>% hour < 11){
        Freq$time[f_idx] <- Freq$time[f_idx] - 7*3600
      }
      if(Freq$time[f_idx[1]] %>% hour > 15){
        Freq$time[f_idx] <- Freq$time[f_idx] + 7*3600
      }
    }
  }

  # Get the sunset time for the given location and time
  dates <- as.Date(Freq$time) %>% unique()
  dates <-seq.Date(dates[1]-1, dates[length(dates)]+1, by = 1)
  sunset_times <- suncalc::getSunlightTimes(date = dates,
                                            lon = Freq$long[1], lat = Freq$lat[1],
                                            keep = c("sunset", "sunrise"))
  Freq$hours_since_sunset <- NA

  for(j in 1:nrow(Freq)){
    sunset_idx <- which.min(abs(Freq$time[j] - sunset_times$sunset))
    Freq$sunset[j] <- sunset_times$sunset[sunset_idx] %>% as.character()
    Freq$hours_since_sunset[j] <- difftime(Freq$time[j], sunset_times$sunset[sunset_idx], units = "hours") %>% as.numeric
  }
  hist(hour(Freq$time))
  hist(Freq$hours_since_sunset)

  # Frequency <- rbind(Frequency, Freq)
  print(length(which(is.na(Frequency$freq))))
  save(Freq, file = paste0("../../../Dropbox/MPI/Wingbeat/data/Frequency_",Freq$species[1], "_", year(Freq$time[1]), ".robj"))
}

# save(Frequency, file = "../../../Dropbox/MPI/Wingbeat/data/Frequency_allspecies.robj")


Frequency$bat %>% table
hist(Frequency$hours_since_sunset)
Frequency[which(Frequency$hours_since_sunset < -2),] %>% View()

# save(Frequency, file = "../../../Dropbox/MPI/Wingbeat/data/Frequency_Phastatus.robj")
save(Frequency, file = "../../../Dropbox/MPI/Wingbeat/data/Frequency_Mvivesi.robj")
# load("../../../Dropbox/MPI/Wingbeat/data/Frequency_Mvivesi.robj")

table(Frequency$bat)
summary(Frequency)
table(Frequency$species)
Frequency[which(Frequency$frequency %>% is.na),]

ggplot(Frequency %>% filter(frequency > 1 & amp > 25), aes(x = frequency, y = amp))+
  geom_point(alpha = 0.1)
ggplot(Frequency %>% filter(frequency > 1 & amp > 25 & hours_since_sunset > -5),
       aes(x = hours_since_sunset, y = frequency))+
  geom_point(aes(cex = amp), alpha = 0.1)+
  geom_smooth(method = "lm")+
  ylim(c(2,4))+
  facet_wrap(~date(time))

ggplot(Frequency[Frequency$frequency > 1 &
                   Frequency$frequency < 20 &
                   Frequency$amp > .5 &
                   Frequency$behavior != "commuting" &
                   Frequency$hours_since_sunset > -2,] %>% na.omit,
       aes(x = frequency, fill = bat))+
         geom_histogram(binwidth = 0.2)+facet_wrap(~species, scales = "free")+
  theme(legend.position = "none")

ggplot(Frequency[Frequency$frequency > 2 &
                 Frequency$frequency < 20 &
                 Frequency$amp > .5 &
                 Frequency$behavior != "commuting" &
                 Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = amp))+geom_histogram()+facet_wrap(~species, scales = "free")

ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 20 &
                   Frequency$amp > .5 &
                   Frequency$behavior != "commuting" &
                   Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = amplitude))+geom_histogram()+facet_wrap(~species, scales = "free")

ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 20 &
                   Frequency$amp > 0.5 &
                   Frequency$behavior != "resting" &
                   Frequency$hours_since_sunset > -2,] %>% na.omit,
       aes(x = amplitude, y = frequency, col = bat, group = species))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm")+
  facet_wrap(~species)+theme(legend.position = "none")

Frequency$power <- Frequency$frequency^3 * Frequency$amplitude^3
ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 20 &
                   Frequency$amp > .5 &
                   Frequency$behavior == "commuting" &
                   Frequency$hours_since_sunset > -12,] %>% na.omit,
       aes(x = power))+geom_histogram(bins = 50)+
  facet_wrap(~species, scales = "free")

################################################################################
# time
################################################################################
Frequency$date
## Frequency
ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 20 &
                   Frequency$amp > .5 &
                   Frequency$behavior != "resting" &
                   Frequency$hours_since_sunset > -2,] %>% na.omit,
       aes(x = hours_since_sunset, y = frequency, #group = paste0(bat, date),
           col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5)+
  facet_wrap(~species, scales = "free")+theme(legend.position = "none")

## Amplitude
ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 20 &
                   Frequency$amp > .5 &
                   Frequency$behavior == "commuting" &
                   Frequency$hours_since_sunset > -2,] %>% na.omit,
       aes(x = hours_since_sunset, y = amp, group = bat, col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5)+
  facet_wrap(~species, scales = "free")+theme(legend.position = "none")
## Power
ggplot(Frequency[Frequency$frequency > 2 &
                   Frequency$frequency < 20 &
                   Frequency$amp > .5 &
                   Frequency$behavior == "commuting" &
                   Frequency$hours_since_sunset > -2,] %>% na.omit,
       aes(x = hours_since_sunset, y = power, group = bat, col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5)+
  facet_wrap(~species, scales = "free")+theme(legend.position = "none")


################################################################################
# model freq
################################################################################
library(glmmTMB)
library(DHARMa)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# all bats
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Frequency$date <- date(Frequency$sunset)

fit_freq_all <- glmmTMB(frequency~round(hours_since_sunset,0)+(1+date|bat+species),
                     #family = Gamma(link = "log"),
                  data = Frequency[Frequency$frequency > 2 &
                                     Frequency$frequency < 20 &
                                     Frequency$amp > .5 &
                                     Frequency$behavior != "resting" &
                                     Frequency$hours_since_sunset > -2,] %>% na.omit)
summary(fit_freq_all)
simres <- simulateResiduals(fit_freq_all)
plot(simres)

#-------------------------------------------------------------------------------
# A jubatus
#-------------------------------------------------------------------------------
fit_freq_aj <- glmmTMB(frequency~round(hours_since_sunset,0)+(1+date|bat),
                        #family = Gamma(link = "log"),
                        data = Frequency[Frequency$frequency > 2 &
                                           Frequency$frequency < 20 &
                                           Frequency$amp > .5 &
                                           Frequency$behavior != "resting" &
                                           Frequency$hours_since_sunset > -2 &
                                         Frequency$species == "Ajubatus",] %>% na.omit)
summary(fit_freq_aj)
simres <- simulateResiduals(fit_freq_aj)
plot(simres)

#-------------------------------------------------------------------------------
# Eidolon helvum
#-------------------------------------------------------------------------------
fit_freq_eh <- glmmTMB(frequency~round(hours_since_sunset,0)+(1+date|bat),
                       #family = Gamma(link = "log"),
                       data = Frequency[Frequency$frequency > 2 &
                                          Frequency$frequency < 20 &
                                          Frequency$amp > .5 &
                                          Frequency$behavior != "resting" &
                                          Frequency$hours_since_sunset > -2 &
                                          Frequency$species == "Ehelvum",] %>% na.omit)
summary(fit_freq_eh)
simres <- simulateResiduals(fit_freq_eh)
plot(simres)

#-------------------------------------------------------------------------------
# A jubatus
#-------------------------------------------------------------------------------
fit_freq_aj <- glmmTMB(frequency~round(hours_since_sunset,0)+(1|bat+date),
                       #family = Gamma(link = "log"),
                       data = Frequency[Frequency$frequency > 2 &
                                          Frequency$frequency < 20 &
                                          Frequency$amp > .5 &
                                          Frequency$behavior != "resting" &
                                          Frequency$hours_since_sunset > -2 &
                                          Frequency$species == "Ajubatus",] %>% na.omit)
summary(fit_freq_aj)
simres <- simulateResiduals(fit_freq_aj)
plot(simres)

#-------------------------------------------------------------------------------
# A jubatus
#-------------------------------------------------------------------------------
fit_freq_aj <- glmmTMB(frequency~round(hours_since_sunset,0)+(1|bat+date),
                       #family = Gamma(link = "log"),
                       data = Frequency[Frequency$frequency > 2 &
                                          Frequency$frequency < 20 &
                                          Frequency$amp > .5 &
                                          Frequency$behavior != "resting" &
                                          Frequency$hours_since_sunset > -2 &
                                          Frequency$species == "Ajubatus",] %>% na.omit)
summary(fit_freq_aj)
simres <- simulateResiduals(fit_freq_aj)
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

