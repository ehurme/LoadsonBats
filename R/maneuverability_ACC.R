# corcoran et al.
# For each accelerometer recording, we quantified the magnitude of flight
# maneuvering (Fig. 3) using the lateral and anterior–posterior components
# of the accelerometer measurements. We first detrended each signal by
# subtracting the mean signal value. We next applied a 4 Hz low-pass filter
# to each signal to remove wingbeat frequency oscillations and highlight
# maneuvers lasting more than one wingbeat. Finally, we took the root-mean-square
# value of the filtered signals as the overall measure of flight maneuvering.


library(pacman)

p_load(zoo, data.table, lubridate, av, seewave, tuneR, phonTools, soundgen)

df <- fread("./../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Acerodonjubatus/tag_1521/tag_1522-annotated-bursts-gps.csv")
load("./../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Acerodonjubatus/tag_1522/accelerateR/tag_1522.robj")
sampling_rate = 18.74

# detrend
idx <- which(ACC$behavior == "commuting")
bursts <- ACC$burst[idx] %>% unique

b1 <- ACC[ACC$burst == bursts[10],]
b1$x0 <- b1$x - mean(b1$x)
b1$y0 <- b1$y - mean(b1$y)
b1$z0 <- b1$z - mean(b1$z)

plot(1:nrow(b1), b1$z0, col = 1, type = "l")
lines(1:nrow(b1), b1$y0, col = 2)
lines(1:nrow(b1), b1$x0, col = 3)


w <- tuneR::Wave(left = b1$z0, samp.rate = sampling_rate, bit = 16)
duration(w)
plot(w)
meanspec(w, wl = length(w))
spectro(w, wl = length(w)/5)
dev.off()
wf <- ffilter(w, f= sampling_rate, from = 0, to = 1, bandpass = TRUE, wl = length(w)/5)
plot(w)
lines((1:nrow(b1))/sampling_rate, wf*100, lwd = 3)

rms(b1$z0)
rms(wf*100)


load("./../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Acerodonjubatus/tag_1521/accelerateR/tag_1521.robj")
files <- list.files(path = "./../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Acerodonjubatus/tag_1521/accelerateR/",
                    pattern = ".robj", full.names = TRUE)
bats <- list.files(path = "./../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Acerodonjubatus/tag_1521/accelerateR/",
                   pattern = ".robj", full.names = FALSE) %>% substr(., 1, 8)
files <- list.files(path = "./../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Pteropuslylei/Model_tag_2268/accelerateR/",
                    pattern = ".robj", full.names = TRUE)
bats <- list.files(path = "./../../../ownCloud - ehurme@ab.mpg.de@owncloud.gwdg.de/Firetail/Pteropuslylei/Model_tag_2268/accelerateR/",
                   pattern = ".robj", full.names = FALSE) %>% substr(., 1, 8)

Freqs <- data.table()
for(i in 1:length(files)){
  load(files[i])
  freqs$bat <- bats[i]
  Freqs <- rbind(Freqs, freqs, fill = TRUE)
}

with(Freqs, plot(frequency, amp, col = behavior %>% factor, cex = 0.5))
with(Freqs, plot(frequency, rms, col = behavior %>% factor, cex = 0.5))
with(Freqs, plot(frequency, rms_filter, col = behavior %>% factor, cex = 0.5))
clean_freqs <- Freqs[Freqs$frequency > 2 & Freqs$frequency < 4 &
                       Freqs$amp > 25 & Freqs$rms > 100 & Freqs$rms_filter < 200,]
# with(clean_freqs, plot(frequency, col = behavior %>% factor))
with(clean_freqs, plot(time, frequency, col = behavior %>% factor))
super_clean_freqs <- clean_freqs[clean_freqs$rms_filter < 200 #& hour(clean_freqs$time) > 8
                                 ,]

super_clean_freqs$bat_date <-  paste(super_clean_freqs$bat, date(super_clean_freqs$time))
bat_date <- table(super_clean_freqs$bat_date)
keep <- bat_date[which(bat_date %>% as.numeric > 20)] %>% names()
super_clean_freqs <- super_clean_freqs[which(super_clean_freqs$bat_date %in% keep),]
table(super_clean_freqs$bat_date)

with(super_clean_freqs, plot(time, frequency, col = behavior %>% factor))
ggplot(super_clean_freqs, aes(hour(time), frequency, col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(aes(group = bat_date),
              alpha = 0.5, se = FALSE)+
  geom_smooth(se = FALSE, lwd = 2, col = 1)+
  theme_classic()

ggsave(filename = "../../../Dropbox/MPI/Wingbeat/plots/Ajubatus_super_clean_freq.png")

ggplot(super_clean_freqs, aes(hour(time), frequency, col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(se = FALSE, lwd = 2, col = 1)+
  facet_wrap(~bat_date)



ggplot(super_clean_freqs, aes(hour(time), frequency, col = bat))+
  geom_point(alpha = 0.1)+
  geom_smooth(aes(group = bat_date),
              alpha = 0.5, se = FALSE,
              method = "lm")+
  geom_smooth(se = FALSE, lwd = 2, col = 1, method = "lm")+
  theme_classic()
ggsave(filename = "../../../Dropbox/MPI/Wingbeat/plots/Ajubatus_super_clean_freq_linear.png")

library(mgcv)
fit <- gamm(frequency ~ s(hour(time)), data = super_clean_freqs,
            random = list(bat_date = ~1))

summary(fit$lme)
summary(fit$gam)
plot(fit$gam)



library(signal)
number_of_cycles = 2
max_y = 40

x = 1:500
a = number_of_cycles * 2*pi/length(x)

y = max_y * sin(x*a)
noise1 = max_y * 1/10 * sin(x*a*10)

plot(x, y, type="l", col="red", ylim=range(-1.5*max_y,1.5*max_y,5))
points(x, y + noise1, col="green", pch=20)
points(x, noise1, col="yellow", pch=20)


bf <- butter(2, 1/50, type="low")
b1 <- filtfilt(bf, y+noise1)
points(x, b1, col="red", pch=20)

table(sum_acc$behavior)
with(sum_acc, plot(Kurtosisz, factor(behavior)))
with(sum_acc, plot(CVz, ICVz))
ggplot(sum_acc, aes(y = ODBA, x = behavior))+geom_violin()
ggplot(sum_acc, aes(y = q, x = behavior))+geom_violin()
ggplot(sum_acc[which(sum_acc$behavior != "resting"),], aes(y = ODBA, x = behavior))+geom_violin()
ggplot(sum_acc[which(sum_acc$behavior != "resting"),], aes(y = CVz, x = behavior))+
  geom_violin()+
  geom_boxplot(width = 0.25)
  #ylim(c(0,30))

library(randomForest)

rf.cf <- randomForest(factor(behavior) ~ .,
                   data = sum_acc[which(sum_acc$behavior != "resting"),],
                   importance = TRUE,
                   proximity = TRUE)
print(rf.cf)
idx <- importance(rf.cf, type = 1) %>% order
importance(rf.cf, type = 1)[idx,]
