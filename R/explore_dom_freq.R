library(pacman)
p_load(data.table, janitor, accelerateR, signal)

load("../../../Dropbox/MPI/Wingbeat/data/Frequency_Ajubatus_2010.robj")
load("../../../Dropbox/MPI/Wingbeat/data/Frequency_Plylei_2012.robj")

Freq %>% names
# freq - fft frequency
# amp - fft amplitude
# wfreq - wave frequency


Freq$bat %>% table
plot(Freq$frequency, Freq$amp, col = Freq$behavior %>% factor)
plot(Freq$wamp, Freq$amp, col = Freq$behavior %>% factor)
plot(Freq$amplitude, Freq$amp, col = Freq$behavior %>% factor)

plot(Freq$rms, Freq$amplitude, col = Freq$behavior %>% factor)
plot(Freq$rms_filter, Freq$frequency, col = Freq$behavior %>% factor)

clean_freq <- Freq %>% dplyr::filter(frequency > 2 & rms > 150 & rms_filter < 200 & wamp > 0.9 &
                                       hours_since_sunset > -1 & behavior != "resting")

plot(clean_freq$hours_since_sunset, clean_freq$speed)
ggplot(clean_freq,#[clean_freq$behavior == "commuting",],
       aes(x = hours_since_sunset, y = speed, col = bat))+
  geom_point()+geom_smooth(se = FALSE)

ggplot(clean_freq, aes(x = behavior, y = frequency))+geom_violin()
ggplot(clean_freq, aes(x = behavior, y = amp))+geom_violin()
ggplot(clean_freq, aes(x = behavior, y = amplitude))+geom_violin()
ggplot(clean_freq, aes(x = behavior, y = speed))+geom_boxplot()
ggplot(clean_freq, aes(x = behavior, y = rms))+geom_violin()
ggplot(clean_freq, aes(x = behavior, y = rms_filter))+geom_violin()+geom_boxplot(width = 0.2)

ggplot(clean_freq, aes(x = hours_since_sunset, y = frequency, col = paste(bat,sunset)))+
  # geom_point()+
  geom_smooth(se = FALSE, method = "lm", alpha = 0.2)+
  theme(legend.position = "none")

library(randomForest)

fit <- randomForest(factor(behavior) ~ frequency+amp+wfreq+wamp+rms+rms_filter+speed+amplitude,
                    data = Freq)#[Freq$behavior != "resting",])
fit
varImpPlot(fit)
Freq$behavior %>% table
