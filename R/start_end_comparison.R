# look at wingbeat frequency from the start and end of the flight
library(pacman)
p_load(data.table, tidyverse, lubridate,
       ggplot2)

# load data
load("../../../ownCloud/Firetail/Acerodonjubatus/tag_1521/accelerateR/tag_1521.robj")
load("../../../ownCloud/Firetail/WhiteStor/Affenburg/")

ggplot(freqs, aes(x = time, y = frequency, size = amp, col = behavior))+
  geom_point()
ggplot(freqs, aes(x = time, y = speed, size = amp, col = behavior))+
  geom_point()

ggplot(freqs, #[freqs$rms > 100 & freqs$frequency > 2 & freqs$frequency < 5,],
       aes(x = rms, y = frequency, size = amp, col = behavior))+
  geom_point()

ggplot(freqs[freqs$rms > 100 & freqs$frequency > 2 & freqs$frequency < 5,],
       aes(x = hour(time), y = frequency))+
  geom_point()+
  geom_smooth(method = 'lm')+
  facet_wrap(~date(time))



# get sunset time

# get departure time

# find flight with high speed and low turn angle
## embc?

# look at corresponding ACC