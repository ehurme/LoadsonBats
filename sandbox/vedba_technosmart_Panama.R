files <- "../waveanalysis/Filename_selection_split#5_0.txt"
sampling_rate <- 25
window_length = sampling_rate*2

df <- fread(files[1])
plot(df$Acc_z)
## resample to match sigfox data
# idx <- seq(1, nrow(df), by = 2)
# df <- df[idx,]
sx <- zoo::rollmean(df$Acc_x, k = window_length, na.pad = TRUE)
sy <- zoo::rollmean(df$Acc_y, k = window_length, na.pad = TRUE)
sz <- zoo::rollmean(df$Acc_z, k = window_length, na.pad = TRUE)
plot(sz)

ax <- df$Acc_x - sx
ay <- df$Acc_y - sy
az <- df$Acc_z - sz
hist(az)

ODBA <- abs(ax)+abs(ay)+abs(az)
VeDBA <- sqrt(ax^2 + ay^2 + az^2)
range(VeDBA, na.rm = T)
plot(ODBA, VeDBA)

av <- zoo::rollmean(VeDBA, k = sampling_rate, na.pad = TRUE)

plot(az, ylim = c(-4,10))
points(ODBA, col = rgb(1,0,0,.1))
points(VeDBA, col = rgb(0,1,0,.1))

points(av, col = rgb(0,0,1,.1))

hist(av, breaks = 1000)

sv <- zoo::rollsum(VeDBA, k = sampling_rate, na.pad = TRUE)
plot(sv, col = rgb(0,0,1,.1))

hist(VeDBA, breaks = 100)
hist(av, breaks = 100)
acc <- data.frame(time = 1:length(VeDBA)/sampling_rate, ODBA, VeDBA, av)
ggplot(acc, aes(time, av))+geom_point(alpha = 0.1)+geom_smooth()