# fleatag
library(pacman)

p_load(data.table, seewave, tuneR)

# load data
df <- fread("../../FleaTagLogs/20230607/Efuscus_bat1/20230607_Ef_bat1_tag3029_lowestweight.csv", header = T, fill = T) #13.25 min
df <- fread("../../FleaTagLogs/20230607/Efuscus_bat2/20230607_Ef_bat2_tag3016_lowestweight.csv", header = T, fill = T) #10.6 min
df <- fread("../../FleaTagLogs/20230608/Mthy_bat1_t1_no_weight/20230608_Mthy_bat1_t1_lowestweight_tag3029.csv", header = T, fill = T) # 10.57 min
df <- fread("../../FleaTagLogs/20230608/Mthy_bat1_t2_weight/20230608_Mthy_bat1_t2_highestweight_tag3029.csv", header = T, fill = T) # 9.28 min
df <- fread("../../FleaTagLogs/20230610/Myotis_californicus_bat1_lightestweight/Myotis_californicus_bat1_lightestweight.csv", header = T, fill = T) # 9.42 min

df <- fread("../../FleaTagLogs/20230611/Mthy_bat1_trial2/Mthy_bat1_trial2_weightedtag.csv", header = T, fill = T) # 10.3 min
df <- fread("../../FleaTagLogs/20230611/Mthy_bat1_trial3/ACC_20230611_Mthy_bat1_trial3.csv", header = T, fill = T) # 10.3 min
df <- fread("../../FleaTagLogs/20230611/Mthy_bat1_trial4/ACC_20230611_Mthy_bat1_trial4.csv", header = T, fill = T) # 10.3 min

df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230614/", header = T, fill = T) # 10.3 min

df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230618/Efus_bat1_trial2/ACC_20230618_Efus_bat1_trial2.csv", header = T, fill = T) # 10.3 min

df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230622/Mthy_bat54_trial2/ACC_20230622_Mthy_bat54.csv", header = T, fill = T)
df <- fread("../../../Dropbox/MPI/Wingbeat/Arizona/flightcage/20230623/Mthy_bat54_trial4/ACC_Mthy_bat54_trial4.csv", header = T, fill = T)
# remove blank rows
df <- df[!is.na(df$burstCount),]

# df1 <- df1[!is.na(df1$burstCount),]
# df2 <- df2[!is.na(df2$burstCount),]

sampling_rate <- 54

# add tag time
df$time <- (1:nrow(df))/sampling_rate
max(df$time)/60

# df1$time <- (1:nrow(df1))/54
# df2$time <- (1:nrow(df2))/54

dev.off()


plot(df$time, df$ColorSensIR_cnt)



with(df,#[1:1000,],#[9500:9700,],
     plot(time/60, accX_mg, type = "l"))
lines(df$time/60, df$accY_mg, col = 2)
lines(df$time/60, df$accZ_mg, col = 3)



w <- tuneR::Wave(left = df$accZ_mg, samp.rate = sampling_rate, bit = 16)
wl <- round(length(w)/20, -1)
duration(w)
plot(w)
spec <- meanspec(w, wl)
fpeaks(spec, nmax = 10)
spectro(w, wl)
dev.off()

flapping <- ffilter(w, f= sampling_rate, from = 7, to = 20, bandpass = TRUE,
                    wl = wl)
spectro(flapping, f = sampling_rate)
layout(1)
spec <- meanspec(flapping, f = sampling_rate, wl)

wf <- {}
wf <- ffilter(w, f= sampling_rate, from = 0, to = 2, bandpass = TRUE, wl)
spectro(wf, f = sampling_rate)
layout(1)

plot(wf, type = "l")

plot(w, col = "gray")
lines((1:nrow(df))/sampling_rate, wf*1000, lwd = 3)

rms(df$accZ_mg)
rms(wf*100)




