# fleatag

library(data.table)

# load data
df <- fread("../../FleaTagLogs/20230607/Efuscus_bat1/20230607_Ef_bat1_tag3029_lowestweight.csv", header = T, fill = T)
df <- fread("../../FleaTagLogs/20230607/Efuscus_bat2/20230607_Ef_bat2_tag3016_lowestweight.csv", header = T, fill = T)
df <- fread("../../FleaTagLogs/20230608/Mthy_bat1_t1_no_weight/20230608_Mthy_bat1_t1_lowestweight_tag3029.csv", header = T, fill = T)
df <- fread("../../FleaTagLogs/20230608/Mthy_bat1_t2_weight/20230608_Mthy_bat1_t2_highestweight_tag3029.csv", header = T, fill = T)

# remove blank rows
df <- df[!is.na(df$burstCount),]

df1 <- df1[!is.na(df1$burstCount),]
df2 <- df2[!is.na(df2$burstCount),]

# add tag time
df$time <- (1:nrow(df))/54
df1$time <- (1:nrow(df1))/54
df2$time <- (1:nrow(df2))/54

dev.off()

plot(df$ColorSensRed_cnt)
plot(df$time, df$ColorSensIR_cnt)



with(df,# [1:3000,],#[9500:9700,],
     plot(time/60, accX_mg, type = "l"))
lines(df$time/60, df$accY_mg, col = 2)
lines(df$time/60, df$accZ_mg, col = 3)


sampling_rate <- 54
w <- tuneR::Wave(left = df$accZ_mg, samp.rate = sampling_rate, bit = 16)
duration(w)
plot(w)
meanspec(w, wl = length(w)/20)
spectro(w, wl = length(w)/5)
dev.off()
wf <- ffilter(w, f= sampling_rate, from = 0, to = 8, bandpass = TRUE, wl = length(w)/5)
plot(w, col = "gray")
lines((1:nrow(df))/sampling_rate, wf*1000, lwd = 3)

rms(df$accZ_mg)
rms(wf*100)




x <- 1:100
y <- x^.5
diff(y)
plot(x,log(y))
