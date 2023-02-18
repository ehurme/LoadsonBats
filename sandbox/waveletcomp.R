library(WaveletComp)
1000/25
x = periodic.series(start.period = 10, length = 20000)
x = x + 0.2*rnorm(length(x)) # add some noise
plot(x, type = "l")



my.data <- data.frame(x = x, date = (1:length(x)/50))
levels <- 1000
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1,
                        dj = 1/levels,
                        lowerPeriod = 8,
                        upperPeriod = 16,
                        make.pval = TRUE,
                        n.sim = 10)
wt.image(my.w, color.key = "quantile", n.levels = levels, col.ridge = "purple",
         legend.params = list(lab = "wavelet power levels", mar = 4.7))
maximum.level = 1.001*max(my.w$Power.avg, my.w$Power.avg)
wt.avg(my.w, maximum.level = maximum.level)

image(my.w$Power)
my_freq <- my.w$Period[apply(my.w$Power,2,which.max)]
plot(my_freq)
# remove points outside the coi
cois <- which(my.w$coi.1 %in% 1:40)
my_coi <- 2^my.w$coi.2[cois]

idx <- which(my_freq < my_coi)
plot(my_freq[idx])


plot(my.w$Period[apply(my.w$Ridge,2,which.max)], pch = 16, col = rgb(0,0,0,.1))
lines(WT$coi.1, 2^WT$coi.2)
points(idx, my_freq[idx])

hist(my.w$Period[apply(my.w$Ridge,2,which.max)][idx])
summary(my.w$Period[apply(my.w$Ridge,2,which.max)][idx])
summary(my_freq[idx])

#### real data
##### Cyprus

w <- Wave(left = db$z, samp.rate = sampling_rate, bit = 16)
wf <- ffilter(w, f= sampling_rate, from = min_freq, to = max_freq, bandpass = TRUE)
spectro(wf, f = sampling_rate, wl = 16, ovlp = 75, fastdisp = TRUE)
# how to filter out points where bats are not flying?

# quick fix - calculate the range of z values
db %>% group_by(burst) %>%
  dplyr::summarise(minz = min(z),
            maxz = max(z),
            meanz = mean(z)) -> db_sum
db_sum$minz %>% plot
db_sum$maxz %>% plot
db_sum$meanz %>% plot
db_idx <- which(db_sum$maxz > 1)

cfreq <- data.frame()
bursts <- unique(db$burst)
i = 210
for(i in db_idx){
  b <- i
  z1 <- db$z[db$burst == b]
  t1 <- db$time[db$burst == b]
  # 100,
  layout(1)
  plot(t1,z1, type = "o")

  my.dataC <- data.frame(x = z1, date = t1)
  levels <- 200
  my.wC <- analyze.wavelet(my.dataC, "x",
                          loess.span = 0,
                          dt = 1,
                          dj = 1/levels,
                          lowerPeriod = 4,
                          upperPeriod = 16,
                          make.pval = TRUE,
                          n.sim = 10)
  wt.image(my.wC, color.key = "quantile", n.levels = levels, col.ridge = "purple",
           legend.params = list(lab = "wavelet power levels", mar = 4.7))
  my_freq <- my.wC$Period[apply(my.wC$Power,2,which.max)]
  # plot(my_freq)
  # remove points outside the coi
  cois <- which(my.wC$coi.1 %in% 1:40)
  my_coi <- 2^my.wC$coi.2[cois]

  idx <- which(my_freq < my_coi)
  # plot(my_freq[idx])


  plot(my.wC$Period[apply(my.wC$Ridge,2,which.max)], pch = 16, col = rgb(0,0,0,.1))
  lines(my.wC$coi.1, 2^my.wC$coi.2)
  points(idx, my_freq[idx])

  # hist(my.wC$Period[apply(my.wC$Ridge,2,which.max)][idx])
  # summary(my.wC$Period[apply(my.wC$Ridge,2,which.max)][idx])
  summary(my_freq[idx])
  c <- data.frame(freqs = my_freq[idx], burst = b, time = t1[1])
  cfreq <- rbind(cfreq, c)
}

plot(cfreq$time, cfreq$freqs, col = rgb(0,0,0,.1))




my.data1 = data.frame(x = db$z[df_stable$from[jj]:df_stable$to[jj]][20001:40000])
                      #date = db$time[df_stable$from[jj]:df_stable$to[jj]][20001:40000])
plot(my.data1$x, type = "l")
levels = 50
my.w1 <- analyze.wavelet(my.data1, "x",
                         loess.span = 0,
                         dt = 1,
                         dj = 1/levels,
                         lowerPeriod = 2,
                         upperPeriod = 8,
                         make.pval = TRUE,
                         n.sim = 10)
wt.image(my.w1, color.key = "quantile", n.levels = levels, legend.params = list(lab = "wavelet power levels", mar = 4.7))
hist(my.w1$Period[apply(my.w1$Ridge,2,which.max)], breaks = 200)
summary(my.w1$Period[apply(my.w1$Ridge,2,which.max)])
