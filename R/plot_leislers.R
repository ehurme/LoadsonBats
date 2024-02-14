path = "../../../../Dropbox/MPI/Wingbeat/Belgium23/Nyctalus_leisleri_fleatag/"
files <- list.files(path, pattern = "*.csv", recursive = TRUE)

sampling_rate <- 105
window_length = sampling_rate*8
burst_length = 156

j = 12
pdf(file = paste0(path, "leisler_acc.pdf"))
for(j in 1:length(files)){
  df <- fread(paste0(path, files[j]))
  # add data by burst
  max_burst <- as.numeric(df$burstCount[nrow(df)])
  d <- data.frame(time = 1:((max_burst+1) * burst_length)/sampling_rate,
                  burst = rep(0:max_burst, each = burst_length),
                  x = NA, y = NA, z = NA)

  bursts <- unique(df$burstCount) %>% as.numeric() %>% na.omit()
  bursts <- bursts[bursts >= 0 & bursts <= max_burst]
  try({
    for(i in 1:length(bursts)){
      idx <- which(df$burstCount == bursts[i])
      didx <- which(d$burst == bursts[i])
      if(length(idx) != length(didx)){
        d$x[didx[1:length(idx)]] <- df$accX_mg[idx] %>% as.numeric
        d$y[didx[1:length(idx)]] <- df$accY_mg[idx] %>% as.numeric
        d$z[didx[1:length(idx)]] <- df$accZ_mg[idx] %>% as.numeric
      }
      if(length(idx) == length(didx)){
        d$x[didx] <- df$accX_mg[idx] %>% as.numeric
        d$y[didx] <- df$accY_mg[idx] %>% as.numeric
        d$z[didx] <- df$accZ_mg[idx] %>% as.numeric
      }
    }
  })

  summary(d)

  d$x[which(abs(d$x) > 8000)] <- 8000
  d$y[which(abs(d$y) > 8000)] <- 8000
  d$z[which(abs(d$z) > 8000)] <- 8000

  # plot data
  layout(1)
  with(d, #[9500:9700,],
       plot(time/60, x/1000, type = "l",
            xlab = "time (min)", ylab = "Acceleration (Gs)",
            ylim = c(-8,8))
       )
  lines(d$time/60, d$y/1000, col = rgb(0,1,0,.4))
  lines(d$time/60, d$z/1000, col = rgb(1,0,0,.4))
  title(main = files[j],  line = 0.5, cex.main = 0.8)

  summary(d)
}
dev.off()


