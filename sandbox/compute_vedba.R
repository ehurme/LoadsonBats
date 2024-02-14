compute_vedba <- function(acc_x, acc_y, acc_z) {
  # Check for invalid arguments
  if (length(acc_x) == 0 || length(acc_y) == 0 || length(acc_z) == 0 || length(acc_x) != length(acc_y) || length(acc_x) != length(acc_z)) {
    stop("Invalid argument error")
  }

  # Compute mean acceleration
  mean_x <- mean(acc_x)
  mean_y <- mean(acc_y)
  mean_z <- mean(acc_z)

  # Compute absolute difference between acceleration and mean
  diff_x <- abs(acc_x - mean_x)
  diff_y <- abs(acc_y - mean_y)
  diff_z <- abs(acc_z - mean_z)

  # Compute VeDBA
  vedba <- sum(sqrt(diff_x^2 + diff_y^2 + diff_z^2))

  # Cap VeDBA at max value (if necessary)
  COMPUTE_VEDBA_MAX <- 1e9  # Example max value, adjust as needed
  if (vedba >= COMPUTE_VEDBA_MAX) {
    vedba <- COMPUTE_VEDBA_MAX
  }

  return(vedba)
}

# with(df[16501:16527,], compute_vedba(accX_mg, accY_mg, accZ_mg))
# with(df[16521:16547,], compute_vedba(accX_mg, accY_mg, accZ_mg))
# with(df[21:47,], compute_vedba(accX_mg, accY_mg, accZ_mg))
# with(df[1:27,], compute_vedba(accX_mg, accY_mg, accZ_mg))
#

files <- "C:/Users/ehumre/Dropbox/MPI/Wingbeat/Arizona/flightcage/20230610/Calibration/CalibrationTag1.csv"
df <- fread(files[1])
# resample to match sigfox data
idx <- seq(1, nrow(df), by = 2)
df <- df[idx,]
with(df[1801:1827,], compute_vedba(accX_mg, accY_mg, accZ_mg))
abline(v = 1800)
