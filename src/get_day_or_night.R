get_day_or_night <- function(latitude, longitude, timestamp) {
  # assumes UTC time
  date <- as.Date(timestamp)
  times <- suncalc::getSunlightTimes(date, lat = latitude, lon = longitude)
  sunrise <- times$sunrise
  sunset <- times$sunset
  is_day <- timestamp >= sunrise & timestamp <= sunset
  return(ifelse(is_day, "day", "night"))
}