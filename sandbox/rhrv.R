# https://cran.r-project.org/web/packages/RHRV/vignettes/RHRV-quickstart.html
# http://rhrv.r-forge.r-project.org/tutorial/tutorial.pdf
# https://link.springer.com/book/10.1007/978-3-319-65355-6
# install.packages("RHRV", dependencies = TRUE)
library(pacman)
p_load(RHRV, data.table, tidyverse, dplyr, lubridate, POT)

# HRVData structure containing the heart beats
data("HRVData")
# HRVData structure storing the results of processing the
# heart beats: the beats have been filtered, interpolated, ...
data("HRVProcessedData")

HRVProcessedData$HR %>% plot(type = "l")

hrv.data  = CreateHRVData()
hrv.data = SetVerbose(hrv.data, TRUE )

hrv.data = BuildNIHR(HRVData = HRVData)

PlotNIHR(hrv.data, main = "niHR")

hrv.data = CreateTimeAnalysis(hrv.data, size = 300,
                              interval = 7.8125)

hrv.data = FilterNIHR(hrv.data)
hrv.data = SetVerbose(hrv.data,TRUE)
hrv.data = CreateTimeAnalysis(hrv.data,size=300,interval = 7.8125)

hrv.data = BuildNIHR(hrv.data)
hrv.data = FilterNIHR(hrv.data)
hrv.data = InterpolateNIHR (hrv.data, freqhr = 4)
hrv.data = CreateFreqAnalysis(hrv.data)
hrv.data = SetVerbose(hrv.data,TRUE)
# Note that it is not necessary to write the boundaries
# for the frequency bands, since they match
# the default values
hrv.data =
  CalculatePowerBand(hrv.data , indexFreqAnalysis = 1,
                     size = 300, shift = 30, type = "fourier",
                     ULFmin = 0, ULFmax = 0.03, VLFmin = 0.03, VLFmax = 0.05,
                     LFmin = 0.05, LFmax = 0.15, HFmin = 0.15, HFmax = 0.4 )
PlotPowerBand(hrv.data, indexFreqAnalysis = 1, ymax = 2000, ymaxratio = 1.7)

hrv.data = BuildNIHR(hrv.data)
hrv.data = FilterNIHR(hrv.data)
hrv.data = InterpolateNIHR (hrv.data, freqhr = 4)
hrv.data = CreateFreqAnalysis(hrv.data)
hrv.data = SetVerbose(hrv.data,TRUE)
# Note that it is not necessary to write the boundaries
# for the frequency bands, since they match the default values
hrv.data =
  CalculatePowerBand( hrv.data , indexFreqAnalysis = 1,
                      type = "wavelet", wavelet = "la8",
                      bandtolerance = 0.01, relative = FALSE,
                      ULFmin = 0, ULFmax = 0.03, VLFmin = 0.03, VLFmax = 0.05,
                      LFmin = 0.05, LFmax = 0.15, HFmin = 0.15, HFmax = 0.4 )


# create structure, load beats, filter and interpolate
hrv.data = CreateFreqAnalysis(hrv.data)
hrv.data = SetVerbose(hrv.data, FALSE)
# use freqAnalysis number 1 for perfoming
# Fourier analysis. This time, we do not
# write the band's boundaries
hrv.data = CalculatePowerBand(hrv.data , indexFreqAnalysis = 1,
                              size = 300, shift = 30, sizesp = 2048,
                              type = "fourier")
# use freqAnalysis number 2 for perfoming
# wavelet analysis. Note the indexFreqAnalysis = 2!!!
hrv.data = CreateFreqAnalysis(hrv.data)
hrv.data = CalculatePowerBand(hrv.data, indexFreqAnalysis= 2,
                              type = "wavelet", wavelet = "la8",
                              bandtolerance = 0.01, relative = FALSE)


# Plotting wavelet analysis
PlotPowerBand(hrv.data, indexFreqAnalysis = 2, ymax = 700, ymaxratio = 50)
