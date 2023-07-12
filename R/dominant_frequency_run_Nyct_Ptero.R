
acc_files <- list.files("C:/Users/Edward/movebank_downloads/study-484630593", pattern = "acc.csv",
                        full.names = TRUE, recursive = TRUE)
gps_files <- list.files("C:/Users/Edward/movebank_downloads/study-484630593", pattern = "gps.csv",
                        full.names = TRUE, recursive = TRUE)
acc_files
gps_files
i = 10
for(i in 1:length(acc_files)){
  dominant_freq(df = fread(acc_files[i]),
                GPS = fread(gps_files[i]),
              save_path = "./../../../ownCloud/Firetail/Nyctaluslasiopterus/GPA-10_8147_S1/Wingbeats/PCA/",
              tag_id = substr(acc_files[i], 56, 64),
              Burst = FALSE,
              PCA = TRUE,
              sampling_rate = 100, # P. hastatus
              min_freq = 2,
              max_freq = 8,
              # wavelet = FALSE,
              saved_cores = 8,
              gps = TRUE,
              min_seg_duration = 10,
              dfreq_threshold = 20,
              use_FFT = TRUE,
              # use_wavelet = FALSE,
              tag_type = "eObs", # "wildfi"
              firetail = FALSE,
              firetail_filter = FALSE,
              solar_time = TRUE,
              save_files = TRUE,
              location = "Portugal",
              var_seg = TRUE,
              sd_adjust = 1)

}
nchar("C:/Users/Edward/movebank_downloads/study-484630593/tag-484637139")




files <- list.files(path = "./../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/",
                    pattern = ".csv", full.names = TRUE)

for(i in 4:length(files)){
  dominant_freq_wvlt(df = fread(files[i]) %>% clean_names(),
                   save_path = "./../../../ownCloud/Firetail/Pteropuslylei/Model_tag_2268/Wingbeats/PCA/",
                   tag_id = substr(files[i], 63, 66),
                   PCA = TRUE,
                   sampling_rate = 18.74,
                   min_freq = 2,
                   max_freq = 8,
                   wavelet = TRUE,
                   saved_cores = 8,
                   gps = TRUE,
                   FFT = TRUE,
                   tag_type = "eObs", # "wildfi"
                   firetail = TRUE,
                   firetail_filter = "c",
                   FFT_amp_thresh = -0.2,
                   solar_time = TRUE,
                   save_files = TRUE,
                   levels = 1000)
}


