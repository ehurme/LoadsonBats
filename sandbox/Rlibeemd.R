# https://github.com/helske/Rlibeemd
# devtools::install_github("helske/Rlibeemd")
library(Rlibeemd)
data(UKgas, package = "datasets")
imfs <- ceemdan(UKgas, ensemble_size = 1000)
plot(imfs, main = "Five IMFs and residual extracted by CEEMDAN algorithm")
plot(imfs[,1])