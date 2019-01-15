library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tuneR)
library(seewave)
options(scipen = 9999)



mfcc <- runMFCC_notrim("audio_data/mngu0_s1_0001.wav")
# save(mfcc,file = "mfcc.RData")
euc <- calc_euc(mfcc)
save(mfcc,euc,file = "euc.RData")