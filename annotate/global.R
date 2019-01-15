library(shiny)
library(DT)
library(tidyverse)

load("static_data.Rdata")
load("data_in_progress.Rdata")
load("functions.Rdata")

windows_to_include <- data.frame(window = c("a m", "a S", "a v", "i m", "i s", "i v", "I v", "I S", "Q n", "V m", "V s"),
                                 vowel=c("a","a","a", "i", "i", "i", "I", "I", "Q", "V", "V"),
                                 consonant=c("m", "S", "v", "m", "s", "v", "v", "S", "n", "m", "s"),use = T)

target_windows_annotation <- target_windows_annotation %>% mutate(b = str_match(target_windows_annotation$window_id,"mngu0_s1_(\\d+)([a-z]?)_.*")[,3]!="")

target_windows_annotation <- target_windows_annotation %>% filter(window %in% windows_to_include$window,!b)



selected_row <- 1

vowel_stability_onset <- NA
consonant_stability_onset <- NA
vowel_stability_offset <- NA
consonant_stability_offset <- NA
vowel_stability_sensor <- NA
consonant_stability_sensor <- NA
