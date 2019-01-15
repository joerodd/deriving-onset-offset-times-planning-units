library(tidyverse)

load("annotation_data_20180821.Rdata")
load("mngu0_untrimmed.RData")

annotated <- target_windows_annotation %>% ungroup()  %>% mutate(vowel = str_split(as.character(window)," ",simplify = T)[,1],
                                                                 consonant = str_split(as.character(window)," ",simplify = T)[,2]) %>%
  filter(!is.na(vowel_stability_onset)&!is.na(consonant_stability_onset))


annotated <- annotated %>% group_by(index,window,
                                    interesting_window_end,interesting_window_end,
                                    row_number,window_id,
                                    vowel_stability_onset,vowel_stability_offset,vowel_stability_sensor,
                                    consonant_stability_onset,consonant_stability_offset,consonant_stability_sensor,
                                    b,vowel,consonant) %>% do(retrieve_acoustic_onsets_offsets(.)) %>% ungroup()


do_rle_overlap_begin_end_times <- function(euc_object){
  # out_file = paste0("griderrors/euclid-worker",Sys.getpid(),".out")
  g <- rle(euc_object$difference_is_positive)$lengths
  longest_run = max(g[2:length(g)-1]) # like this to exclude runs that cross the window boundary
  euc_object$longest_run <- NA
  if(length(which(g==longest_run))==1){ #if valid longest run
    longest_run_start = sum(g[1:which(g==longest_run)-1])+1
    longest_run_end = longest_run_start + longest_run - 1
    euc_object$longest_run <- NA
    euc_object$longest_run[longest_run_start:longest_run_end] <- euc_object$interpolated_positive_difference[longest_run_start:longest_run_end]
    euc_object
    return(data.frame(acoustic_consonant_onset=euc_object$time[longest_run_start],
                      acoustic_vowel_offset=euc_object$time[longest_run_end],
                      overlap_duration=euc_object$time[longest_run_end]-euc_object$time[longest_run_start]))
  } else {
    return(data.frame(acoustic_consonant_onset=NA,
                      acoustic_vowel_offset=NA,
                      overlap_duration=NA))
  }
}


retrieve_acoustic_onsets_offsets <- function(x){
  
  mngu0[[x$index]]$euc %>% filter(time >= x$interesting_window_start,time <= x$interesting_window_end) %>% do_rle_overlap_begin_end_times()
  
}


save(annotated,file="annotation_data_20180821_with_acoustic.Rdata")