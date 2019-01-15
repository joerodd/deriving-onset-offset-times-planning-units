library(tidyverse)
load("intermediate_data/annotation_data_20180821_with_acoustic.Rdata")
load("grid_results/grid_result_multi_within_definitiveA.RData")



# ----------- correlation -----

joined_grid_result <- grid_result %>% filter(score_function == "multi_within",interpolation_n==40)

load("grid_results/grid_result_multi_within_definitiveB.RData")

joined_grid_result <- bind_rows(joined_grid_result,
                                grid_result %>% filter(score_function == "multi_within",interpolation_n==40))
  
get_multidimensional_peak <- function(vowel_offset,consonant_onset,scores,index=NA){
    # if(!is.na(index[1])){
    #   cat(index[1])
    # }
    if(length(scores[which(scores<1)]) > 5){
      scores = 1 - scores[which(scores<1)]
      scores <-  scores /  sum(scores)
      vowel_offset <- vowel_offset[which(scores<1)]
      consonant_onset <- consonant_onset[which(scores<1)]
      
      ppp <- kde(matrix(c(vowel_offset,consonant_onset),ncol=2),w = scores)
      
      maxpos = which(ppp$estimate==max(ppp$estimate))
      maxrow = maxpos%%151
      maxcol = 1 + (maxpos-(maxpos%%151))/151
      
      return(list(vowel_offset = ppp$eval.points[[1]][maxrow],
                  consonant_onset = ppp$eval.points[[2]][maxcol],
                  dim1 = dim(ppp$estimate)[1],
                  dim2 = dim(ppp$estimate)[2]))
    } else{
      return(list(vowel_offset = NA,
                  consonant_onset = NA,
                  dim1 = NA,
                  dim2 = NA
      ))  
    }
  }

indexcatter <- function(index){
  cat(index[1])
  cat(" ")
  return(index[1])
}

rsq <- function (x, y){
  summary(lm(y~x))$r.squared
  # cor(x, y,use = "complete.obs") ^ 2
} 

library(multidplyr)
cl <- multidplyr::create_cluster(4)
cluster_library(cl,c("tidyverse","ks"))
cluster_copy(cl,get_multidimensional_peak)

multi_summary <- joined_grid_result %>% partition(cluster = cl,index,window,window_id,interpolation_n,score_function) %>% summarise(
  articulatory_consonant_onset = get_multidimensional_peak(vowel_offset,consonant_onset,score,score_function)$consonant_onset,
  articulatory_vowel_offset = get_multidimensional_peak(vowel_offset,consonant_onset,score)$vowel_offset) %>%
  mutate(articulatory_overlap_duration = articulatory_vowel_offset-articulatory_consonant_onset) %>% collect()

comparison_multi <- left_join(multi_summary, annotated) %>% 
  mutate(articulatory_consonant_onset = articulatory_consonant_onset - interesting_window_end,
         articulatory_vowel_offset = articulatory_vowel_offset - interesting_window_end,
         acoustic_consonant_onset  = acoustic_consonant_onset - interesting_window_end,
         acoustic_vowel_offset = acoustic_vowel_offset - interesting_window_end,
         peak_detection="multidimensional peak detection")


collapsed_comparison <- comparison_multi %>% gather(event,time,articulatory_consonant_onset,articulatory_vowel_offset,acoustic_consonant_onset,acoustic_vowel_offset) %>% mutate(
  metric=recode(event,articulatory_consonant_onset="articulatory",articulatory_vowel_offset="articulatory",
                acoustic_consonant_onset="acoustic",acoustic_vowel_offset="acoustic"),
  event=recode(event,articulatory_consonant_onset="consonant_onset",articulatory_vowel_offset="vowel_offset",
               acoustic_consonant_onset="consonant_onset",acoustic_vowel_offset="vowel_offset")) %>% spread(metric,time)

collapsed_comparison <- collapsed_comparison %>% ungroup()

annotations <- collapsed_comparison %>% group_by(interpolation_n,score_function,peak_detection) %>% summarise(rsq = rsq(articulatory,acoustic))

library(boot)

rsq2 <- function(original,sample){original[sample,] %>% summarise(rsq = rsq(articulatory,acoustic)) %>% unlist()}

collapsed_comparison %>% rsq2()

boot.ci(boot(data = collapsed_comparison,statistic = rsq2,R = 10000))

ff <- Hmisc::rcorr(collapsed_comparison$acoustic, collapsed_comparison$articulatory,use = "complete.obs")
rsq(collapsed_comparison$acoustic, collapsed_comparison$articulatory)

correlationLM <- lm(articulatory~acoustic,data=collapsed_comparison)

save(collapsed_comparison,annotations,correlationLM,file="Manuscript/definitive_correlation_result.RData")
