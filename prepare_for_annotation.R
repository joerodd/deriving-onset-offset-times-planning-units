options(scipen=5)
# load("mngu0_trimmed.RData")
# windows <- plyr::ldply(mngu0_trimmed,function(x){x$windows})
# mngu0 <- mngu0_trimmed 
# rm(mngu0_trimmed)
load("mngu0_untrimmed.RData")
windows <- plyr::ldply(mngu0,function(x){x$windows})
library(tidyverse)


# prepare windows

get_windows <- function(x){
  index <- x$index
  windows %>% filter(filename==mngu0[[index]]$ema$filename[1]) %>% mutate(row_number = row_number(),window_id = paste0(filename,"_",row_number,"_", window))
}

target_windows <- data.frame(index = seq(1,length(mngu0),1)) %>% group_by(index) %>% do(get_windows(.))

# prepare raw ema


prepare_raw_ema <- function(rr){
  mngu0[[rr$index]]$ema %>% filter(time > rr$interesting_window_start,time < rr$interesting_window_end)
}

raw_ema <- target_windows %>% group_by(window,window_id,index) %>% do(prepare_raw_ema(.)) %>% ungroup()


# pcas

pca_north_west <- function(data,sensor,...){
  pca <- prcomp(data,retx = F,...)
  
  # PC1
  if(abs(pca$rotation["y",1]) > abs(pca$rotation["z",1])){ # closer to y axis than to z axis
    if(pca$rotation["y",1]<0){ # should be rotated
      cat(paste("flipped PC1 y axis for ", sensor,"\n"))
      pca$rotation[,"PC1"] <- pca$rotation[,"PC1"] %*% matrix(c(-1,0,0,-1),ncol=2)
      # pca$x[,"PC1"] <- pca$x[,"PC1"] * -1
    }
  } else if(abs(pca$rotation["y",1]) < abs(pca$rotation["z",1])){ # closter to z axis than to y axis
    if(pca$rotation["z",1]<0){ # should be rotated
      cat(paste("flipped PC1 z axis for ", sensor,"\n"))
      pca$rotation[,"PC1"] <- pca$rotation[,"PC1"] %*% matrix(c(-1,0,0,-1),ncol=2)
      # pca$x[,"PC1"] <- pca$x[,"PC1"] * -1
    }
  }
  
  # PC2
  if(abs(pca$rotation["y",2]) > abs(pca$rotation["z",2])){ # closer to y axis than to z axis
    if(pca$rotation["y",2]<0){ # should be rotated
      cat(paste("flipped PC2 y axis for ", sensor,"\n"))
      pca$rotation[,"PC2"] <- pca$rotation[,"PC2"] %*% matrix(c(-1,0,0,-1),ncol=2)
      # pca$x[,"PC2"] <- pca$x[,"PC2"] * -1
    }
  } else if(abs(pca$rotation["y",2]) < abs(pca$rotation["z",2])){ # closter to z axis than to y axis
    if(pca$rotation["z",2]<0){ # should be rotated
      cat(paste("flipped PC2 z axis for ", sensor,"\n"))
      pca$rotation[,"PC2"] <- pca$rotation[,"PC2"] %*% matrix(c(-1,0,0,-1),ncol=2)
      # pca$x[,"PC2"] <- pca$x[,"PC2"] * -1
    }
  }
  return(pca)
}

TONGUE_pca <- pca_north_west(na.omit(bind_cols(ema%>% filter(sensor=="T1") %>% select(y,z),ema %>% filter(sensor=="T2") %>% select(y,z),ema %>% filter(sensor=="T3") %>% select(y,z))),sensor = "TONGUE",scale. = T)
LIPS_pca <- pca_north_west(na.omit(bind_cols(raw_ema%>% filter(sensor=="UL") %>% select(y,z),raw_ema %>% filter(sensor=="LL") %>% select(y,z))),sensor = "LIPS",scale. = T)

T1_pca_2 <- prcomp(raw_ema %>% filter(sensor=="T1") %>% select(y,z) %>% na.omit(),scale. = F)

T1_pca <- pca_north_west(raw_ema %>% filter(sensor=="T1") %>% select(y,z) %>% na.omit(),sensor = "T1",scale. = F)
T2_pca <- pca_north_west(raw_ema %>% filter(sensor=="T2") %>% select(y,z) %>% na.omit(),sensor = "T2",scale. = F)
T3_pca <- pca_north_west(raw_ema %>% filter(sensor=="T3") %>% select(y,z) %>% na.omit(),sensor = "T3",scale. = F)
JAW_pca <- pca_north_west(raw_ema %>% filter(sensor=="JAW") %>% select(y,z) %>% na.omit(),sensor = "JAW",scale. = F)
UL_pca <- pca_north_west(raw_ema %>% filter(sensor=="UL") %>% select(y,z) %>% na.omit(),sensor = "UL",scale. = F)
LL_pca <- pca_north_west(raw_ema %>% filter(sensor=="LL") %>% select(y,z) %>% na.omit(),sensor = "LL",scale. = F)


PCA_summaries <- bind_rows(as.data.frame(summary(TONGUE_pca)$importance) %>% mutate(sensor="TONGUE",metric=c("Standard deviation","Proportion of Variance","Cumulative Proportion")),
                           as.data.frame(summary(T1_pca)$importance) %>% mutate(sensor="T1",metric=c("Standard deviation","Proportion of Variance","Cumulative Proportion")),
                           as.data.frame(summary(T3_pca)$importance) %>% mutate(sensor="T2",metric=c("Standard deviation","Proportion of Variance","Cumulative Proportion")),
                           as.data.frame(summary(JAW_pca)$importance) %>% mutate(sensor="JAW",metric=c("Standard deviation","Proportion of Variance","Cumulative Proportion")),
                           as.data.frame(summary(UL_pca)$importance) %>% mutate(sensor="UL",metric=c("Standard deviation","Proportion of Variance","Cumulative Proportion")),
                           as.data.frame(summary(LL_pca)$importance) %>% mutate(sensor="LL",metric=c("Standard deviation","Proportion of Variance","Cumulative Proportion")),
                           as.data.frame(summary(LIPS_pca)$importance) %>% mutate(sensor="LIPS",metric=c("Standard deviation","Proportion of Variance","Cumulative Proportion")))

PCA_summaries %>% filter(metric=="Proportion of Variance")

# apply pcas to data -------------------------------------------------------

apply_PCAS <- function(index){
  
  w <- mngu0[[index]]$ema
  
  y <- list()
  y$labels <- mngu0[[index]]$labels
  y$wf <- mngu0[[index]]$wf
  
  y$ema_PCs <- bind_rows(
    bind_cols(w%>% filter(sensor=="T1") %>% select(time,file,filename,sensor,y,z,yz_delta,yz_deltadelta),
              as.data.frame(predict(T1_pca, w%>% filter(sensor=="T1") %>% select(y,z))) %>% select(PC1,PC2)),
    
    bind_cols(w%>% filter(sensor=="T2") %>% select(time,file,filename,sensor,x,y,z,yz_delta,yz_deltadelta),
              as.data.frame(predict(T2_pca, w%>% filter(sensor=="T2") %>% select(y,z))) %>% select(PC1,PC2)),
    
    bind_cols(w%>% filter(sensor=="T3") %>% select(time,file,filename,sensor,x,y,z,yz_delta,yz_deltadelta),
              as.data.frame(predict(T3_pca, w%>% filter(sensor=="T3") %>% select(y,z))) %>% select(PC1,PC2)),
    
    bind_cols(w%>% filter(sensor=="JAW") %>% select(time,file,filename,sensor,x,y,z,yz_delta,yz_deltadelta),
              as.data.frame(predict(JAW_pca, w%>% filter(sensor=="JAW") %>% select(y,z))) %>% select(PC1,PC2)),
    
    bind_cols(w%>% filter(sensor=="UL") %>% select(time,file,filename,sensor,x,y,z,yz_delta,yz_deltadelta),
              as.data.frame(predict(UL_pca, w%>% filter(sensor=="UL") %>% select(y,z))) %>% select(PC1,PC2)),
    
    bind_cols(w%>% filter(sensor=="LL") %>% select(time,file,filename,sensor,x,y,z,yz_delta,yz_deltadelta),
              as.data.frame(predict(LL_pca, w%>% filter(sensor=="LL") %>% select(y,z))) %>% select(PC1,PC2))
  )
  
  return(y)
}

mngu0_for_annotation <- lapply(seq(1,length(mngu0),1),apply_PCAS)


# functions for plots ------------------------------------------------------

prepare_ema <- function(rr){
  mngu0_for_annotation[[rr$index]]$ema_PCs %>% filter(time > rr$interesting_window_start - 100,time < rr$interesting_window_end + 100)
}

prepare_ema_wf <- function(rr){
  bind_rows(mngu0_for_annotation[[rr$index]]$ema_PCs %>% filter(time > rr$interesting_window_start - 100,time < rr$interesting_window_end + 100),
            mngu0_for_annotation[[rr$index]]$wf %>% filter(time > rr$interesting_window_start - 100,time < rr$interesting_window_end + 100) %>% mutate(sensor="waveform"))
}

prepare_labels <- function(rr){
  mngu0_for_annotation[[rr$index]]$labels %>% filter((end_time - dur > rr$interesting_window_start - 30&end_time -dur < rr$interesting_window_end + 30)
                                                     |
                                                       (end_time > rr$interesting_window_start - 30&end_time < rr$interesting_window_end + 30))
}

rotation_calc <- function(pca,target_sensor,suffix=""){
  g <- raw_ema %>% group_by(sensor) %>% filter(sensor==target_sensor)
  as.data.frame(matrix(c(1,0,0,1),ncol = 2) %*% t(pca$rotation)*-1) %>% mutate(y_end = mean(g$y,na.rm=T) - y,
                                                                            z_end = mean(g$z,na.rm=T) - z,
                                                                            y_origin=mean(g$y,na.rm=T),
                                                                            z_origin=mean(g$z,na.rm=T),PC=c("PC1","PC2"))
}


data_for_pcas_key <- bind_rows(rotation_calc(T1_pca,"T1") %>% mutate(sensor="T1"),
                               rotation_calc(T2_pca,"T2") %>% mutate(sensor="T2"),
                               rotation_calc(T3_pca,"T3") %>% mutate(sensor="T3"),
                               rotation_calc(UL_pca,"UL") %>% mutate(sensor="UL"),
                               rotation_calc(LL_pca,"LL") %>% mutate(sensor="LL"),
                               rotation_calc(JAW_pca,"JAW") %>% mutate(sensor="JAW"))

ema_for_rugplots <- target_windows %>% rowwise() %>% do(prepare_ema(.)) %>% gather("PC","value",PC1,PC2) %>% group_by(PC,sensor) %>%
  summarise(max=max(value,na.rm = T),min=min(value,na.rm = T)) %>% gather("metric","value",max,min)

# save functions for annotation -------------------------------------------

save(prepare_ema,prepare_ema_wf,prepare_labels,rotation_calc,file = "annotate/functions.Rdata")

# save data for annotation ------------------------------------------------

target_windows_annotation <- target_windows %>% mutate(vowel_stability_onset=as.numeric(NA),vowel_stability_offset=as.numeric(NA),vowel_stability_sensor=as.numeric(NA),
                                                       consonant_stability_onset=as.numeric(NA),consonant_stability_offset=as.numeric(NA),consonant_stability_sensor=as.numeric(NA))
save(target_windows_annotation,file = "annotate/data_in_progress.Rdata")

save(data_for_pcas_key,T1_pca,T2_pca,T3_pca,LL_pca,UL_pca,JAW_pca,mngu0_for_annotation,ema_for_rugplots,file = "annotate/static_data.Rdata")
