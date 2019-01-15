library(tidyverse)
library(HDInterval)
load("annotate/static_data.Rdata")
load("annotate/functions.Rdata")

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# identify targets --------------------------------------------------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# load("phase1data.Rdata")
# 
# phase1 <- target_windows_annotation %>% ungroup()  %>% mutate(vowel = str_split(as.character(window)," ",simplify = T)[,1],
#                                                                    consonant = str_split(as.character(window)," ",simplify = T)[,2]) %>%
#   filter(!is.na(vowel_stability_onset)&!is.na(consonant_stability_onset)) %>% mutate(b = str_match(window_id,"mngu0_s1_(\\d+)([a-z]?)_.*")[,3]!="")

load("intermediate_data/annotation_data_20180821.Rdata")

annotated <- target_windows_annotation %>% ungroup()  %>% mutate(vowel = str_split(as.character(window)," ",simplify = T)[,1],
                                                                consonant = str_split(as.character(window)," ",simplify = T)[,2]) %>%
  filter(!is.na(vowel_stability_onset)&!is.na(consonant_stability_onset))

annotated_for_targets <- annotated %>% unique() %>% filter(!b)

annotated_for_targets <- annotated_for_targets %>% filter(window %in% c("a m", "a S", "a v", "i m", "i s", "i v", "I v", "I S", "Q n", "V m", "V s"))

# retrieve ema and labels -------------------------------------------------

ema_PCs <- annotated_for_targets %>% group_by(window,window_id,index) %>% do(prepare_ema(.)) %>% ungroup()
labels <-  annotated_for_targets %>% group_by(window,window_id,index) %>% do(prepare_labels(.)) %>% ungroup()

good_consonants <- unlist(annotated_for_targets %>% group_by(consonant) %>% summarise(n=n()) %>% filter(n >= 8) %>% select(consonant))
good_vowels <- unlist(annotated_for_targets %>% group_by(vowel) %>% summarise(n=n()) %>% filter(n >= 8) %>% select(vowel))

sample_consonant <- annotated %>% filter(consonant %in% good_consonants) %>% select(window_id) %>% unlist()

sample_vowel <- annotated %>% filter(vowel %in% good_vowels)  %>% select(window_id) %>% unlist()


# calculate relative times ------------------------------------------------

ema_PCs_relative_consonant <- left_join(ema_PCs,annotated_for_targets) %>% mutate(average_consonant_time = (consonant_stability_offset + consonant_stability_onset)/2,
                                                                     time = time - average_consonant_time)

labels_relative_consonant <- left_join(labels,annotated_for_targets) %>% mutate(average_consonant_time = (consonant_stability_offset + consonant_stability_onset)/2,
                                                                   time = time - average_consonant_time,
                                                                   end_time = end_time - average_consonant_time)

r_relative_consonant <- annotated %>% mutate(average_consonant_time = (consonant_stability_offset + consonant_stability_onset)/2,
                                            vowel_stability_onset = vowel_stability_onset - average_consonant_time,
                                            vowel_stability_offset = vowel_stability_offset - average_consonant_time,
                                            consonant_stability_onset = consonant_stability_onset - average_consonant_time,
                                            consonant_stability_offset = consonant_stability_offset - average_consonant_time)

ema_PCs_relative_vowel <- left_join(ema_PCs,annotated_for_targets) %>% mutate(average_vowel_time = (vowel_stability_offset + vowel_stability_onset)/2,
                                                                 time = time - average_vowel_time)

labels_relative_vowel <- left_join(labels,annotated_for_targets) %>% mutate(average_vowel_time = (vowel_stability_offset + vowel_stability_onset)/2,
                                                               time = time - average_vowel_time,
                                                               end_time = end_time - average_vowel_time)

r_relative_vowel <- annotated %>% mutate(average_vowel_time = (vowel_stability_offset + vowel_stability_onset)/2,
                                         vowel_stability_onset = vowel_stability_onset - average_vowel_time,
                                         vowel_stability_offset = vowel_stability_offset - average_vowel_time,
                                         consonant_stability_onset = consonant_stability_onset - average_vowel_time,
                                         consonant_stability_offset = consonant_stability_offset - average_vowel_time)




# establish targets -------------------------------------------------------



vowel_targets <- ema_PCs_relative_vowel %>% filter(window_id %in% sample_vowel) %>% gather("PC","value",PC1,PC2) %>% filter(time>vowel_stability_onset-average_vowel_time,
                                                                                                                            time<vowel_stability_offset-average_vowel_time) %>%
  group_by(vowel,sensor,PC) %>% summarise(upper=quantile(value,0.95),
                                          target=quantile(value,0.5),
                                          lower=quantile(value,0.05))

consonant_targets <- ema_PCs_relative_consonant %>% filter(window_id %in% sample_consonant) %>% gather("PC","value",PC1,PC2) %>% filter(time>consonant_stability_onset-average_consonant_time,
                                                                                                                                        time<consonant_stability_offset-average_consonant_time) %>%
  group_by(consonant,sensor,PC) %>% summarise(upper=quantile(value,0.95),
                                              target=quantile(value,0.5),
                                              lower=quantile(value,0.05))

get_HDI <- function(x){
  ee <- as.data.frame(hdi(density(x$value),0.9,allowSplit=T)) %>% transmute(lower=begin,upper=end,difference=end-begin,HDI_index=row_number())
  
  if(nrow(ee)>1){
    ef <- ee %>% filter(difference>0.1)
    if(nrow(ef)<1){
      ee
    } else {
      ef
    }
  } else {
    ee
  }
  
}

consonant_targets <- ema_PCs_relative_consonant %>% filter(window_id %in% sample_consonant) %>% gather("PC","value",PC1,PC2) %>% filter(time>consonant_stability_onset-average_consonant_time,
                                                                                                                                            time<consonant_stability_offset-average_consonant_time) %>%
  group_by(consonant,sensor,PC) %>% do(get_HDI(.))
  
  
vowel_targets <- ema_PCs_relative_vowel %>% filter(window_id %in% sample_vowel) %>% gather("PC","value",PC1,PC2) %>% filter(time>vowel_stability_onset-average_vowel_time,
                                                                                                                            time<vowel_stability_offset-average_vowel_time) %>%
  group_by(vowel,sensor,PC) %>% do(get_HDI(.))

  


save(vowel_targets,consonant_targets,file="intermediate_data/targets.Rdata")



# plot of curves with targets ---------------------------------------------

ema_PCs_relative_consonant %>% filter(window_id %in% sample_consonant) %>% gather("PC","value",PC1,PC2) %>% ggplot()+
  # geom_line(aes(x=time,y=amplitude,group=window_id))+
  geom_rect(aes(xmin=vowel_stability_onset,ymin=-Inf,ymax=Inf,xmax=vowel_stability_offset),data=r_relative_consonant %>% filter(window_id %in% sample_consonant),fill="royalblue1",alpha=0.2)+
  geom_rect(aes(xmin=consonant_stability_onset,ymin=-Inf,ymax=Inf,xmax=consonant_stability_offset),data=r_relative_consonant %>% filter(window_id %in% sample_consonant),fill="palegreen2",alpha=0.2)+
  
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=lower,ymax=upper,fill=PC),data=consonant_targets,alpha=0.4)+
  
  geom_rug(aes(y=value,colour=PC),data=ema_for_rugplots %>% filter(sensor%in%c("LL","UL","JAW","T1","T2","T3")))+
  geom_line(aes(x=time,y=value,colour=PC,group=paste(window_id,PC)))+
  # geom_vline(aes(xintercept=end_time-dur),data=labels)+
  # geom_vline(aes(xintercept=end_time),data=labels)+
  geom_text(aes(x=end_time - (dur / 2),y=Inf,label=label),data=labels_relative_consonant %>% filter(window_id %in% sample_consonant),vjust=2,family="Gentium Plus")+
  facet_grid(sensor~consonant,scales = "free_y")


ema_PCs_relative_vowel %>% filter(window_id %in% sample_vowel) %>% gather("PC","value",PC1,PC2) %>% ggplot()+
  # geom_line(aes(x=time,y=amplitude,group=window_id))+
  geom_rect(aes(xmin=vowel_stability_onset,ymin=-Inf,ymax=Inf,xmax=vowel_stability_offset),data=r_relative_vowel %>% filter(window_id %in% sample_vowel),fill="royalblue1",alpha=0.2)+
  geom_rect(aes(xmin=consonant_stability_onset,ymin=-Inf,ymax=Inf,xmax=consonant_stability_offset),data=r_relative_vowel %>% filter(window_id %in% sample_vowel),fill="palegreen2",alpha=0.2)+
  
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=lower,ymax=upper,fill=PC),data=vowel_targets,alpha=0.4)+
  
  geom_rug(aes(y=value,colour=PC),data=ema_for_rugplots %>% filter(sensor%in%c("LL","UL","JAW","T1","T2","T3")))+
  geom_line(aes(x=time,y=value,colour=PC,group=paste(window_id,PC)),size=0.2)+
  geom_line(aes(x=time,y=value,colour=PC,group=paste(window_id,PC)),data=
              ema_PCs_relative_vowel %>% filter(window_id %in% sample_vowel) %>% gather("PC","value",PC1,PC2) %>% filter(time>vowel_stability_onset-average_vowel_time,
                                                                                                                         time<vowel_stability_offset-average_vowel_time))+
  # geom_vline(aes(xintercept=end_time-dur),data=labels)+
  # geom_vline(aes(xintercept=end_time),data=labels)+
  geom_text(aes(x=end_time - (dur / 2),y=Inf,label=label),data=labels_relative_vowel %>% filter(window_id %in% sample_vowel),vjust=2,family="Gentium Plus")+
  facet_grid(sensor~vowel,scales = "free_y")



# density plots -----------------------------------------------------------

library(ggridges)

ema_PCs_relative_consonant %>% filter(time>consonant_stability_onset-average_consonant_time,
                                      time<consonant_stability_offset-average_consonant_time,
                                      window_id %in% sample_consonant) %>% gather("PC","value",PC1,PC2) %>% ggplot()+
  # geom_line(aes(x=time,y=amplitude,group=window_id))+
  geom_density_ridges(aes(y=PC,x=value,colour=PC))+
  geom_rect(aes(ymin=-Inf,ymax=Inf,xmin=lower,xmax=upper,fill=PC),data=consonant_targets_HDI,alpha=0.4)+
  geom_rect(aes(ymin=-Inf,ymax=Inf,xmin=lower,xmax=upper,colour=PC),fill=NA,data=consonant_targets,alpha=0.4)+
  # geom_vline(aes(xintercept=end_time-dur),data=labels)+
  # geom_vline(aes(xintercept=end_time),data=labels)+
  # geom_text(aes(x=end_time - (dur / 2),y=Inf,label=label),data=labels_relative_consonant %>% filter(window_id %in% sample_consonant),vjust=2,family="Gentium Plus")+
  facet_grid(sensor~consonant,scales = "free_y") %>% suppressMessages()

ema_PCs_relative_vowel %>% filter(time>vowel_stability_onset-average_vowel_time,
                                  time<vowel_stability_offset-average_vowel_time,
                                  window_id %in% sample_vowel) %>% gather("PC","value",PC1,PC2) %>% ggplot()+
  # geom_line(aes(x=time,y=amplitude,group=window_id))+
  geom_density_ridges(aes(y=PC,x=value,colour=PC))+
  geom_rect(aes(ymin=-Inf,ymax=Inf,xmin=lower,xmax=upper,fill=PC),data=vowel_targets_HDI,alpha=0.4)+
  geom_rect(aes(ymin=-Inf,ymax=Inf,xmin=lower,xmax=upper,colour=PC),fill=NA,data=vowel_targets,alpha=0.4)+
  # geom_vline(aes(xintercept=end_time-dur),data=labels)+
  # geom_vline(aes(xintercept=end_time),data=labels)+
  # geom_text(aes(x=end_time - (dur / 2),y=Inf,label=label),data=labels_relative_consonant %>% filter(window_id %in% sample_consonant),vjust=2,family="Gentium Plus")+
  facet_grid(sensor~vowel,scales = "free_y")


# target comparison transition --------------------------------------------

ggplot()+
  geom_rect(aes(xmin=0,xmax=2,ymin=lower,ymax=upper,fill=PC),alpha=0.3,data=vowel_targets)+
  geom_rect(aes(xmin=1,xmax=3,ymin=lower,ymax=upper,fill=PC),alpha=0.3,data=consonant_targets)+
  geom_text(aes(x=-Inf,y=Inf,label=n,alpha=use),windows,vjust="inward",hjust="inward")+
  scale_alpha_manual(values = c(`TRUE`=1,`FALSE`=0.4))+
  facet_grid(consonant~vowel+sensor)
