library(tidyverse)
load("annotate/static_data.Rdata")
load("annotate/functions.Rdata")


library(nloptr)

max_density <- function(x){
  density(x)$x[which.max(density(x)$y)] 
}


mround <- function(x,base){ 
  base*round(x/base) 
} 

windows_to_include_phase2 <- data.frame(window = c("a m", "a S", "a v", "i m", "i s", "i v", "I v", "I S", "Q n", "V m", "V s"),
                                        vowel=c("a","a","a", "i", "i", "i", "I", "I", "Q", "V", "V"),
                                        consonant=c("m", "S", "v", "m", "s", "v", "v", "S", "n", "m", "s"),use = T)

windows_to_include_phase3 <- data.frame(window = c("a m", "a S", "a v", "i m", "i s", "i v", "I v", "I S", "Q n", "V m", "V s"),
                                        vowel=c("a","a","a", "i", "i", "i", "I", "I", "Q", "V", "V"),
                                        consonant=c("m", "S", "v", "m", "s", "v", "v", "S", "n", "m", "s"),use = T)


load("intermediate_data/annotation_data_20180821.Rdata")


annotated <- target_windows_annotation %>% ungroup()  %>% mutate(vowel = str_split(as.character(window)," ",simplify = T)[,1],
                                                                 consonant = str_split(as.character(window)," ",simplify = T)[,2]) %>%
  filter(!is.na(vowel_stability_onset)&!is.na(consonant_stability_onset))

load("intermediate_data/targets.Rdata")


# function definitions ----------------------------------------------------

interpolate_targets <- function(target_vowel="V", target_consonant="v", vowel_onset=0,vowel_offset=120,consonant_onset=110,consonant_offset=240,n=50){
  
  interpolate_targets_inner <- function(x){
    
    vowel_target <- vowel_targets %>% filter(vowel==target_vowel, sensor==x$sensor,PC==x$PC) %>% mutate(HDI_index = paste0("v",HDI_index))
    consonant_target <- consonant_targets %>% filter(consonant==target_consonant, sensor==x$sensor,PC==x$PC) %>% mutate(HDI_index = paste0("c",HDI_index))
    
    smooth_series <- function(x,n){
      averaged = x %>% group_by(time) %>% summarise(upper=mean(upper),
                                                    lower=mean(lower))
      # smooth = 
      # smooth = pracma::movavg(averaged$target,n = 50,type = "t")
      # smooth::es
      # smooth = loess.smooth(x=x$time,y=smooth,span = 2/3)
      return(data.frame(time=averaged$time,
                        smoothed_upper = pracma::movavg(averaged$upper,n = n,type = "t"),
                        smoothed_lower = pracma::movavg(averaged$lower,n = n,type = "t")
      ))
    }
    
    interpolate_targets_inner2 <- function(x) {
      vowel_sub_target <- vowel_target %>% filter(HDI_index == x$Var1)
      consonant_sub_target <- consonant_target %>% filter(HDI_index == x$Var2)
      
      combined_targets <- bind_rows(full_join(data.frame(time=seq(mround(vowel_onset,1),vowel_offset,1),vowel=target_vowel,consonant=""),vowel_sub_target,by="vowel"),
                                    full_join(data.frame(time=seq(mround(consonant_onset,1),consonant_offset,1),vowel="",consonant=target_consonant),consonant_sub_target,by="consonant"))
      
      left_join(combined_targets %>% group_by(PC,sensor) %>% do(smooth_series(.,n=n)),
                combined_targets, by = c("PC", "sensor", "time"))
      
    }
    
    expand.grid(names(table(vowel_target$HDI_index)),names(table(consonant_target$HDI_index))) %>% mutate(pair=row_number()) %>%
      group_by(pair) %>% do(interpolate_targets_inner2(.))
    
  }
  
  
  expand.grid(c("JAW","LL","T1","T2","T3","UL"),c("PC1","PC2")) %>% transmute(sensor=Var1,PC=Var2) %>% rowwise() %>% do(interpolate_targets_inner(.))
  
  
}

make_starting_points <- function(observed,meta){
  
  starting_vowel_offset = observed[1,]$consonant_stability_onset
  starting_consonant_onset = observed[1,]$vowel_stability_offset
  consonant_stability_offset = observed[1,]$consonant_stability_offset
  vowel_stability_onset = observed[1,]$vowel_stability_onset
  
  data.frame(
    starting_vowel_offset = rnorm(mean=starting_vowel_offset,sd=meta$starting_points_sd,meta$starting_points_n*10),
    starting_consonant_onset = rnorm(mean=starting_consonant_onset,sd=meta$starting_points_sd,meta$starting_points_n*10)) %>%
    filter(starting_consonant_onset < starting_vowel_offset,
           starting_vowel_offset <= consonant_stability_offset,
           starting_consonant_onset >= vowel_stability_onset) %>% sample_n(meta$starting_points_n) %>% mutate(sample = seq(1,meta$starting_points_n))
}

bootstrapped_settings_test <- function(observed,meta){
  if(meta$optimiser == "'bobyqa'"){
    
    optimise_one_bobyqa <- function(vowel_offset,consonant_onset,observed,meta){
      tt <- observed[1,]
      optimisation <- bobyqa(x0=c(vowel_offset,consonant_onset),
                             upper=c(tt$consonant_stability_offset,Inf),
                             lower=c(-Inf,tt$vowel_stability_onset),
                             fn = function(par){test_settings(observed,par,meta)$score},
                             control = list(ftol_abs=0.1))
      return(data.frame(vowel_offset = optimisation$par[1],consonant_onset = optimisation$par[2],
                        score=optimisation$value,
                        iter=optimisation$iter,
                        convergence=optimisation$convergence,
                        message=optimisation$message))
    }
    
    starting_points <- make_starting_points(observed,meta=meta)
    
    if(meta$parallel){
      cluster_copy(cl,obj=observed)
      cluster_copy(cl,obj=meta)
      cluster_copy(cl,obj = optimise_one_bobyqa)
      results <- starting_points %>% partition(cluster = cl) %>% group_by(sample,starting_vowel_offset,starting_consonant_onset) %>% do(optimise_one_bobyqa(.$starting_vowel_offset,.$starting_consonant_onset,observed,meta=meta)) %>% collect()
    } else {
      results <- starting_points %>% group_by(sample,starting_vowel_offset,starting_consonant_onset) %>% do(optimise_one_bobyqa(.$starting_vowel_offset,.$starting_consonant_onset,observed,meta=meta))
    }
    
  } else if (meta$optimiser == "'hydroPSO'") {
    
    starting_points <- make_starting_points(observed,meta=meta) %>%
      select(starting_vowel_offset,starting_consonant_onset)
    
    optimisation <- hydroPSO(par=starting_points,
                             upper=c(tt$consonant_stability_offset,tt$consonant_stability_offset),
                             lower=c(tt$vowel_stability_onset,tt$vowel_stability_onset),
                             fn = function(par){test_settings_C(observed,par,meta)$score},
                             control = list(parallel="parallel",
                                            par.nnodes=3,
                                            maxit=30,
                                            npart=12)
    )
  
    
    return(data.frame(vowel_offset = optimisation$par[1],consonant_onset = optimisation$par[2],
                      score=optimisation$value,
                      iter=optimisation$counts[2],
                      convergence=optimisation$convergence,
                      message=optimisation$message))  
  return(optimisation)
  }
}



test_settings <- function(observed,par,meta){
  
  tt <- observed[1,]
  vowel_offset = par[1]
  consonant_onset = par[2]
  
  if(mround(vowel_offset,1) < mround(consonant_onset,1)){ # i.e. no overlap
    return(list(score=4,distance_score=4,reason="no overlap"))
  }
  
  if (mround(vowel_offset,1) >= mround(tt$consonant_stability_offset,1)){ # i.e. 100% overlap
    return(list(score=4,distance_score=4,reason="100% overlap in consonant"))
  }
  
  if (mround(consonant_onset,1) <= mround(tt$vowel_stability_onset,1)){ # i.e. 100% overlap
    return(list(score=4,distance_score=4,reason="100% overlap in vowel"))
  }
  
  targets <- interpolate_targets(target_vowel = tt$vowel,
                                 target_consonant = tt$consonant,
                                 vowel_onset = tt$vowel_stability_onset,
                                 vowel_offset = par[1],
                                 consonant_onset = par[2],
                                 consonant_offset = tt$consonant_stability_offset,
                                 n=meta$interpolation_n) %>% mutate(time=round(time)) %>% filter(time%%5==0)
  
  outcome <- full_join(observed %>% ungroup() %>% select(-vowel,-consonant),targets,by=c("time","sensor","PC")) %>%
    mutate(in_bounds = ifelse(value>=smoothed_lower&value<=smoothed_upper,T,F),
           distance_from_boundary = ifelse(value<smoothed_lower,abs(value-smoothed_lower),ifelse(value>smoothed_upper,abs(value-smoothed_upper),0))
    )
  
  if(meta$score_function %in% c("'multi_within'", "'multi_distance'")){
    outcome_for_score <- outcome %>% filter(time>=(tt$vowel_stability_onset+tt$vowel_stability_offset)/2,
                                            time<=(tt$consonant_stability_onset+tt$consonant_stability_offset)/2) %>%
      rowwise() %>% mutate(target_centre = mean(c(smoothed_upper,smoothed_lower)),
                           distance_sq = (value-target_centre)^2) %>% ungroup()
    
    outcome_for_score_summarised <- outcome_for_score %>% group_by(time) %>% summarise(in_bounds=!FALSE%in%in_bounds,
                                                                                       distance = sqrt(sum(distance_sq)))
  } else {
    
  }
  
  score = NA
  
  if(meta$score_function == "'multi_within'"){
    score = 1-nrow(outcome_for_score_summarised[which(outcome_for_score_summarised$in_bounds==T),]) / nrow(outcome_for_score_summarised)
  } else if( meta$score_function == "'multi_distance'") {
    score = mean(outcome_for_score_summarised$distance,na.rm = T)
  }
  
  plot = NA
  
  if(meta$make_plot){
    plot = ggplot(outcome %>% mutate(panel="fit"))+
      geom_ribbon(aes(x=time,ymax=upper,ymin=lower,colour=PC,group=paste(PC,sensor,vowel,consonant,pair)),fill=NA,alpha=0.4)+
      geom_ribbon(aes(x=time,ymax=smoothed_upper,ymin=smoothed_lower,fill=PC,fill=paste(PC,pair)),alpha=0.4)+
      geom_line(aes(x=time,y=value,colour=PC),size=2)+
      geom_vline(xintercept=(tt$vowel_stability_onset+tt$vowel_stability_offset)/2,alpha=0.5)+
      geom_text(aes(x=(tt$vowel_stability_onset+tt$vowel_stability_offset)/2,y = Inf,label=tt$vowel),vjust="inward",hjust="center")+
      geom_vline(xintercept=(tt$consonant_stability_onset+tt$consonant_stability_offset)/2,alpha=0.5)+
      geom_text(aes(x=(tt$consonant_stability_onset+tt$consonant_stability_offset)/2,y = Inf,label=tt$consonant),vjust="inward",hjust="center")+
      facet_grid(panel~sensor,scales = "free_y")
  } 
  
  return(list(outcome=outcome,
              outcome_for_score=outcome_for_score,
              outcome_for_score_summarised=outcome_for_score_summarised,
              targets=targets,
              score = score,
              plot = plot,
              overlap_duration = par[1] - par[2]
  ))
}


# run a single optimisation


# grid task  --------------------------------------------------------------


# sample <- annotated %>% group_by(index,window,window_id,vowel,consonant,consonant_stability_onset,consonant_stability_offset,vowel_stability_onset,vowel_stability_offset)
# observed <-sample %>% do(prepare_ema(.)) %>% gather("PC","value",PC1,PC2)  %>% mutate(time=round(time))
# 
# 
# 
# observed_n <- observed %>% filter(PC=="PC1",sensor=="JAW") %>% summarise(n=n(),shard=base::sample(seq(1,4),1))
#  
# already_done <- names(table(grid_result$window_id))
#  
# observed_n[which(!observed_n$window_id%in%already_done),"shard"] <- sample(seq(4,10),218,T)
# 
# observed <- left_join(observed,observed_n)
# 
# 
# save(sample,observed,file = "intermediate_data/testing_sample.RData")

load(file = "intermediate_data/testing_sample.RData")


library(multidplyr)
cl <- multidplyr::create_cluster(70)
cluster_library(cl,packages = c("nloptr","tidyverse","hydroPSO"))
cluster_copy(cl,obj = test_settings)
cluster_copy(cl,obj = interpolate_targets)
cluster_copy(cl,obj = vowel_targets)
cluster_copy(cl,obj = consonant_targets)
cluster_copy(cl,obj = mround)


# parallelisation happens inside, that's why not obvious here...

grid_result <- data.frame()

# metaparameters <- list('optimiser' = "'bobyqa'",
#                        'interpolation_n' = 20,
#                        'starting_points_n' = 40,
#                        'starting_points_sd' = 25,
#                        'score_function' = "'multi_within'",
#                        'make_plot'=F,
#                        'parallel'=T)
# outcome <- observed %>% do(bootstrapped_settings_test(.,meta=metaparameters))
# grid_result <- bind_rows(grid_result,
#                          outcome %>% mutate_(.dots=metaparameters))
# save(grid_result,file="grid_results/grid_result_multi_within.RData")
# 
# metaparameters <- list('optimiser' = "'bobyqa'",
#                        'interpolation_n' = 20,
#                        'starting_points_n' = 40,
#                        'starting_points_sd' = 25,
#                        'score_function' = "'multi_within'",
#                        'make_plot'=F,
#                        'parallel'=T)
# outcome <- observed %>% do(bootstrapped_settings_test(.,meta=metaparameters))
# grid_result <- bind_rows(grid_result,
#                          outcome %>% mutate_(.dots=metaparameters))
# save(grid_result,file="grid_results/grid_result_multi_within.RData")
# 
# 
# metaparameters <- list('optimiser' = "'bobyqa'",
#                        'interpolation_n' = 25,
#                        'starting_points_n' = 40,
#                        'starting_points_sd' = 25,
#                        'score_function' = "'multi_within'",
#                        'make_plot'=F,
#                        'parallel'=T)
# outcome <- observed %>% do(bootstrapped_settings_test(.,meta=metaparameters))
# grid_result <- bind_rows(grid_result,
#                          outcome %>% mutate_(.dots=metaparameters))
# save(grid_result,file="grid_results/grid_result_multi_within.RData")


# metaparameters <- list('optimiser' = "'bobyqa'",
#                        'interpolation_n' = 30,
#                        'starting_points_n' = 40,
#                        'starting_points_sd' = 25,
#                        'score_function' = "'multi_within'",
#                        'make_plot'=F,
#                        'parallel'=T)
# outcome <- observed %>% do(bootstrapped_settings_test(.,meta=metaparameters))
# grid_result <- bind_rows(grid_result,
#                          outcome %>% mutate_(.dots=metaparameters))
# save(grid_result,file="grid_results/grid_result_multi_within2.RData")


# metaparameters <- list('optimiser' = "'bobyqa'",
#                        'interpolation_n' = 35,
#                        'starting_points_n' = 40,
#                        'starting_points_sd' = 25,
#                        'score_function' = "'multi_within'",
#                        'make_plot'=F,
#                        'parallel'=T)
# outcome <- observed %>% do(bootstrapped_settings_test(.,meta=metaparameters))
# grid_result <- bind_rows(grid_result,
#                          outcome %>% mutate_(.dots=metaparameters))
# save(grid_result,file="grid_results/grid_result_multi_within3.RData")


metaparameters <- list('optimiser' = "'bobyqa'",
                       'interpolation_n' = 40,
                       'starting_points_n' = 200,
                       'starting_points_sd' = 25,
                       'score_function' = "'multi_within'",
                       'make_plot'=F,
                       'parallel'=T)
# outcome <- observed %>% filter(shard==1,n>metaparameters$interpolation_n+5) %>% do(bootstrapped_settings_test(.,meta=metaparameters))
# grid_result <- bind_rows(grid_result,
#                          outcome %>% mutate_(.dots=metaparameters))
# save(outcome,file="grid_results/grid_result_multi_within_shard1.RData")
# save(grid_result,file="grid_results/grid_result_multi_within_definitive.RData")
# 
# outcome <- observed %>% filter(shard==2,n>metaparameters$interpolation_n+5) %>% do(bootstrapped_settings_test(.,meta=metaparameters))
# grid_result <- bind_rows(grid_result,
#                          outcome %>% mutate_(.dots=metaparameters))
# save(outcome,file="grid_results/grid_result_multi_within_shard2.RData")
# save(grid_result,file="grid_results/grid_result_multi_within_definitive.RData")
# 
# outcome <- observed %>% filter(shard==3,n>metaparameters$interpolation_n+5) %>% do(bootstrapped_settings_test(.,meta=metaparameters))
# grid_result <- bind_rows(grid_result,
#                          outcome %>% mutate_(.dots=metaparameters))
# save(outcome,file="grid_results/grid_result_multi_within_shard3.RData")
# save(grid_result,file="grid_results/grid_result_multi_within_definitive.RData")

outcome <- observed %>% filter(shard==4,n>metaparameters$interpolation_n+5) %>% do(bootstrapped_settings_test(.,meta=metaparameters))
grid_result <- bind_rows(grid_result,
                         outcome %>% mutate_(.dots=metaparameters))
save(outcome,file="grid_results/grid_result_multi_within_shard4.RData")
save(grid_result,file="grid_results/grid_result_multi_within_definitive.RData")

outcome <- observed %>% filter(shard==5,n>metaparameters$interpolation_n+5) %>% do(bootstrapped_settings_test(.,meta=metaparameters))
grid_result <- bind_rows(grid_result,
                         outcome %>% mutate_(.dots=metaparameters))
save(outcome,file="grid_results/grid_result_multi_within_shard5.RData")
save(grid_result,file="grid_results/grid_result_multi_within_definitive.RData")

outcome <- observed %>% filter(shard==6,n>metaparameters$interpolation_n+5) %>% do(bootstrapped_settings_test(.,meta=metaparameters))
grid_result <- bind_rows(grid_result,
                         outcome %>% mutate_(.dots=metaparameters))
save(outcome,file="grid_results/grid_result_multi_within_shard6.RData")
save(grid_result,file="grid_results/grid_result_multi_within_definitive.RData")

outcome <- observed %>% filter(shard==7,n>metaparameters$interpolation_n+5) %>% do(bootstrapped_settings_test(.,meta=metaparameters))
grid_result <- bind_rows(grid_result,
                         outcome %>% mutate_(.dots=metaparameters))
save(outcome,file="grid_results/grid_result_multi_within_shard7.RData")
save(grid_result,file="grid_results/grid_result_multi_within_definitive.RData")

outcome <- observed %>% filter(shard==8,n>metaparameters$interpolation_n+5) %>% do(bootstrapped_settings_test(.,meta=metaparameters))
grid_result <- bind_rows(grid_result,
                         outcome %>% mutate_(.dots=metaparameters))
save(outcome,file="grid_results/grid_result_multi_within_shard8.RData")
save(grid_result,file="grid_results/grid_result_multi_within_definitive.RData")
