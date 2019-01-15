library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tuneR)
library(seewave)
library(stringr)
options(scipen = 9999)
cat(getwd())

spectro_data <- function(filename,window_onset=0,window_length=0.025,ovlp=90,rate="unknown",trial="unknown",ceiling=6000,dyn_range=-60) {
  # temp_root = "/tmp/joerod/euc/"
  # wav_file = paste0(temp_root,Sys.getpid(),".wav")
  # out_file = paste0(temp_root,Sys.getpid(),".out")
  # source_file = paste0("source_data/trial_recordings_speech_channel_only/",filename)
  # sox_call =  paste0("sox ", source_file," ", wav_file, " trim ", (window_onset)/1000," ", ((window_offset) - (window_onset))/1000)
  # write(sox_call,file = out_file,append = T)
  # system(sox_call,intern = F,ignore.stdout = T,ignore.stderr = T,wait = T)
  wave = readWave(filename)
  wl = wave@samp.rate * window_length
  spectrogram <- spectro(wave, plot = FALSE, ovlp=ovlp,wl=wl)
  frequency <- rep(spectrogram$freq*1000, times = ncol(spectrogram$amp))
  time <- window_onset+rep(spectrogram$time*1000, each = nrow(spectrogram$amp))
  amplitude <- as.vector(spectrogram$amp)
  df <- data.frame(time, frequency, amplitude) %>%
    filter(frequency<=ceiling,amplitude>=dyn_range) %>%
    mutate(experimental_trial_id = trial)
  return(df)
}

runMFCC_notrim = function(filename,window_period_ms=10,window_length_ms=25,samplefreq=16000){
  tempid = paste(sample(c(0:9, letters, LETTERS),8, replace=TRUE),collapse="")
  temp_root = "/data/clusterfs/pol/joerod/overlap_validate/temp/"
  wav_file = filename
  config_file = paste0(temp_root,"config/",tempid,'.config')
  mfc_file = paste0(temp_root,Sys.getpid(),".mfc")
  out_file = paste0(temp_root,"/euclid-worker",Sys.getpid(),".out")
  config = paste0('SOURCERATE = ',10000000/samplefreq,'\nSOURCEFORMAT = WAV\nSOURCEKIND = WAVEFORM\nTARGETRATE = ',
                  window_period_ms * 10000,
                  '\nTARGETFORMAT = HTK\nTARGETKIND = MFCC_E_D_A_Z\nZMEANSOURCE = T\nWINDOWSIZE = ',
                  window_length_ms * 10000,
                  '\nUSEHAMMING = T\nREEMCOEF = 0.97\nNUMCHANS = 26\nNUMCEPS = 12\nCEPLIFTER = 22\nENORMALISE = T\nDELTAWINDOW = 2\nSAVEWITHCRC = F\nBYTEORDER = VAX')
  write(config,file = config_file)
  
  hcopy_call = paste0("/usr/local/apps/htk/HCopy -C ",config_file," ",
                      wav_file, " ", mfc_file)
  write(hcopy_call,file = out_file,append = T)
  system(hcopy_call,intern = F,ignore.stdout = T,ignore.stderr = T,wait = T)
  while (!file.exists(mfc_file)) {
    write("waiting for HCopy",file = out_file,append = T)
    Sys.sleep(0.5)
  }
  hlist_call = paste0("/usr/local/apps/htk/HList -r -i 39 ",
                      mfc_file)
  write(hlist_call,file = out_file,append = T)
  q <- data.table::fread(hlist_call)
  write("HList done",file = out_file,append = T)
  colnames(q) <- c("MFCC1","MFCC2","MFCC3","MFCC4","MFCC5","MFCC6","MFCC7","MFCC8","MFCC9","MFCC10","MFCC11","MFCC12","E","Del1","Del2","Del3","Del4","Del5","Del6","Del7","Del8","Del9","Del10","Del11","Del12","DelE","Acc1","Acc2","Acc3","Acc4","Acc5","Acc6","Acc7","Acc8","Acc9","Acc10","Acc11","Acc12","AccE")
  q$sample <- seq(0,nrow(q),1)
  q$time <- ((q$sample) * window_period_ms) + window_length_ms/2
  q$filename <- rep(filename,nrow(q))
  return(q)
}

calc_euc <- function(mfcc_object,window_start_time=-Inf,window_end_time=Inf,short_filter_length = 30,long_filter_length = 90,summarised=F,scores=F,interpolation_interval=0.1){
  temp_root = "/data/clusterfs/pol/joerod/overlap_validate/temp/"
  out_file = paste0(temp_root,"/euclid-worker",Sys.getpid(),".out")
  write(paste0(short_filter_length," - ",long_filter_length),file = out_file,append = T)  
  euc_object <- mfcc_object %>%
    transmute(spectral_distance =
                sqrt(
                  (lag(MFCC1)-MFCC1)^2+
                    (lag(MFCC2)-MFCC2)^2+
                    (lag(MFCC3)-MFCC3)^2+
                    (lag(MFCC4)-MFCC4)^2+
                    (lag(MFCC5)-MFCC5)^2+
                    (lag(MFCC6)-MFCC6)^2+
                    (lag(MFCC7)-MFCC7)^2+
                    (lag(MFCC8)-MFCC8)^2+
                    (lag(MFCC9)-MFCC9)^2+
                    (lag(MFCC10)-MFCC10)^2+
                    (lag(MFCC11)-MFCC11)^2+
                    (lag(MFCC12)-MFCC12)^2
                ),
              energy_weighting = 
                (E+lag(E)/2), # already logged...
              instability_euclidean_back = abs(spectral_distance * energy_weighting),
              short_filter = smoother::smth.gaussian(instability_euclidean_back,window=short_filter_length/10),
              long_filter = smoother::smth.gaussian(instability_euclidean_back,window=long_filter_length/10),
              time=time)
  
  raw_func = splinefun(euc_object$time,euc_object$instability_euclidean_back)
  short_func = splinefun(euc_object$time,euc_object$short_filter)
  long_func = splinefun(euc_object$time,euc_object$long_filter)
  timeseq = seq(from = min(euc_object$time),to=max(euc_object$time),by=interpolation_interval)
  
  h = data.frame(time=timeseq,
                 interpolated_raw = raw_func(timeseq),
                 interpolated_short = short_func(timeseq),
                 interpolated_long = long_func(timeseq)) %>%
    mutate(interpolated_difference = interpolated_short - interpolated_long,
           interpolated_positive_difference = ifelse(interpolated_difference>=0,interpolated_difference+interpolated_long,NA),
           difference_is_positive = ifelse(interpolated_difference>=0,1,NA)) %>%
    filter(time>=window_start_time,time<window_end_time)
  
  return(h)
}

waveform_data <- function(filename,window_onset=0,window_length=0.025,ovlp=90,rate="unknown",trial="unknown",ceiling=6000,dyn_range=-60) {
  # temp_root = "/tmp/joerod/euc/"
  # wav_file = paste0(temp_root,Sys.getpid(),".wav")
  # out_file = paste0(temp_root,Sys.getpid(),".out")
  # source_file = paste0("source_data/trial_recordings_speech_channel_only/",filename)
  # sox_call =  paste0("sox ", source_file," ", wav_file, " trim ", (window_onset)/1000," ", ((window_offset) - (window_onset))/1000)
  # write(sox_call,file = out_file,append = T)
  # system(sox_call,intern = F,ignore.stdout = T,ignore.stderr = T,wait = T)
  wave = readWave(filename)
  wl = wave@samp.rate * window_length
  amplitude <- as.vector(wave@left)
  time <- (seq(1,length(amplitude)) / wave@samp.rate)*1000
  df <- data.frame(time,amplitude) %>%
    mutate(experimental_trial_id = trial)
  return(df)
}

gather_ema_data <- function(filename="ema_data/mngu0_s1_0001.ema",onset_time=0){
  g <- jsonlite::fromJSON(system2("python3",args=c("estfile.py",filename,"--json"),stdout = T)) %>% mutate(file=filename,time=onset_time+(timestamp*1000))
  
  p <- bind_cols(g %>% transmute(T3=T3_px,T2=T2_px,T1=T1_px,JAW=jaw_px,UL=upperlip_px,LL=lowerlip_px,REF=ref_px,time,file) %>% gather(sensor,x,T3,T2,T1,JAW,JAW,UL,LL,REF),
                 g %>% transmute(T3=T3_py,T2=T2_py,T1=T1_py,JAW=jaw_py,UL=upperlip_py,LL=lowerlip_py,REF=ref_py) %>% gather(sensor,y,T3,T2,T1,JAW,JAW,UL,LL,REF),
                 g %>% transmute(T3=T3_pz,T2=T2_pz,T1=T1_pz,JAW=jaw_pz,UL=upperlip_pz,LL=lowerlip_pz,REF=ref_pz) %>% gather(sensor,z,T3,T2,T1,JAW,JAW,UL,LL,REF)) %>% select(-sensor1,-sensor2)
  q <- p %>% group_by(sensor) %>% mutate(yz_delta =
                                  sqrt(
                                    ((lag(y, order_by=time)-y)^2)+
                                    ((lag(z, order_by=time)-z)^2)
                                  ),
                                yz_deltadelta = lag(yz_delta,order_by=time)-yz_delta
                                )
  return(ungroup(q))
}

gather_label_data <- function(filename="label_data/mngu0_s1_0001.lab"){
  temp_root = "/data/clusterfs/pol/joerod/overlap_validate/temp/"
  out_file = paste0(temp_root,"/euclid-worker",Sys.getpid(),".out")
  
  monophthongs = c("A","u","i","Q","a","O","I","V","U")
  fricatives_nasals = c("m","n","N","f","v","T","s","z","S","Z","l")
  
  target_sounds <- expand.grid(monophthongs,fricatives_nasals)
  
  labels <- data.table::fread(filename,header = F,drop=c(1),sep = "\t") %>% transmute(
    end_time = as.numeric(str_extract_all(V2,"\\d+\\.\\d+"))*1000,
    dur = end_time - lag(end_time,default = 0),
    time = lag(end_time,default = 0) + dur / 2,
    label = V3
  )
  
  windows <- data.frame(window = paste(labels$label,lead(labels$label)))
  windows$interesting <- paste(labels$label,lead(labels$label)) %in% paste(target_sounds$Var1,target_sounds$Var2) 
  windows$interesting_window_start <- labels$time
  windows$interesting_window_end <- lead(labels$time)
  
  windows <- filter(windows,interesting)
  
  return(list(labels=labels,
              windows=windows,
              begin_speech=labels[1,1]))
  
}

process_file <- function(filename="mngu0_s1_0001"){
  out <- list()
  labels <- gather_label_data(paste0("label_data/",filename,".lab"))
  out$labels <- labels$labels %>% mutate(filename=filename)
  out$windows <- labels$windows %>% mutate(filename=filename)
  out$ema <- gather_ema_data(paste0("ema_data/",filename,".ema"),onset_time = 0) %>% mutate(filename=filename)
  
  out$spec <- spectro_data(paste0("audio_data/",filename,".wav")) %>% mutate(filename=filename)
  out$wf <- waveform_data(paste0("audio_data/",filename,".wav")) %>% mutate(filename=filename)
  out$mfcc <- runMFCC_notrim(paste0("/data/clusterfs/pol/joerod/overlap_validate/mfcc/overlap-validate/audio_data/",filename,".wav")) %>% mutate(filename = filename)
  out$euc <- calc_euc(out$mfcc) %>% mutate(filename =filename)
  
  out$model_data <- left_join(out$euc %>% filter(round(time,1)%in%round(out$ema$time,1)),
                          left_join(
                            out$ema %>% select(sensor,yz_delta,time)%>%mutate(time=round(time,1)) %>% ungroup()  %>% mutate(sensor=paste0(sensor,"_delta")) %>% spread(sensor,yz_delta,),
                            out$ema %>% select(sensor,yz_deltadelta,time)%>%mutate(time=round(time,1))%>% ungroup() %>% mutate(sensor=paste0(sensor,"_deltadelta")) %>% spread(sensor,yz_deltadelta)))  %>% mutate(filename =filename)
  return(out)
}

library(doParallel)
library(foreach)
library(multidplyr)

cl <- makeCluster(19)
cluster_library(cl,c("dplyr","tidyr","jsonlite","tuneR","seewave","stringr"))
cluster_eval(cl,options(scipen = 9999))
cluster_copy(cl,process_file)
cluster_copy(cl,spectro_data)
cluster_copy(cl,waveform_data)
cluster_copy(cl,gather_ema_data)
cluster_copy(cl,gather_label_data)
cluster_copy(cl,runMFCC_notrim)
cluster_copy(cl,calc_euc)
registerDoParallel(cl)

mngu0 <- plyr::llply(na.omit(str_match(list.files("audio_data/"),pattern = "(.*)\\.wav")[,2])[1:100],process_file,.parallel = TRUE,.paropts = list(.errorhandling="remove"))
model_data <- plyr::ldply(mngu0,function(x){x$model_data})
windows <- plyr::ldply(mngu0,function(x){x$windows})
save(model_data,windows,file="model_data_untrimmed.RData")
save(mngu0,file="mngu0_untrimmed.RData")