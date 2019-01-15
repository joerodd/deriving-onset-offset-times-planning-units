#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

shinyServer(function(input, output, session) {
  
  target_data <- reactive({
    input$storeConsonant
    input$storeVowel
    
    target_windows_annotation %>% select(window_id,window,vowel_stability_onset,vowel_stability_offset,vowel_stability_sensor,
                                         consonant_stability_onset,consonant_stability_offset,consonant_stability_sensor)
    
  })
  
  output$tbl <- renderDataTable(isolate({
    target_data()
      }),
  selection = list(mode = 'single'),
  options=list(searching=FALSE)
  )
  
  observeEvent(input$storeConsonant,
               {
                 if(!is.null(input$plot_brush)){
                   s <- target_windows_annotation %>% ungroup()
                   selected_row <<- input$tbl_rows_selected[1]
                   consonant_stability_onset <<- input$plot_brush$xmin
                   consonant_stability_offset <<- input$plot_brush$xmax
                   consonant_stability_sensor <<- input$plot_brush$panelvar1
                   
                   
                   s[selected_row,"consonant_stability_onset"] <-  consonant_stability_onset
                   s[selected_row,"consonant_stability_offset"] <- consonant_stability_offset
                   s[selected_row,"consonant_stability_sensor"] <- consonant_stability_sensor
                   
                   target_windows_annotation <<- s
                   
                   dataTableProxy('tbl') %>% replaceData(target_data()) %>% selectRows(selected_row) %>% selectPage((selected_row %/% 10) + 1)
                 }
               })
  
  observeEvent(input$storeVowel,
               {
                 if(!is.null(input$plot_brush)){
                   s <- target_windows_annotation %>% ungroup()
                   selected_row <<- input$tbl_rows_selected[1]
                   
                   vowel_stability_onset <<- input$plot_brush$xmin
                   vowel_stability_offset <<- input$plot_brush$xmax
                   vowel_stability_sensor <<- input$plot_brush$panelvar1
                   
                   s[selected_row,"vowel_stability_onset"] <-  vowel_stability_onset
                   s[selected_row,"vowel_stability_offset"] <- vowel_stability_offset
                   s[selected_row,"vowel_stability_sensor"] <- vowel_stability_sensor
                   
                   target_windows_annotation <<- s
                   dataTableProxy('tbl') %>% replaceData(target_data()) %>% selectRows(selected_row) %>% selectPage((selected_row %/% 10) + 1)
                 }
               })
  
  observeEvent(input$saveButton, {
    save(target_windows_annotation,file = "data_in_progress.Rdata")
  })
  observeEvent(input$nextButton, {
    selected_row <- max(which(!is.na(target_windows_annotation$vowel_stability_onset))) + 1
    dataTableProxy('tbl') %>% replaceData(target_data()) %>% selectRows(selected_row) %>% selectPage((selected_row %/% 10) + 1) 
    # cat("\ninput$tbl:\n")
    # cat(str(input$tbl_length))
    cat("\npage number:\n")
  })
  
  output$keyPlot <- renderPlot({
    
    input$storeConsonant
    input$storeVowel
    
    selected_window_id <- unlist(target_windows_annotation[unlist(input$tbl_rows_selected),"window_id"])
    
    rr <- target_windows_annotation %>% filter(window_id==selected_window_id)
    
    ema_PCs_vowel <- rr %>% group_by(window,window_id,index) %>% do(prepare_ema(.)) %>% filter(sensor%in%c("LL","UL","JAW","T1","T2","T3"),time>=rr$vowel_stability_onset,time<=rr$vowel_stability_offset) %>% mutate(sample_number = as.numeric(row_number()))
    ema_PCs_consonant <- rr %>% group_by(window,window_id,index) %>% do(prepare_ema(.)) %>% filter(sensor%in%c("LL","UL","JAW","T1","T2","T3"),time>=rr$consonant_stability_onset,time<=rr$consonant_stability_offset) %>% mutate(sample_number = as.numeric(row_number()))
    ema_PCs <- rr %>% group_by(window,window_id,index) %>% do(prepare_ema(.)) %>% filter(sensor%in%c("LL","UL","JAW","T1","T2","T3"),time>=rr$interesting_window_start - 100,time<=rr$interesting_window_end + 100) %>% mutate(sample_number = as.numeric(row_number()))
    
    data_for_pcas_key %>% ggplot()+
      # geom_point(aes(x=y,y=z),data=ema,size=0.3,alpha=0.3)+
      geom_path(aes(x=y,y=z,group=sensor),colour="yellow",data=ema_PCs,size=1)+
      geom_path(aes(x=y,y=z,group=sensor),colour="royalblue1",data=ema_PCs_vowel,size=2,arrow = arrow())+
      geom_path(aes(x=y,y=z,group=sensor),colour="palegreen2",data=ema_PCs_consonant,size=2,arrow = arrow())+
      # geom_rug(aes(x=y,y=z),data=sample_n(ema,5000),size=0.3,alpha=0.1)+
      geom_segment(aes(x=y_origin,y=z_origin,xend=y_end,yend=z_end,colour=PC),arrow = arrow())+
      geom_text(aes(x=y_origin,y=z_origin,label=sensor),colour="black")
  })
  
  output$mainPlot <- renderPlot({
    
    input$storeConsonant
    input$storeVowel
    
    selected_window_id <- unlist(target_windows_annotation[unlist(input$tbl_rows_selected),"window_id"])
    
    rr <- target_windows_annotation %>% filter(window_id==selected_window_id) %>% ungroup()
    
    ema_PCs <- rr %>% group_by(window,window_id,index) %>% do(prepare_ema_wf(.)) %>% ungroup()
    labels <-  rr %>% group_by(window,window_id,index) %>% do(prepare_labels(.)) %>% ungroup()
    
    ema_PCs %>% gather("PC","value",PC1,PC2) %>% ggplot()+
      geom_line(aes(x=time,y=amplitude))+
      geom_rect(aes(xmin=vowel_stability_onset,ymin=-Inf,ymax=Inf,xmax=vowel_stability_offset),data=rr,fill="royalblue1",alpha=0.4)+
      geom_rect(aes(xmin=consonant_stability_onset,ymin=-Inf,ymax=Inf,xmax=consonant_stability_offset),data=rr,fill="palegreen2",alpha=0.4)+
      geom_rug(aes(y=value,colour=PC),data=ema_for_rugplots %>% filter(sensor%in%c("LL","UL","JAW","T1","T2","T3")))+
      geom_line(aes(x=time,y=value,colour=PC))+
      geom_vline(aes(xintercept=end_time-dur),data=labels)+
      geom_vline(aes(xintercept=end_time),data=labels)+
      geom_text(aes(x=end_time - (dur / 2),y=Inf,label=label),data=labels,vjust=2,family="Gentium Plus")+
      facet_grid(sensor~.,scales = "free_y")
    
  })
  
})