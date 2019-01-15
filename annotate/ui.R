#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(tags$script(src = "message-handler.js")),
  
  # Application title
  
  # Sidebar with a slider input for number of bins 
  
  fluidPage(
    fluidRow(
      column(7,
             
             plotOutput('mainPlot',
                        brush = brushOpts(
                          id = "plot_brush",
                          direction = "x",
                          resetOnNew = T
                        ),
                        height = "50vh",
                        width = "56vw"),
             actionButton("storeVowel", "Store vowel stable period",style="background-color: rgb(72,188,255);"),
             actionButton("storeConsonant", "Store consonant stable period",style="background-color: rgb(144,238,144);"),
             actionButton("nextButton", "Proceed to next unannotated"),
             actionButton("saveButton", "Save data")
             
      ),
      column(5,
             fluidRow(
               plotOutput('keyPlot',
                          height = "50vh",
                          width = "38vw"
               )
             )
             # fluidRow(
             # verbatimTextOutput('temp')
             # )
      )),
    fluidRow(
      DTOutput('tbl')
    ),
    fluidRow(
      column(6, verbatimTextOutput('x5'))
    )
  )
)