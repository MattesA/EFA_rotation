

library(shiny)
library(tidyverse)
library(psych)
source("EFA_functions.R")


data <- bfi %>% 
  mutate_at(vars(A1, C4, C5, E1, E2, O2, O5), list(~ 7 - .)) %>% 
  select(starts_with(c("C", "E", "N")), -education) %>% 
  filter(complete.cases(.))
m.cor <- cor(data)


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Exploratorische Faktorenanalyse"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(width = 2,
      
      selectInput("method", h2("Methode"),
                  choices = list("PFA" = "pfa", "PCA" = "pca")),
      
      
      selectInput("rotationtype", h2("Rotation"),
                  choices = list("Varimax" = "varimax", "Equamax" = "equamax")),
      
      sliderInput("comp.pct", h2("Einfachstruktur"),
                  min = 0, max = 100, value = 0)
      
    ),
    
    
    # Main panel for displaying outputs ----
    mainPanel(

      fluidRow(
        
        column(3),

        column(3,
               h4("Unrotierte Ladungen"),
               plotOutput("load.unrotated")),
        
        column(3,
               h4("Annäherung an Einfachstruktur"),
               plotOutput("load.random")),
        
        column(3,
               h4("Rotierte Ladungen"),
               plotOutput("load.rotated"))

      ),

      fluidRow(

        column(3,
               h4("Beobachtete Korrelationen"),
               plotOutput("cor.observed")),
        
        column(3,
               h4("Modellbasierte Korrelationen"),
               plotOutput("cor.unrotated")),
        
        column(3,
               h4("Modellbasierte Korrelationen"),
               plotOutput("cor.random")),
        
        column(3,
               h4("Modellbasierte Korrelationen"),
               plotOutput("cor.rotated"))

      )

    )
  )
)




# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  output$load.unrotated <- renderPlot({
    out <- run.efa.unrotated(method = input$method, data = data, m.cor = m.cor)
    plot(out$load.unrotated)
  })
  
  output$load.random <- renderPlot({
    out <- run.simulations(method = input$method, rotation = input$rotationtype, desired.comp = input$comp.pct, data = data, m.cor = m.cor)
    out$load.random
  })
  
  output$load.rotated <- renderPlot({
    out <- run.efa.rotated(method = input$method, rotation = input$rotationtype, data = data, m.cor = m.cor)
    out$load.rotated
  })
  
  
  output$cor.observed <- renderPlot({
    cplot(m.cor)
  })
  
  
  output$cor.unrotated <- renderPlot({
    out <- run.efa.unrotated(method = input$method, data = data, m.cor = m.cor)
    out$cor.unrotated
  })
  
  output$cor.random <- renderPlot({
    out <- run.simulations(method = input$method, rotation = input$rotationtype, desired.comp = input$comp.pct, data = data, m.cor = m.cor)
    out$cor.random
  })
  
  output$cor.rotated <- renderPlot({
    out <- run.efa.rotated(method = input$method, rotation = input$rotationtype, data = data, m.cor = m.cor)
    out$cor.rotated
  })
  
}

