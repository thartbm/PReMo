library(shiny)
library(optimx)
library(graphics)

# Define UI: data/schedule and model parameters
ui <- fluidPage(
    
    theme = bslib::bs_theme(bootswatch = "materia"),
    
    # https://bootswatch.com/
    # Lux
    # Materia
    
    
    # Application title
    titlePanel("PReMo"),
    
    tags$a(href="https://doi.org/10.1101/2021.12.21.473747", 
           "Implementation of the model in Tsay et al. (2021; bioRxiv)."),
    p("Currently only implements fully implicit adaptation."),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("upload", "Upload a schedule and data file", accept = c(".csv", ".tsv")),
            
            sliderInput("K_slider", "K", min = 0.001, max = 0.999, step=0.001, value = c(0.001, 0.999)),
            sliderInput("eta_p_slider", "eta_p", min = 0.001, max = 0.999, step=0.001, value = c(0.001, 0.999)),
            sliderInput("eta_v_slider", "eta_v", min = 0.001, max = 0.999, step=0.001, value = c(0.001, 0.999)),
            sliderInput("var_u_slider", "var_u", min = 0.100, max = 20.00, step=0.100, value = c(0.100, 20.00)),
            sliderInput("var_p_slider", "var_p", min = 0.100, max = 20.00, step=0.100, value = c(0.100, 20.00)),
            sliderInput("var_v_slider", "var_v", min = 0.100, max = 5.000, step=0.100, value = c(0.100, 5.000)),
            sliderInput("beta_vast_slider", "beta_vsat", min = 0.100, max = 10.00, step=0.100, value = c(0.100, 10.00)),
            sliderInput("beta_psat_slider", "beta_psat", min = 0.100, max = 15.00, step=0.100, value = c(0.100, 15.00)),
            
            actionButton("fit", "Fit model to data", class = "btn-lg btn-success", icon=icon("gears"))
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            p("File (.csv or .tsv), where rows=trials, 2 columns define the experiment: `rotation` and `trialtype` (1: normal, 2: error-clamped). Also one of two columns: `reachdev` (reach deviations) or `proprecal` (proprioceptive estimates). Measures relative to target in the same unit."),
            p("Parameters: free (fit in the range), or fixed (when minimum == maximum)."),
            plotOutput("fittedDataPlot"),
            tableOutput("fittedParameters")
        )
    )
)

# Define server logic 
server <- function(input, output) {
    
    source("R/PReMo.R")
    
    # uploaded data: schedule + actual data
    # dataframe of parameters
    
    parlimits <- data.frame('parameter'=c('K',  'eta_p', 'eta_v', 'var_u', 'var_p', 'var_v', 'beta_vsat', 'beta_psat'),
                            'lo'       =c(0.001, 0.001,  0.001,    0.001,    0.001,   0.001,   0.001,       0.001),
                            'hi'       =c(0.999, 0.999,  0.999,   20,       20,       5,      10,          15))
    
    model <- reactiveValues(dataschedule=NULL, parlimits=parlimits, parfitted=NULL)
    
    observeEvent(input$upload, {
        ext <- tools::file_ext(input$upload$name)
        model$dataschedule <- switch(ext,
                                     csv = read.csv(input$upload$datapath, sep = ","),
                                     tsv = read.csv(input$upload$datapath, sep = "\t"),
                                     validate("Invalid file; Please upload a .csv or .tsv file")
        )
    })
    
    
    observeEvent(input$fit, {
        validate(
            need(!is.null(model$dataschedule), "Please load some data to select.")
        )
        
        # this takes some time:
        fit <- PReMo_fit(dataschedule = model$dataschedule,
                         parlimits = model$parlimits)
        model$parfitted <- data.frame('parameter' = names(fit),
                                      'value' = as.numeric(fit))
        
        
        print(model$parfitted)
    })
    
    observeEvent(input$K_slider, {
        idx <- which(model$parlimits$parameter == 'K')
        model$parlimits$lo[idx] <- input$K_slider[1]
        model$parlimits$hi[idx] <- input$K_slider[2]
    })
    
    observeEvent(input$eta_p_slider, {
        idx <- which(model$parlimits$parameter == 'eta_p')
        model$parlimits$lo[idx] <- input$eta_p_slider[1]
        model$parlimits$hi[idx] <- input$eta_p_slider[2]
    })
    
    observeEvent(input$eta_v_slider, {
        idx <- which(model$parlimits$parameter == 'eta_v')
        model$parlimits$lo[idx] <- input$eta_v_slider[1]
        model$parlimits$hi[idx] <- input$eta_v_slider[2]
    })
    
    observeEvent(input$var_u_slider, {
        idx <- which(model$parlimits$parameter == 'var_u')
        model$parlimits$lo[idx] <- input$var_u_slider[1]
        model$parlimits$hi[idx] <- input$var_u_slider[2]
    })
    
    observeEvent(input$var_p_slider, {
        idx <- which(model$parlimits$parameter == 'var_p')
        model$parlimits$lo[idx] <- input$var_p_slider[1]
        model$parlimits$hi[idx] <- input$var_p_slider[2]
    })
    
    observeEvent(input$var_v_slider, {
        idx <- which(model$parlimits$parameter == 'var_v')
        model$parlimits$lo[idx] <- input$var_v_slider[1]
        model$parlimits$hi[idx] <- input$var_v_slider[2]
    })
    
    observeEvent(input$beta_vsat_slider, {
        idx <- which(model$parlimits$parameter == 'beta_vsat')
        model$parlimits$lo[idx] <- input$beta_vsat_slider[1]
        model$parlimits$hi[idx] <- input$beta_vsat_slider[2]
    })
    
    observeEvent(input$beta_psat_slider, {
        idx <- which(model$parlimits$parameter == 'beta_psat')
        model$parlimits$lo[idx] <- input$beta_psat_slider[1]
        model$parlimits$hi[idx] <- input$beta_psat_slider[2]
    })
    
    
    output$fittedDataPlot <- renderPlot({
        
        req(model$parfitted)
        
        # convert fitted par data frame to vector:
        parfitted <- model$parfitted$value
        names(parfitted) <- model$parfitted$parameter
        
        # run model on schedule & fitted parameters
        modelfit <- PReMo_model(parfree=parfitted,
                                parset=c(),
                                schedule=model$dataschedule)
        
        # determine range of Y values in graph:
        rangevalues <- c(range(dataschedule$rotation),
                         range(modelfit$reachdev),
                         range(modelfit$proprecal))
        if ('reachdev' %in% names(model$dataschedule))  { rangevalues <- c(rangevalues, range(model$dataschedule$reachdev))}
        if ('proprecal' %in% names(model$dataschedule)) { rangevalues <- c(rangevalues, range(model$dataschedule$proprecal))}
        
        # make plot:
        plot(-1000,-1000,
             main='',xlab='trial',ylab='deviation',
             xlim=c(0,(dim(dataschedule)[1]+1)),ylim=range(rangevalues),
             bty='n'
             )
        # add schedule:
        lines(dataschedule$rotation)
        
        # add data:
        if ('reachdev' %in% names(model$dataschedule))  { lines(model$dataschedule$reachdev,  col='#0fd2e2ff') }
        if ('proprecal' %in% names(model$dataschedule)) { lines(model$dataschedule$proprecal, col='#ff8200ff') }
        
        # add model processes:
        lines(modelfit$reachdev,  col='#005de4ff')
        lines(modelfit$proprecal, col='#e51636ff')
        
    })
    
    output$fittedParameters <- renderTable({
        model$parfitted
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
