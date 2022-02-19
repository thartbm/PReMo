library(shiny)
library(optimx)
library(graphics)

# Define UI: data/schedule and model parameters
ui <- fluidPage(
    
    theme = bslib::bs_theme(bootswatch = "lumen"),
    
    # https://bootswatch.com/
    # Lux
    # Materia
    # Yeti
    # Lumen
    
    
    # Application title
    titlePanel("PReMo"),
    
    tags$a(href="https://doi.org/10.1101/2021.12.21.473747", 
           "Implementation of the model in Tsay et al. (2021; bioRxiv)."),
    p("Currently only implements fully implicit adaptation."),
    tags$a(href="https://github.com/thartbm/PReMo", "Source code on GitHub."),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            fluidRow(
                
                column(8,
                       sliderInput("K_slider", "K", min = 0.001, max = 0.999, step=0.001, value = 0.50),
                ),
                column(4,
                       checkboxInput("K_fixed", "fix K", value = FALSE),
                ),
                column(8,
                       sliderInput("eta_p_slider", "eta_p", min = 0.001, max = 0.999, step=0.001, 0.75),
                ),
                column(4,
                       checkboxInput("eta_p_fixed", "fix eta_p", value = FALSE),
                ),
                column(8,
                       sliderInput("eta_v_slider", "eta_v", min = 0.001, max = 0.999, step=0.001, 0.75),
                ),
                column(4,
                       checkboxInput("eta_v_fixed", "fix eta_v", value = FALSE),
                ),
                column(8,
                       sliderInput("var_u_slider", "var_u", min = 0.100, max = 20.00, step=0.100, 2.00),
                ),
                column(4,
                       checkboxInput("var_u_fixed", "fix var_u", value = FALSE),
                ),
                column(8,
                       sliderInput("var_p_slider", "var_p", min = 0.100, max = 20.00, step=0.100, 7.00),
                ),
                column(4,
                       checkboxInput("var_p_fixed", "fix var_p", value = FALSE),
                ),
                column(8,
                       sliderInput("var_v_slider", "var_v", min = 0.100, max = 5.000, step=0.100, 0.50),
                ),
                column(4,
                       checkboxInput("var_v_fixed", "fix var_v", value = FALSE),
                ),
                column(8,
                       sliderInput("beta_vsat_slider", "beta_vsat", min = 0.100, max = 10.00, step=0.100, value = 0.5),
                ),
                column(4,
                       checkboxInput("beta_vsat_fixed", "fix beta_vsat", value = FALSE),
                ),
                column(8,
                       sliderInput("beta_psat_slider", "beta_psat", min = 0.100, max = 15.00, step=0.100, value = 5.0),
                ),
                column(4,
                       checkboxInput("beta_psat_fixed", "fix beta_psat", value = FALSE),
                ),
                column(8,
                       sliderInput("rotation_slider", "rotation", min = -90.00, max = 90.00, step=1.00, value = -30.0),
                ),
            ),
            
            actionButton("reset", "Reset", class = "btn-danger", icon=icon("sync")),
            hr(),
            fileInput("upload", "Upload a schedule and data file", accept = c(".csv", ".tsv")),
            actionButton("unload", "Unload file", class = "btn-danger", icon=icon("trash")),
            p("File (.csv or .tsv), where rows=trials, 2 columns define the experiment: `rotation` and `trialtype` (1: normal, 2: error-clamped). At least one of two columns: `reachdev` (reach deviations) or `proprecal` (proprioceptive estimates). Measures relative to target in the same unit."),
            actionButton("fit", "Fit model to data", class = "btn-lg btn-success", icon=icon("gears")),
            p("Parameters: free (starting value from sliders), or fixed (when ticked)."),
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("fittedDataPlot"),
            tableOutput("fittedParameters"),
        )
    )
)

# Define server logic 
server <- function(input, output) {
    
    source("R/PReMo.R")
    
    # uploaded data: schedule + actual data
    # dataframe of parameters
    
    storage <- reactiveValues(dataschedule=NULL, parfitted=NULL)
    
    modelsettings <- reactive({
        parameter <- c('K', 'eta_p', 'eta_v', 'var_u', 'var_v', 'var_p', 'beta_vsat', 'beta_psat')
        value     <- c(input$K_slider,
                       input$eta_p_slider,
                       input$eta_v_slider,
                       input$var_u_slider,
                       input$var_v_slider,
                       input$var_p_slider,
                       input$beta_vsat_slider,
                       input$beta_psat_slider)
        fixed     <- c(input$K_fixed,
                       input$eta_p_fixed,
                       input$eta_v_fixed,
                       input$var_u_fixed,
                       input$var_v_fixed,
                       input$var_p_fixed,
                       input$beta_vsat_fixed,
                       input$beta_psat_fixed)
        data.frame(parameter, value, fixed)
    })
    
    observeEvent(input$upload, {
        ext <- tools::file_ext(input$upload$name)
        storage$dataschedule <- switch(ext,
                                     csv = read.csv(input$upload$datapath, sep = ","),
                                     tsv = read.csv(input$upload$datapath, sep = "\t"),
                                     validate("Invalid file; Please upload a .csv or .tsv file")
        )
    })
    
    observeEvent(input$fit, {
        validate(
            need(!is.null(storage$dataschedule), "Please load some data.")
        )
        
        if (is.null(storage$dataschedule)) {
            dataschedule <- data.frame('rotation'=c(rep(0,20),rep(input$rotation_slider,80)),
                                       'trialtype'=rep(1,100))
        } else {
            dataschedule <- storage$dataschedule
        }
        
        fit <- PReMo_fastfit(dataschedule = dataschedule,
                             settings     = modelsettings())
        
        idx <- which(names(fit) == 'value') - 1
        
        for (parname in names(fit)[c(1:idx)]) {
            parvalue <- as.numeric(fit[1,parname])
            if (parname == 'K') {updateSliderInput(inputId = "K_slider", value = parvalue)}
            if (parname == 'eta_p') {updateSliderInput(inputId = "eta_p_slider", value = parvalue)}
            if (parname == 'eta_v') {updateSliderInput(inputId = "eta_v_slider", value = parvalue)}
            if (parname == 'var_u') {updateSliderInput(inputId = "var_u_slider", value = parvalue)}
            if (parname == 'var_p') {updateSliderInput(inputId = "var_p_slider", value = parvalue)}
            if (parname == 'var_v') {updateSliderInput(inputId = "var_v_slider", value = parvalue)}
            if (parname == 'beta_vsat') {updateSliderInput(inputId = "beta_vsat_slider", value = parvalue)}
            if (parname == 'beta_psat') {updateSliderInput(inputId = "beta_psat_slider", value = parvalue)}
        }
    })
    
    observeEvent(input$reset, {
        updateSliderInput(inputId = "K_slider", value = 0.50)
        updateCheckboxInput(inputId = "K_fixed", value = FALSE)
        updateSliderInput(inputId = "eta_p_slider", value = 0.75)
        updateCheckboxInput(inputId = "eta_p_fixed", value = FALSE)
        updateSliderInput(inputId = "eta_v_slider", value = 0.75)
        updateCheckboxInput(inputId = "eta_v_fixed", value = FALSE)
        updateSliderInput(inputId = "var_u_slider", value = 2.00)
        updateCheckboxInput(inputId = "var_u_fixed", value = FALSE)
        updateSliderInput(inputId = "var_p_slider", value = 7.00)
        updateCheckboxInput(inputId = "var_p_fixed", value = FALSE)
        updateSliderInput(inputId = "var_v_slider", value = 0.50)
        updateCheckboxInput(inputId = "var_v_fixed", value = FALSE)
        updateSliderInput(inputId = "beta_vsat_slider", value = 0.50)
        updateCheckboxInput(inputId = "beta_vsat_fixed", value = FALSE)
        updateSliderInput(inputId = "beta_psat_slider", value = 5.00)
        updateCheckboxInput(inputId = "beta_psat_fixed", value = FALSE)
    })
    
    observeEvent(input$reset, {
        storage$dataschedule <- NULL
    })

    
    output$fittedDataPlot <- renderPlot({
        
        if (is.null(storage$dataschedule)) {
            schedule <- data.frame('rotation'=c(rep(0,20),rep(input$rotation_slider,80)),
                                   'trialtype'=rep(1,100))
        } else {
            schedule <- storage$dataschedule[,c('rotation','trialtype')]
        }
        
        # get parameters:
        settings      <- modelsettings()
        parset        <- settings$value
        names(parset) <- settings$parameter
        
        #print(parset)
        
        # run model on schedule & fitted parameters
        modelfit <- PReMo_model(parfree=c(),
                                parset=parset,
                                schedule=schedule)
        
        # determine range of Y values in graph:
        rangevalues <- c(range(schedule$rotation),
                         range(modelfit$reachdev),
                         range(modelfit$proprecal))
        if ('reachdev' %in% names(storage$dataschedule))  { rangevalues <- c(rangevalues, range(storage$dataschedule$reachdev))}
        if ('proprecal' %in% names(storage$dataschedule)) { rangevalues <- c(rangevalues, range(storage$dataschedule$proprecal))}
        
        # make plot:
        plot(-1000,-1000,
             main='',xlab='trial',ylab='deviation',
             xlim=c(0,(dim(schedule)[1]+1)),ylim=range(rangevalues),
             bty='n'
             )
        # add schedule:
        lines(schedule$rotation, col='black')
        
        # add data:
        if (!is.null(storage$dataschedule)) {
            if ('reachdev' %in% names(storage$dataschedule))  { lines(storage$dataschedule$reachdev,  col='#0fd2e2ff') }
            if ('proprecal' %in% names(storage$dataschedule)) { lines(storage$dataschedule$proprecal, col='#ff8200ff') }
        }
        
        # add model processes:
        lines(modelfit$reachdev,  col='#005de4ff')
        lines(modelfit$proprecal, col='#e51636ff')
        lines(modelfit$handest,   col='#ff80ceff')
        
        legend(0,max(rangevalues),
               c('adaptation', 'proprioception', 'integrated hand position'),
               col=c('#005de4ff','#e51636ff','#ff80ceff'),
               lty=c(1,1,1),bty='n')
        
    })
    
    # output$fittedParameters <- renderTable({
    #     storage$parfitted
    # })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
