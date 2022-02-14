library(shiny)
library(optimx)

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
            
            actionButton("fit", "Fit model to data", class = "btn-lg btn-success"),
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            p("Data input needs to be a .csv or .tsv file, where rows are trials, and with 2 columns defining the experiment: `rotation` (with degrees rotation) and `trialtype` either 1: normal cursor, or 2: error-clamped feedback. It also need to have at least one of two columns: `reachdev` with reach deviations (in degrees) or `proprecal` with proprioceptive estimates. All measures relative to the target."),
            p("Parameters can be set by making the minimum and maximum the same. If not set, they will be fit using the minimum and maximum as lower and upper limits.")
            # plotOutput("distPlot")
            # tableOutput("fittedParameters")
        )
    )
)

# Define server logic 
server <- function(input, output) {
    
    dataschedule <- reactive({
        req(input$upload)
        
        ext <- tools::file_ext(input$upload$name)
        switch(ext,
               csv = read.csv(input$upload$datapath, sep = ","),
               tsv = read.csv(input$upload$datapath, sep = "\t"),
               validate("Invalid file; Please upload a .csv or .tsv file")
        )
    })
    
    
    
    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
