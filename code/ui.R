# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Projecting plankton gene clusters"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            h4("- Data selection -"),
            radioButtons("pathway",
                         label="Which metabolic pathway to plot to plot?",
                         choices=c("Phosphorylative Oxydation"="00190", "Photosynthesis"="00195")
                         ),
            actionButton("do", "Run the model"),
            h4("- Output parameters -"),
            checkboxInput("do_cor", "Spatial correlation", value = TRUE),
            uiOutput("SliderWidget")
        ), #sidebarPanel

        # Show a plot of the generated distribution
        mainPanel(
           strong("Dimension of the selected data:", style = "font-family : 'arial'; font-si24pt"),
           br(),
           htmlOutput("input_size"),
           br(),
           fluidRow(
               splitLayout(cellWidths = c("50%", "50%"),
                           plotOutput("legend_plot", height = "300px"),
                           plotOutput("correlation_plot", height = "300px"))
           ),
           br(),
           plotOutput("target_map", height = "300px")
        ) # mainPanel
    ) # sidebarLayout
) # ui
