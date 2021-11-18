# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Projecting plankton gene clusters"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            h4("1. Description :"),
            p("Based on a machine learning framework, this ready-to-use tool aims to
              explore the biogeography of plankton genetic diversity related to key metabolic pathways."),
            h4("2. Data selection :"),
            p("Please, select a metabolic pathway among the following choices. The extracted data correspond
              to the abundance of corresponding gene clusters and the environmnetal conditions at each TARA Ocean stations."),
            radioButtons("pathway",
                         label="Which metabolic pathway to model?",
                         choices=c("Phosphorylative Oxydation"="00190", 
                                   "Photosynthesis"="00195",
                                   "Carbon fixation"="00710",
                                   "Methane metabolism"="00680",
                                   "Nitrogen metabolism"="00910",
                                   "Sulfur metabolism"="00920")
                         ),
            p("To model the selected metabolic pathway, please click below. Note that this may take a few minuts to process"),
            actionButton("do", "Run the model"),
            h4("3. Output parameters :"),
            p("To display the spatial correlation among gene clusters corresponding to the selected pathway, tick the mark below:"),
            checkboxInput("do_cor", "Spatial correlation", value = TRUE),
            p("The gene clusters corresponding to the selected pathway can be selected using the sliding bar below.
              Note that a functional description of the cluster is provided below the map."),
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
           plotOutput("target_map", height = "300px"),
           br(),
           htmlOutput("target_description")
        ) # mainPanel
    ) # sidebarLayout
) # ui
