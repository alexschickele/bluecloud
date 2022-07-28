
ui <- fluidPage(

    # Application title
    titlePanel("Exploring C4-enzyme projections in picoeukaryotes"),

    # Layout
    sidebarLayout(
      
        # --- Side bar
        sidebarPanel(

            selectInput("plot_type",
                        label = "Type of enzyme-level plot:",
                        choices = list(Potential = "unscaled", Realised = "scaled", Difference = "diff")),
            
            selectInput("enz_name",
                        label = "Select the enzyme to plot:",
                        choices = as.list(setNames(c(1:10), names(plot_list)))),
            
            uiOutput("SliderWidget")
        ), #sidebarPanel

        # --- Main panel
        mainPanel(
           fluidRow(
               splitLayout(cellWidths = c("200px", "500px"),
                           plotOutput("legend_plot", height = "300px"),
                           plotOutput("enzyme_plot", height = "300px"))
           ),
           hr(style = "border-top: 3px solid #000000;"),
           fluidRow(
             splitLayout(cellWidths = c("200px", "500px"),
                         plotOutput("CC_legend_plot", height = "300px"),
                         plotOutput("CC_scale_plot", height = "150px"))
           ),
           plotOutput("CC_level_map", height = "400px", width = "700px"),
           htmlOutput("CC_desc")

        ) # mainPanel
    ) # sidebarLayout
) # ui
