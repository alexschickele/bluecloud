
ui <- fixedPage(theme = shinytheme("sandstone"),
    # ---Application title
    titlePanel("Exploring C4-enzyme projections in picoeukaryotes"),

    # --- Layout type
    sidebarLayout(
      
        # --- Side bar
        sidebarPanel(width = 4,
            selectInput("plot_type",
                        label = "Type of enzyme-level plot:",
                        choices = list(Potential = "unscaled", Realised = "scaled", Difference = "diff")),
            HTML("<h5> <b>Potential:</b> constructed from Connected Component-level where we alleviated the signal corresponding to their relative abundance in the observed data (i.e. alleviating taxonomic effect)",
                 "<br/>",
                 "<b>Realised:</b> constructed from Connected Component-level where we considered the signal corresponding to their relative abundance in the observed data (i.e. considering taxonomic effect)",
                 "<br/>",
                 "<b>Difference:</b> Potential - Realised </h5>",
                 "<br/>"),
            
            tags$hr(style="border-color: black;"),
            selectInput("enz_name",
                        label = "Select the enzyme to plot:",
                        choices = as.list(setNames(c(1:10), names(plot_list)))),
            HTML("<h5> <b>More informations?</b> <br/> A detailed description of the above mentionned enzymes is available at <a>https://www.genome.jp/pathway/map00710</a> </h5> <br/>"),
            
            tags$hr(style="border-color: black;"),
            uiOutput("SliderWidget"),
            HTML("<h5> The projection corresponding to the underlying Connected Components is available in the corresponding tab. They can be browsed using the slidebar above </h5>"),
            
            tags$hr(style="border-color: black;"),
            HTML("<h5> <b>References:</b> <br/> <a>BlueCloud documentation</a> </h5>")
            
        ), # End Side Bar

        # --- Main panel
        mainPanel(
          tabsetPanel(
            # --- Enzyme level tab
            tabPanel("Enzyme - level",
                     fixedRow(
                       splitLayout(cellWidths = c("580px", "200px"),
                                   plotOutput("enzyme_plot", height = "350px"),
                                   plotOutput("legend_plot", height = "250px"))
                     ),
                     HTML("<h5> <b> Description : </b> </h5>",
                          "<h5> Distributional pattern of the selected C4-carbon concentration related enzyme. The uncertainty related to biogeographical projections, 
                                corresponding to algorithm and bootstrap-related variability, is included in the form of a 2D colour scale. 
                                The number of Connected Components (CC) used to construct in each enzyme projection is indicated by the length of the slide bar. 
                                </h5>")
                     
            ), # end enzyme level tab
            
            # --- CC level tab
            tabPanel("Connected component - level",
                     plotOutput("CC_scale_plot", height = "150px", width = "450px"),
                     fixedRow(
                       splitLayout(cellWidths = c("580px", "200px"),
                                   plotOutput("CC_level_map_OLD", height = "350px"),
                                   plotOutput("CC_legend_plot", height = "250px"))
                     ),
                     HTML("<h5> <b> Description : </b> </h5>",
                          "<h5> Distributional pattern of the underlying Connected Component. The uncertainty related to biogeographical projections, 
                                corresponding to algorithm and bootstrap-related variability, is included in the form of a 2D colour scale. The scale
                                bar above refers to the connected component importance in the observed data (i.e., in terms of relative abundance).
                                </h5>"),
                     tags$hr(style="border-color: black;"),
                     htmlOutput("CC_desc")
            ), # CC enzyme level tab
            
          ) # end tabset panel
        ) # mainPanel
    ) # sidebarLayout
) # ui
