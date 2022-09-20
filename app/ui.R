
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
            p(tags$b("Potential"), "constructed from Connected Component-level where we alleviated the signal corresponding to their relative abundance in the observed data (i.e. alleviating taxonomic effect)"),
            p(tags$b("Realised:"), "constructed from Connected Component-level where we considered the signal corresponding to their relative abundance in the observed data (i.e. considering taxonomic effect)")),
            p(tags$b("Difference:"), "Potential - Realised")),
            
            tags$hr(style="border-color: black;"),
            selectInput("enz_name",
                        label = "Select the enzyme to plot:",
                        choices = as.list(setNames(c(1:10), names(plot_list)))),
            p(tags$b("More information?"), tags$br(), "A detailed description of the above mentionned enzymes is available in", a(href="https://www.genome.jp/pathway/map00710", "the KEGG browser")),
            
            tags$hr(style="border-color: black;"),
            # NB JOI: je pense que ce morceau devrait n'être que dans la tab correspondante. Tu peux mettre des bouts d'UI dans une tab, tout n'a pas à être dans cette colonne de gauche.
            uiOutput("SliderWidget"),
            p("The projection corresponding to the underlying Connected Components is available in the corresponding tab. They can be browsed using the slidebar above"),
            
            tags$hr(style="border-color: black;"),
            h5("Going further")
            p("You can generate maps for other functions using the PlanktonGenomics virtual lab of the ", a(href="https://blue-cloud.org", "Blue Cloud project") ". To do so, ", a(href="https://blue-cloud.d4science.org/web/planktongenomics", "register to the VLab"), " and ", a(href="?? ] -> link to the handbook PDF", "read the documentation"), ".")
            
        ), # End Side Bar

        # --- Main panel
        mainPanel(
          tabsetPanel(
            # --- Enzyme level tab
            tabPanel("Enzyme - level",
                     plotOutput("enzyme_plot", height = "430px"),
                     fixedRow(
                       column(width = 4, plotOutput("legend_plot", height = "250px")
                              ),
                       column(width = 7, HTML("<h5> <b> Description : </b> </h5>",
                                              "<h5> Distributional pattern of the selected C4-carbon concentration related enzyme. The uncertainty related to biogeographical projections, 
                                corresponding to algorithm and bootstrap-related variability, is included in the form of a 2D colour scale. 
                                The number of Connected Components (CC) used to construct in each enzyme projection is indicated by the length of the slide bar. 
                                </h5>")
                              )
                     ),
            ), # end enzyme level tab
            
            # --- CC level tab
            tabPanel("Connected component - level",
                     plotOutput("CC_scale_plot", height = "100px"),
                     plotOutput("CC_level_map_OLD", height = "430px"),
                     fixedRow(
                       column(width = 4, plotOutput("CC_legend_plot", height = "250px")
                              ),
                       column(width = 7, HTML("<h5> <b> Description : </b> </h5>",
                                              "<h5> Distributional pattern of the underlying Connected Component. The uncertainty related to biogeographical projections, 
                                corresponding to algorithm and bootstrap-related variability, is included in the form of a 2D colour scale. The scale
                                bar above refers to the connected component importance in the observed data (i.e., in terms of relative abundance).
                                </h5>")
                              )
                     ),
                     tags$hr(style="border-color: black;"),
                     htmlOutput("CC_desc")
            ), # CC enzyme level tab
            
          ) # end tabset panel
        ) # mainPanel
    ) # sidebarLayout
) # ui
