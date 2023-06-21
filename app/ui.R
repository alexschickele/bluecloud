
ui <- fixedPage(theme = shinytheme("sandstone"),
    # ---Application title
    titlePanel("Exploring the genomic potential for C4-photosynthesis in picoeukaryotes"),

    # --- Layout type
    sidebarLayout(
      
        # --- Side bar
        sidebarPanel(width = 4,
                     # --- Enzyme name
                     selectInput("enz_name",
                                 label = "Select the enzyme to plot:",
                                 choices = as.list(setNames(c(1:10), names(plot_list)))),
                     p(tags$b("More information?"), tags$br(), "A detailed description of the above mentionned enzymes is available in", a(href="https://www.genome.jp/pathway/map00710", "the KEGG browser")),
                     
                     tags$hr(style="border-color: black;"),
                     
                     # --- Conditional enzyme-level plot type
                     conditionalPanel(condition = "input.tabselected==1",
                                      selectInput("plot_type",
                                                  label = "Type of enzyme-level plot:",
                                                  choices = list(Standardized = "unscaled", Weighted = "scaled", Difference = "diff")),
                                      p(tags$b("Standardized"), "constructed from protein functional cluster-level where we alleviated the signal corresponding to their relative abundance in the observed data (i.e. alleviating taxonomic effect)"),
                                      p(tags$b("Weighted:"), "constructed from protein functional cluster-level where we considered the signal corresponding to their relative abundance in the observed data (i.e. considering taxonomic effect)"),
                                      p(tags$b("Difference:"), "Standardized - Weighted"),
                                      
                                      tags$hr(style="border-color: black;")),
                     
                     # --- Conditional CC-level selection bar
                     conditionalPanel(condition = "input.tabselected==2",
                             uiOutput("SliderWidget"),
                             p("The projection corresponding to the underlying Protein Functional Cluster is available in the corresponding tab. They can be browsed using the slidebar above"),
                             tags$hr(style="border-color: black;")),
                     
                     # --- Going further
                     p(tags$b("Going further:"), tags$br(),
                       "You can generate maps for other metabolic functions using the PlanktonGenomics virtual lab of the ", 
                       a(href="https://blue-cloud.org", "Blue Cloud project"), ". To do so, ", 
                       a(href="https://blue-cloud.d4science.org/web/planktongenomics", "register to the VLab"), ", read the ", 
                       a(href="https://data.d4science.net/qa7Z", "catalogue documentation"), ", and find corresponding files",
                       a(href="https://data.d4science.net/3hxw", "here"), ".")

        ), # End Side Bar

        # --- Main panel
        mainPanel(
          tabsetPanel(id = "tabselected",
                      
            # --- Enzyme level tab
            tabPanel("Enzyme - level",  value = 1,
                     plotOutput("enzyme_plot", height = "430px"),
                     fixedRow(
                       column(width = 4, plotOutput("legend_plot", height = "250px")
                              ),
                       column(width = 7, HTML("<h5> <b> Description : </b> </h5>",
                                              "<h5> Distributional pattern of the selected C4-enzyme. The uncertainty related to projections, 
                                corresponding to algorithm and bootstrap-related variability, is included in the form of a 2D colour scale. 
                                The number of Protein Functional Clusters used to construct in each enzyme projection is indicated by the length of the slide bar. 
                                </h5>")
                              )
                     ),
            ), # end enzyme level tab
            
            # --- CC level tab
            tabPanel("Protein Functional Cluster - level", value = 2,
                     plotOutput("CC_scale_plot", height = "100px"),
                     plotOutput("CC_level_map_OLD", height = "430px"),
                     fixedRow(
                       column(width = 4, plotOutput("CC_legend_plot", height = "250px")
                              ),
                       column(width = 7, HTML("<h5> <b> Description : </b> </h5>",
                                              "<h5> Distribution pattern of the underlying Protein Functional cluster. The uncertainty related to projections, 
                                corresponding to algorithm and bootstrap-related variability, is included in the form of a 2D colour scale. The scale
                                bar above refers to the protein functional cluster importance in the observed data (i.e., in terms of relative abundance).
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
