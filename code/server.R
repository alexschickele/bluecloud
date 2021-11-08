# Define server logic required to project the maps
# 1) query on the database - triggered by input pathway + execute button
# 2) output = number of station, features and targets selected
# 3) run model learning - triggered by input pathway + execute button
# 4) run bootstrap - triggered by input pathway + execute button
# 5) output = legend plot in a panel
# 6) output = reactive plot of maps - triggered by the chosen target number (slidebar)

# Example runs :
# A) Change the pathway -> step 1-5 done when execute button triggered
# B) Change the target to plot -> step 6 immediately

# TO DO LIST:
# - fix file path dynamically if possible
# - double check reactivity
# - check how to display console for the python progression bar on run and proj

server <- function(input, output) {
  # --- PART 1 : querying the database and displaying data size ---
  query <- reactive({
    message("--- querying database ---")
    query_data(KEGG_p = input$pathway)
  }) # end event
  
  output$input_size <- renderUI({
    df <- req(query())
    message("--- display output df ---")
    nstations = paste("Number of stations :", nrow(df$X))
    nfeatures = paste("Number of features :", ncol(df$X))
    ntargets = paste("Number of targets :", ncol(df$Y))
    HTML(paste(nstations, nfeatures, ntargets, sep = '<br/>'))
  }) #end output
  
  # --- PART 2 : run the model learning ---
  query_time <- eventReactive(input$do,{
    message("--- get query time ---")
    file.info("/home/aschickele/workspace/bluecloud descriptor/data/X.feather")$ctime
  })
  
  run <- observeEvent(query_time(), {
    req(query_time())
    message("--- running model learning ---")
    model_run()
  }) # end event
  
  # --- PART 3 : run the projection calculations ---
  run_time <- eventReactive(input$do,{
    message("--- get run time ---")
    file.info("/home/aschickele/workspace/bluecloud descriptor/data/m")$ctime
  })
  
  proj <- eventReactive(run_time(), {
    req(run_time())
    message("--- building projections ---")
    model_proj()
  }) # end event
  
  # --- FINAL OUTPUTS : fix legend + dynamic maps
  output$legend_plot <- renderPlot({
    p <- req(proj())
    message("--- plotting legend")
    legend_proj(col_matrix = p$col_matrix, cutx = p$cutx, cuty = p$cuty)
  })
  
  output$correlation_plot <- renderPlot({
    req(input$do_cor)
    p <- req(proj())
    message("--- plotting spatial correlations ---")
    cor_proj(proj = p$proj)
  })
  
  output$SliderWidget <- renderUI({
    p <- req(proj())
    message("--- Add slidebar to select target map")
    sliderInput("map_nb", "Target number to plot:", 
                min = 1, max = nlayers(p$proj), value = 1, step = 1)
  })
  
  output$target_map <- renderPlot({
    p <- req(proj())
    message("--- selecting target map")
    map_proj(proj = p$proj, col = p$col, targetID = input$map_nb)
  })

} # end server

