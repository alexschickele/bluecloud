
server <- function(input, output) {
  
  # =============== ENZYME LEVEL ===================
  # --- Selecting aggregated map data (raster)
  enz_r <- reactive({
    message("--- updating raster data")
    get(paste0("func_r_", input$plot_type))
  }) # end reactive
  
  # --- Selecting aggregated map data (color palette)
  enz_pal <- reactive({
    message("--- updating color palette")
    get(paste0("func_r_pal_", input$plot_type))
  }) # end reactive

  # --- Plotting color scale legend
  output$legend_plot <- renderPlot({
    message(" --- plotting legend")
    legend_proj(col_matrix = proj$col_matrix)
  }) # end output
  
  # --- Plotting aggregated map
  output$enzyme_plot <- renderPlot({
    message(" --- plotting aggregated map")
    enz_data <- list(r = enz_r(), pal = enz_pal())
    plot(enz_data$r[[as.numeric(input$enz_name)]], col = enz_data$pal[[as.numeric(input$enz_name)]], legend = FALSE)
  }) # end output

  # ========= CONNECTED COMPONENT LEVEL =============
  # --- Get raster stack ID of corresponding CC
  id_r <- reactive({
    message("  --- updating corresponding raster ID")
    CC_desc_e$pos_nn_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[as.numeric(input$enz_name)]])==TRUE)]
  }) # end reactive
  
  # --- Creating slide bar
  output$SliderWidget <- renderUI({
    message("  --- update SliderWidget")
    tmp <- id_r()
    sliderInput("CC_to_plot", "Explore underlying Connected Component projection:",
                min = 1, max = length(tmp), value = 1, step = 1)
  }) # end UI output
  
  # --- Plot CC_legend
  output$CC_legend_plot <- renderPlot({
    message(" --- plotting CC_legend")
    legend_proj(col_matrix = proj$col_matrix)
  }) # end output
  
  # --- Plot CC_scale
  output$CC_scale_plot <- renderPlot({
    message(" --- plotting CC_scale")
    tmp <- CC_desc_e$sum_CC[which(str_detect(CC_desc_e$kegg_ko, plot_list[[as.numeric(input$enz_name)]])==TRUE)]
    map_scale <- tmp[input$CC_to_plot]/max(tmp)
    barplot(map_scale, xlim = c(0,1), border = NA, col = "black", horiz = TRUE,
            main = "Scale :")
  }) # end output
  
  # --- Plot CC_level_map
  output$CC_level_map <- renderPlot({
    message("  --- plotting CC_level_map")
    tmp <- id_r()
    tmp <- tmp[input$CC_to_plot]
    plot(proj$proj[[tmp]], col = proj$col[[tmp]], legend = FALSE)
  }) # end output
  
  # --- Add CC_desc
  output$CC_desc <- renderUI({
    tmp <- which(str_detect(CC_desc_e$kegg_ko, plot_list[[as.numeric(input$enz_name)]])==TRUE)[input$CC_to_plot]
    HTML(paste("<u><b>Connected Component description:</b></u>", "<br/>"),
         paste("<b>KEGG KO:</b>", CC_desc_e$kegg_ko[tmp], "<br/>"),
         paste("<b>Unknown rate:</b>", round(CC_desc_e$unknown_rate[tmp],2), "%","<br/>"),
         paste("<b>Taxonomic Class:</b>", CC_desc_e$class[tmp], "<br/>"),
         paste("<b>Taxonomic Genus:</b>", "<i>",CC_desc_e$genus[tmp], "</i>", "<br/>")
    ) # end HTML
  }) # end output
  
} # end server

