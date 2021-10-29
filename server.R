# Define server logic required to draw a histogram
server <- function(input, output) {

  get_data <- reactive({
    message("get data")
    # generate bins based on input$bins from ui.R
    my_data[,input$column]
  })

  define_bins <- reactive({
    message("define bins")
    x <- get_data()
    seq(min(x), max(x), length.out = input$bins + 1)
  })

  output$distPlot <- renderPlot({
    message("draw the histogram")
    x <- get_data()
    bins <- define_bins()
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = input$clr, border = 'white')
  })


}
