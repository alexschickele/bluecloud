# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            radioButtons("column",
                         label="Which data to plot?",
                         choices=c("eruption time"="eruptions", "waiting time"="waiting")
                         ),
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30),
            radioButtons("clr",
                         label="Colour for the bars",
                         choices=c("grey", "blue", "red"),
                         selected="grey")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)
