# Thresholds tab --------------------------------------------------------------
# Table of marker thresholds. Selecting a marker plots its expression; clicking
# the scatter plot moves the threshold and re-assigns the clusters.

mod_thresholds_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(width = 6,
           box(width = NULL, title = "Marker Expression:",
               plotOutput(ns("scatter"), click = ns("scatter_click"))),
           box(width = NULL, title = "Histogram:",
               plotOutput(ns("histogram")))
    ),
    column(width = 6,
           box(width = NULL, DTOutput(ns("table"))))
  )
}

mod_thresholds_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    th_columns <- c("cell", "threshold", "bi_mod")

    # Re-render the table only when the set of markers changes; threshold edits
    # are pushed through the proxy so that selection and scroll position stay.
    th_cells <- reactiveVal(NULL)
    observe(th_cells(rownames(state$th)))

    output$table <- DT::renderDT({
      req(th_cells())
      DT::datatable(
        isolate(state$th)[, th_columns],
        editable = F,
        extensions = c('Buttons', 'Scroller'),
        selection = 'single',
        options = list(
          dom = 'Bfrtip',
          deferRender = TRUE,
          scrollY = 500,
          scroller = TRUE,
          buttons = list(
            list(extend = 'csv', filename = "MarkerThresholds")
          )
        )
      )
    })

    proxy <- DT::dataTableProxy("table")
    observeEvent(state$th, {
      DT::replaceData(proxy, state$th[, th_columns], resetPaging = FALSE, clearSelection = "none")
    })

    selected <- reactive({
      req(state$th, input$table_rows_selected)
      state$th[input$table_rows_selected, ]
    })

    # 0-1 expression of the selected marker with a fixed random jitter
    points <- reactive({
      req(state$expr)
      x <- state$expr[, selected()$cell]
      set.seed(1)
      data.frame(x = x, y = rnorm(length(x)))
    })

    output$scatter <- renderPlot({
      plot_threshold_scatter(points(), selected()$threshold, selected()$color)
    })

    output$histogram <- renderPlot({
      plot_threshold_histogram(points(), selected()$threshold, selected()$color)
    })

    observeEvent(input$scatter_click, {
      req(input$table_rows_selected)
      state$th[input$table_rows_selected, "threshold"] <- round(input$scatter_click$x, 3)
      rebuild_annotation(state)
    })
  })
}
