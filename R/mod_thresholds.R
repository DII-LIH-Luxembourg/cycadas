# Thresholds tab --------------------------------------------------------------
# Table of marker thresholds. Selecting a marker plots its expression; clicking
# the scatter plot moves the threshold and re-assigns the clusters.

mod_thresholds_ui <- function(id) {
  ns <- NS(id)
  tags$div(
    class = "container-fluid py-3",
    uiOutput(ns("notice")),
    layout_columns(
      col_widths = c(5, 7),
      card(card_header("Markers",
                       tags$span(class = "small text-body-secondary fw-normal ms-2",
                                 "Select a marker to adjust its threshold")),
           DTOutput(ns("table"))),
      tagList(
        card(full_screen = TRUE,
             card_header(uiOutput(ns("title"), inline = TRUE)),
             plotOutput(ns("scatter"), click = ns("scatter_click"), height = "300px"),
             card_footer(class = "small text-body-secondary",
                         "Click the plot to move the threshold. Phenotypes are re-assigned immediately.")),
        card(full_screen = TRUE, card_header("Histogram"),
             plotOutput(ns("histogram"), height = "260px"))
      )
    )
  )
}

mod_thresholds_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    th_columns <- c("cell", "threshold", "bi_mod")

    output$notice <- renderUI(if (is.null(state$expr)) no_data_notice())

    output$title <- renderUI({
      if (is.null(input$table_rows_selected) || is.null(state$th)) return("Marker expression")
      sel <- selected()
      tagList(sel$cell,
              tags$span(class = "badge text-bg-light ms-2", sprintf("threshold %.3f", sel$threshold)),
              if (sel$color == "red")
                tags$span(class = "badge text-bg-warning ms-1", "not bimodal"))
    })

    # Re-render the table only when the set of markers changes; threshold edits
    # are pushed through the proxy so that selection and scroll position stay.
    th_cells <- reactiveVal(NULL)
    observe(th_cells(rownames(state$th)))

    output$table <- DT::renderDT({
      req(th_cells())
      DT::datatable(
        isolate(state$th)[, th_columns],
        rownames = FALSE,
        colnames = c("Marker", "Threshold", "Bimodality"),
        class = "compact hover",
        editable = F,
        extensions = c('Buttons', 'Scroller'),
        selection = 'single',
        options = list(
          dom = 'Bfrtip',
          deferRender = TRUE,
          scrollY = 500,
          scroller = TRUE,
          buttons = list(
            list(extend = 'csv', text = "Download CSV", filename = "MarkerThresholds",
                 className = "btn-sm btn-outline-secondary")
          )
        )
      ) %>%
        DT::formatRound("threshold", 3) %>%
        DT::formatStyle("bi_mod", color = DT::styleInterval(BIMODALITY_CUTOFF, c("#e8590c", "inherit")))
    })

    proxy <- DT::dataTableProxy("table")
    observeEvent(state$th, {
      DT::replaceData(proxy, state$th[, th_columns], resetPaging = FALSE, clearSelection = "none",
                      rownames = FALSE)
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
