# Thresholds tab --------------------------------------------------------------
# Table of marker thresholds. Selecting a marker plots its expression; clicking
# the scatter plot moves the threshold and re-assigns the clusters. The
# estimated threshold is kept, shown next to the current one and can be
# restored.

mod_thresholds_ui <- function(id) {
  ns <- NS(id)
  tags$div(
    class = "container-fluid py-3",
    uiOutput(ns("notice")),
    layout_columns(
      col_widths = c(5, 7),
      card(card_title(tagList("Markers", uiOutput(ns("n_modified"), inline = TRUE)),
                      actionButton(ns("resetAll"), "Reset all", icon = icon("rotate-left"),
                                   class = "btn-outline-secondary btn-sm"),
                      downloadButton(ns("download"), "CSV", class = "btn-outline-secondary btn-sm")),
           DTOutput(ns("table")),
           card_footer(class = "small text-body-secondary",
                       "Bold thresholds were changed by hand. Select a marker to adjust its threshold.")),
      tagList(
        card(full_screen = TRUE,
             card_title(uiOutput(ns("title"), inline = TRUE),
                        actionButton(ns("resetOne"), "Reset to estimate", icon = icon("rotate-left"),
                                     class = "btn-outline-secondary btn-sm")),
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

    output$notice <- renderUI(if (is.null(state$expr)) no_data_notice())

    # Table data: `modified` is a hidden column used for styling
    table_data <- function(th) {
      data.frame(cell = th$cell, threshold = th$threshold, estimated = th$estimated,
                 bi_mod = th$bi_mod, modified = threshold_modified(th))
    }

    # Re-render the table only when the set of markers changes; threshold edits
    # are pushed through the proxy so that selection and scroll position stay.
    th_cells <- reactiveVal(NULL)
    observe(th_cells(rownames(state$th)))

    output$table <- DT::renderDT({
      req(th_cells())
      DT::datatable(
        table_data(isolate(state$th)),
        rownames = FALSE,
        colnames = c("Marker", "Threshold", "Estimated", "Bimodality", "modified"),
        class = "compact hover",
        editable = F,
        extensions = "Scroller",
        selection = 'single',
        options = list(
          dom = 'frtip',
          deferRender = TRUE,
          scrollY = 500,
          scroller = TRUE,
          columnDefs = list(list(targets = 4, visible = FALSE))
        )
      ) %>%
        DT::formatRound(c("threshold", "estimated"), 3) %>%
        DT::formatStyle("threshold", valueColumns = "modified",
                        fontWeight = DT::styleEqual(c(FALSE, TRUE), c("normal", "bold")),
                        color = DT::styleEqual(c(FALSE, TRUE), c("inherit", COLOR_SELECTED))) %>%
        DT::formatStyle("estimated", color = COLOR_ESTIMATE) %>%
        DT::formatStyle("bi_mod", color = DT::styleInterval(BIMODALITY_CUTOFF, c("#e8590c", "inherit")))
    })

    proxy <- DT::dataTableProxy("table")
    observeEvent(state$th, {
      DT::replaceData(proxy, table_data(state$th), resetPaging = FALSE, clearSelection = "none",
                      rownames = FALSE)
    })

    n_modified <- reactive(if (is.null(state$th)) 0 else sum(threshold_modified(state$th)))

    output$n_modified <- renderUI({
      if (n_modified() > 0)
        tags$span(class = "badge text-bg-light ms-2 fw-normal", sprintf("%d changed", n_modified()))
    })

    selected <- reactive({
      req(state$th, input$table_rows_selected)
      state$th[input$table_rows_selected, ]
    })

    observe({
      shinyjs::toggleState("resetAll", condition = n_modified() > 0)
      sel_modified <- !is.null(input$table_rows_selected) && !is.null(state$th) &&
        isTRUE(threshold_modified(state$th[input$table_rows_selected, ]))
      shinyjs::toggleState("resetOne", condition = sel_modified)
    })

    output$title <- renderUI({
      if (is.null(input$table_rows_selected) || is.null(state$th)) return("Marker expression")
      sel <- selected()
      tagList(sel$cell,
              tags$span(class = "badge text-bg-light ms-2", sprintf("threshold %.3f", sel$threshold)),
              if (!is.na(sel$estimated))
                tags$span(class = "badge text-bg-light ms-1 fw-normal",
                          sprintf("estimate %.3f", sel$estimated)),
              if (sel$color == "red")
                tags$span(class = "badge text-bg-warning ms-1", "not bimodal"))
    })

    # 0-1 expression of the selected marker with a fixed random jitter
    points <- reactive({
      req(state$expr)
      x <- state$expr[, selected()$cell]
      set.seed(1)
      data.frame(x = x, y = rnorm(length(x)))
    })

    output$scatter <- renderPlot({
      sel <- selected()
      plot_threshold_scatter(points(), sel$threshold, sel$color, sel$estimated)
    })

    output$histogram <- renderPlot({
      sel <- selected()
      plot_threshold_histogram(points(), sel$threshold, sel$color, sel$estimated)
    })

    observeEvent(input$scatter_click, {
      req(input$table_rows_selected)
      state$th[input$table_rows_selected, "threshold"] <- round(input$scatter_click$x, 3)
      rebuild_annotation(state)
    })

    # Reset to estimates ------------------------------------------------------
    observeEvent(input$resetOne, {
      state$th <- reset_thresholds(state$th, selected()$cell)
      rebuild_annotation(state)
    })

    observeEvent(input$resetAll, {
      showModal(modalDialog(
        title = "Reset all thresholds",
        sprintf("Set %d changed thresholds back to their estimates? Phenotypes are re-assigned.",
                n_modified()),
        easyClose = TRUE,
        footer = tagList(modalButton("Cancel"),
                         actionButton(session$ns("confirmResetAll"), "Reset", class = "btn-primary"))
      ))
    })

    observeEvent(input$confirmResetAll, {
      removeModal()
      state$th <- reset_thresholds(state$th)
      rebuild_annotation(state)
    })

    output$download <- downloadHandler(
      filename = function() paste0("MarkerThresholds_", Sys.Date(), ".csv"),
      content = function(file) write_thresholds(state$th, file)
    )
  })
}
