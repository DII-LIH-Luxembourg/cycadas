# Explore Data tab ------------------------------------------------------------
# UMAP coloured by marker expression, and a brushable UMAP that shows the
# selected clusters as heatmap and table.

mod_explore_ui <- function(id) {
  ns <- NS(id)
  tags$div(
    class = "container-fluid py-3",
    uiOutput(ns("notice")),
    layout_columns(
      col_widths = c(6, 6),
      card(full_screen = TRUE,
           card_header(class = "d-flex align-items-center justify-content-between",
                       "Marker expression",
                       tags$div(style = "width: 14rem; font-weight: normal;",
                                pickerInput(ns("markerSelect"), NULL, choices = NULL,
                                            options = picker_options(), width = "100%"))),
           plotOutput(ns("umap_marker"), height = "460px")),
      card(full_screen = TRUE,
           card_header("Select clusters",
                       tags$span(class = "small text-body-secondary fw-normal ms-2",
                                 "Drag on the UMAP to select an area")),
           plotOutput(ns("umap"), brush = ns("umap_brush"), height = "460px"))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(full_screen = TRUE, card_header("Heatmap of the selected clusters"),
           plotOutput(ns("heatmap"), height = "460px")),
      card(full_screen = TRUE, card_header("Selected clusters"),
           DTOutput(ns("umap_data")))
    )
  )
}

mod_explore_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    output$notice <- renderUI(if (is.null(state$expr)) no_data_notice())

    observeEvent(state$markers, {
      updatePickerInput(session, "markerSelect", choices = state$markers)
    }, ignoreNULL = FALSE)

    brushed <- reactive({
      req(state$umap)
      brushedPoints(state$umap, input$umap_brush)
    })

    output$umap_marker <- renderPlot({
      req(state$umap, input$markerSelect %in% state$markers)
      plot_umap_marker(state$umap, input$markerSelect)
    })

    output$umap <- renderPlot({
      req(state$umap)
      plot_umap(state$umap)
    })

    output$heatmap <- renderPlot({
      sel <- brushed()
      if (nrow(sel) > 0) {
        plot_cluster_heatmap(sel %>% dplyr::select(-c("cluster_number", "u1", "u2")))
      } else {
        plot_message("Select an area on the UMAP")
      }
    })

    output$umap_data <- DT::renderDT(server = FALSE, {
      sel <- brushed()
      if (nrow(sel) == 0) {
        return(DT::datatable(
          data.frame(Message = "No clusters selected"),
          rownames = FALSE,
          options = list(dom = 't', paging = FALSE)
        ))
      }

      DT::datatable(
        round(sel, 2),
        filter = 'top',
        extensions = 'Buttons',
        options = list(
          scrollY = 360,
          scrollX = TRUE,
          dom = '<"float-left"l><"float-right"f>rt<"row"<"col-sm-4"B><"col-sm-4"i><"col-sm-4"p>>',
          lengthMenu = list(c(10, 25, 50, -1), c('10', '25', '50', 'All')),
          scrollCollapse = TRUE,
          lengthChange = TRUE,
          widthChange = TRUE,
          rownames = TRUE
        )
      )
    })
  })
}
