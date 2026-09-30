# Explore Data tab ------------------------------------------------------------
# UMAP coloured by marker expression, and a brushable UMAP that shows the
# selected clusters as heatmap and table.

mod_explore_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(column(width = 10,
                    box(width = NULL, plotOutput(ns("umap_marker")))),
             column(width = 2,
                    box(width = NULL, selectInput(ns("markerSelect"), "Select:", choices = NULL)))),
    fluidRow(column(width = 6,
                    box(width = NULL, plotOutput(ns("umap"), brush = ns("umap_brush")))),
             column(width = 6,
                    box(width = NULL, plotOutput(ns("heatmap"))))),
    fluidRow(column(width = 12,
                    box(width = NULL, DTOutput(ns("umap_data")))))
  )
}

mod_explore_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    observeEvent(state$markers, {
      updateSelectInput(session, "markerSelect", "Select:", state$markers)
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
        pheatmap(round(sel, 2) %>% dplyr::select(-c("cluster_number", "u1", "u2")), cluster_cols = F)
      } else {
        ggplot() + theme_void() + ggtitle("Select area on Umap to plot Heatmap")
      }
    })

    output$umap_data <- DT::renderDT(server = FALSE, {
      sel <- brushed()
      if (nrow(sel) == 0) {
        return(DT::datatable(
          data.frame(Message = "No points selected"),
          rownames = FALSE,
          options = list(dom = 't', paging = FALSE)
        ))
      }

      DT::datatable(
        round(sel, 2),
        filter = 'top',
        extensions = 'Buttons',
        options = list(
          scrollY = 600,
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
