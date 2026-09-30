# DA interactive Tree tab -----------------------------------------------------
# Click a node of the annotation tree to test the proportion of the node and
# all its descendants between conditions.

mod_da_tree_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    tags$head(tags$style(sprintf("#%s{height:700px !important;}", ns("treebox")))),
    column(width = 8,
           box(id = ns("treebox"), width = NULL, title = "Interactive DA Tree",
               visNetworkOutput(ns("tree"), width = "100%", height = "700px"))),
    column(width = 4,
           box(width = NULL, tableOutput(ns("DA_table"))),
           box(width = NULL, plotOutput(ns("boxplot"))))
  )
}

mod_da_tree_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    output$tree <- renderVisNetwork({
      req(state$graph)
      tree_network(state$graph, session$ns("tree_click"))
    })

    node_id <- reactiveVal(NULL)
    observeEvent(input$tree_click, node_id(input$tree_click[1]))
    observeEvent(state$dataset_version, node_id(NULL), ignoreInit = TRUE)

    node_label <- reactive({
      req(node_id())
      state$graph$nodes$label[state$graph$nodes$id == node_id()]
    })

    props <- reactive({
      req(node_id(), state$counts, state$md, state$expr)
      node_proportions(state$counts, state$expr$cell, state$graph, node_id(), state$md)
    })

    output$DA_table <- renderTable({
      # the root holds all clusters, there is nothing to compare
      req(node_id() != 1)
      node_da(props(), node_label())
    })

    output$boxplot <- renderPlot({
      if (is.null(node_id())) {
        plot_message("No Data Available")
      } else {
        plot_da_boxplot(props(), node_label())
      }
    })
  })
}
