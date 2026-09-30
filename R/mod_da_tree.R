# DA interactive Tree tab -----------------------------------------------------
# Click a node of the annotation tree to test the proportion of the node and
# all its descendants between conditions.

mod_da_tree_ui <- function(id) {
  ns <- NS(id)
  tags$div(
    class = "container-fluid py-3",
    uiOutput(ns("notice")),
    layout_columns(
      col_widths = c(8, 4),
      card(full_screen = TRUE,
           card_header("Annotation tree",
                       tags$span(class = "small text-body-secondary fw-normal ms-2",
                                 "Click a node to compare its proportion between conditions")),
           card_body(class = "p-1", visNetworkOutput(ns("tree"), width = "100%", height = "640px"))),
      tagList(
        card(full_screen = TRUE, card_header(textOutput(ns("title"), inline = TRUE)),
             plotOutput(ns("boxplot"), height = "380px")),
        card(card_header("Pairwise tests"), tableOutput(ns("DA_table")))
      )
    )
  )
}

mod_da_tree_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    output$notice <- renderUI({
      if (is.null(state$expr)) no_data_notice()
      else if (!has_df(state$md) || !has_df(state$counts))
        tags$div(class = "alert alert-secondary",
                 "Load sample metadata and cluster counts on the Workspace tab to compare conditions.")
    })

    output$tree <- renderVisNetwork({
      req(state$graph, state$expr)
      tree_network(state$graph, session$ns("tree_click"), expr = state$expr)
    })

    output$title <- renderText(if (is.null(node_id())) "Proportion per condition" else node_label())

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
        plot_message("Click a node in the tree")
      } else {
        plot_da_boxplot(props(), node_label())
      }
    })
  })
}
