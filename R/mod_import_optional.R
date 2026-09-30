# Optional uploads ------------------------------------------------------------
# Marker thresholds, annotation tree, metadata and counts table. Each can be
# loaded on top of an imported dataset.

mod_import_optional_ui <- function(id) {
  ns <- NS(id)
  csv <- c("text/csv", ".csv")
  tags$form(
    id = ns("upload_form"),
    accordion(
      open = FALSE,
      accordion_panel(
        "Marker thresholds", icon = icon("sliders"),
        tags$p(class = "small text-body-secondary",
               "Replaces the estimated thresholds. Existing phenotypes are re-assigned."),
        fileInput(ns("fTH"), NULL, multiple = FALSE, accept = csv)),
      accordion_panel(
        "Annotation tree", icon = icon("sitemap"),
        tags$p(class = "small text-body-secondary",
               "Nodes and edges files from a previous annotation export."),
        fileInput(ns("fNodes"), "Nodes file", multiple = FALSE, accept = csv),
        fileInput(ns("fEdges"), "Edges file", multiple = FALSE, accept = csv),
        actionButton(ns("btnImportTree"), "Import tree", class = "btn-outline-primary btn-sm")),
      accordion_panel(
        "Sample metadata", icon = icon("table-list"),
        tags$p(class = "small text-body-secondary",
               "CSV with the columns sample_id and condition. Needed for differential abundance."),
        fileInput(ns("metadata"), NULL, multiple = FALSE, accept = csv)),
      accordion_panel(
        "Cluster counts", icon = icon("table-cells"),
        tags$p(class = "small text-body-secondary",
               "Cells per cluster (rows) and sample (columns). Needed for differential abundance."),
        fileInput(ns("counts_table"), NULL, multiple = FALSE, accept = csv))
    )
  )
}

mod_import_optional_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    # Marker thresholds ----
    observeEvent(input$fTH, {
      th <- prepare_thresholds(read.csv(input$fTH$datapath))
      if (!is.null(state$expr)) th <- add_estimates(th, state$expr, state$markers)
      state$th <- th
      # re-calculate the tree after threshold upload
      if (!is.null(state$expr)) rebuild_annotation(state)
    })

    # Annotation tree ----
    observeEvent(input$btnImportTree, {
      req(input$fNodes, input$fEdges)
      if (is.null(state$expr)) {
        showNotification("Import expression data before the annotation tree.", type = "error")
        return()
      }
      state$graph <- getGraphFromLoad(read.csv(input$fNodes$datapath),
                                      read.csv(input$fEdges$datapath))
      rebuild_annotation(state)
    })

    # Metadata and counts ----
    observeEvent(input$metadata, {
      state$md <- drop_index_column(read.csv(input$metadata$datapath))
    })
    observeEvent(input$counts_table, {
      state$counts <- drop_index_column(read.csv(input$counts_table$datapath))
    })

    observeEvent(state$reset_version, shinyjs::reset("upload_form"), ignoreInit = TRUE)
  })
}
