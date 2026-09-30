# Optional uploads ------------------------------------------------------------
# Marker thresholds, annotation tree, metadata and counts table. Each can be
# loaded on top of an imported dataset.

mod_import_optional_ui <- function(id) {
  ns <- NS(id)
  csv <- c("text/csv", ".csv")
  tags$form(
    id = ns("upload_form"),
    settings_box("Optional - Marker-thresholds", status = "warning", collapsed = TRUE,
                 fileInput(ns("fTH"), "Choose CSV File", multiple = FALSE, accept = csv)),
    settings_box("Optional - Annotation Tree", status = "warning", collapsed = TRUE,
                 fileInput(ns("fNodes"), "Choose Nodes File", multiple = FALSE, accept = csv),
                 fileInput(ns("fEdges"), "Choose Edges File", multiple = FALSE, accept = csv),
                 actionButton(ns("btnImportTree"), "Import")),
    settings_box("Optional - Metadata", status = "warning", collapsed = TRUE,
                 fileInput(ns("metadata"), "Choose CSV File", multiple = F, accept = csv)),
    settings_box("Optional - Count Table", status = "warning", collapsed = TRUE,
                 fileInput(ns("counts_table"), "Choose CSV File", multiple = F, accept = csv))
  )
}

mod_import_optional_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    # Marker thresholds ----
    observeEvent(input$fTH, {
      state$th <- prepare_thresholds(read.csv(input$fTH$datapath))
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
