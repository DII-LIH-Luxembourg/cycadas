# CATALYST import -------------------------------------------------------------
# Reads a SingleCellExperiment saved as RDS. Only offered when CATALYST is
# installed. The meta-cluster level is chosen on the Tree-Annotation tab.

mod_import_catalyst_ui <- function(id) {
  ns <- NS(id)
  settings_box("Catalyst import", uiOutput(ns("upload")))
}

mod_import_catalyst_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    output$upload <- renderUI({
      state$reset_version  # re-render to clear the file input
      if (catalyst_available()) {
        fileInput(session$ns("sce"), "Upload RDS",
                  placeholder = "Choose RDS File",
                  multiple = FALSE,
                  accept = c(".rds"))
      } else {
        p("The CATALYST package is not available.")
      }
    })

    observeEvent(input$sce, {
      req(has_file(input$sce))

      sce <- readRDS(input$sce$datapath)
      load_dataset(state, read_catalyst(sce))
      state$sce <- sce
      state$meta_level <- colnames(sce@metadata$cluster_codes)[1]
    })
  })
}
