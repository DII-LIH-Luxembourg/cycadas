# GigaSOM / FlowSOM import ----------------------------------------------------
# Reads the marker expression and cluster frequency CSVs into the shared state.

mod_import_gigasom_ui <- function(id) {
  ns <- NS(id)
  settings_box(
    "GigaSOM / FlowSOM import",
    tags$form(id = ns("upload_form"),
              fileInput(ns("fMarkerExpr"), "Upload Marker Expressions", accept = c("text/csv", ".csv")),
              fileInput(ns("cluster_freq"), "Upload Cluster Frequencies", accept = c("text/csv", ".csv"))
    ),
    tags$hr(),
    uiOutput(ns("status")),
    actionButton(ns("btnImport"), "Import", class = "btn btn-success")
  )
}

mod_import_gigasom_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    files <- file_tracker(input, c("fMarkerExpr", "cluster_freq"))

    observe(shinyjs::toggleState("btnImport", condition = files$ready()))

    output$status <- renderUI({
      file_status_ui(files$present(),
                     c("Marker Expressions (CSV)", "Cluster Frequencies (CSV)"),
                     "Both files present — ready to import.",
                     "Please select both files to enable Import.")
    })

    observeEvent(input$btnImport, {
      req(files$ready())

      ok <- tryCatch({
        load_dataset(state, read_gigasom_csv(input$fMarkerExpr$datapath,
                                             input$cluster_freq$datapath))
        TRUE
      }, error = function(e) FALSE)

      if (ok) {
        showNotification("GigaSOM / FlowSOM data imported.", type = "message")
      } else {
        showNotification("Import failed. Please check your CSVs.", type = "error", duration = 8)
        shinyjs::reset("upload_form")
      }
      # the button stays disabled until new files are chosen
      files$clear()
    })

    observeEvent(state$reset_version, {
      shinyjs::reset("upload_form")
      files$clear()
    }, ignoreInit = TRUE)
  })
}
