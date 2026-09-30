# RemoteSOM import ------------------------------------------------------------
# Reads the feature names, cluster counts and median expression JSON files.

mod_import_remotesom_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$p(class = "text-body-secondary small",
           "Feature names, cluster counts and median expression exported by RemoteSOM, as JSON files."),
    tags$form(id = ns("upload_form"),
              fileInput(ns("features"), "Features Names", accept = ".json"),
              fileInput(ns("counts"), "Cluster Counts", accept = ".json"),
              fileInput(ns("medians"), "Median Expression", accept = ".json")
    ),
    uiOutput(ns("status")),
    actionButton(ns("btnImport"), "Import", icon = icon("file-import"), class = "btn-primary mt-3")
  )
}

mod_import_remotesom_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    files <- file_tracker(input, c("features", "counts", "medians"))

    observe(shinyjs::toggleState("btnImport", condition = files$ready()))

    output$status <- renderUI({
      file_status_ui(files$present(),
                     c("Feature Names", "Cluster Counts", "Median Expression"),
                     "All files present, ready to import.",
                     "Select all three files to enable Import.")
    })

    observeEvent(input$btnImport, {
      req(files$ready())

      tryCatch({
        load_dataset(state, read_remotesom_json(input$features$datapath,
                                                input$counts$datapath,
                                                input$medians$datapath))
        showNotification("RemoteSOM data imported.", type = "message")
      }, error = function(e) {
        showNotification(paste("RemoteSOM import failed:", conditionMessage(e)),
                         type = "error", duration = 8)
      })
      files$clear()
    })

    observeEvent(state$reset_version, {
      shinyjs::reset("upload_form")
      files$clear()
    }, ignoreInit = TRUE)
  })
}
