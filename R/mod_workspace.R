# Workspace -------------------------------------------------------------------
# Save / load / clear the workspace, load the demo data and show what is loaded.

mod_workspace_ui <- function(id) {
  ns <- NS(id)
  tagList(
    settings_box(
      "Workspace", status = "success",
      splitLayout(cellWidths = c("50%", "50%"),
                  downloadButton(ns("btnSaveWorkspace"), "Save Workspace"),
                  tags$div(id = ns("load_form"),
                           fileInput(ns("btnLoadWorkspace"), label = NULL,
                                     placeholder = "Choose Workspace File",
                                     multiple = FALSE,
                                     accept = c(".rds", ".RDS")))
      ),
      tags$hr(),
      actionButton(ns("btnClearWorkspace"), "Clear Workspace")
    ),
    settings_box(
      "Load Demo Data", status = "success",
      splitLayout(cellWidths = c("50%", "50%"),
                  actionButton(ns("btnLoadDemoData"), "Unannotated"),
                  actionButton(ns("btnLoadAnnoData"), "Annotated")
      )
    ),
    box(title = "Workspace status", width = NULL, solidHeader = TRUE, status = "primary",
        uiOutput(ns("status")))
  )
}

mod_workspace_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    # Demo data ----
    observeEvent(input$btnLoadDemoData, {
      load_dataset(state, demo_dataset())
    })

    observeEvent(input$btnLoadAnnoData, {
      th <- prepare_thresholds(getExportedValue("cycadas", "marker_th_demo_data"))
      load_dataset(state, demo_dataset(), th = th)

      state$md <- drop_index_column(getExportedValue("cycadas", "meta_demo_data"))
      state$counts <- drop_index_column(getExportedValue("cycadas", "cluster_counts_demoData"))
      state$graph <- getGraphFromLoad(getExportedValue("cycadas", "nodes_demo_data"),
                                      getExportedValue("cycadas", "edges_demo_data"))
      rebuild_annotation(state)
    })

    # Save / load / clear ----
    output$btnSaveWorkspace <- downloadHandler(
      filename = function() {
        paste0("cycadas-workspace_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".rds")
      },
      content = function(file) {
        save_workspace(file, workspace_from_state(state))
      }
    )

    observeEvent(input$btnLoadWorkspace, {
      req(input$btnLoadWorkspace$datapath)
      restore_workspace(state, load_workspace(input$btnLoadWorkspace$datapath))
      showNotification("Workspace loaded.", type = "message")
    })

    observeEvent(input$btnClearWorkspace, {
      reset_app_state(state)
      shinyjs::reset("load_form")
      showNotification("Workspace cleared.", type = "message")
    })

    # Status ----
    output$status <- renderUI({
      dims <- function(x) if (has_df(x)) sprintf("%d×%d", nrow(x), ncol(x)) else "-"
      df <- data.frame(
        Item    = c("Median expression", "Metadata", "Thresholds", "Cell Frequencies",
                    "Cluster Counts table", "CATALYST object"),
        Present = c(has_df(state$expr), has_df(state$md), has_df(state$th),
                    has_val(state$cell_freq), has_df(state$counts), has_val(state$sce)),
        Details = c(dims(state$expr), dims(state$md), dims(state$th), dims(state$cell_freq),
                    dims(state$counts),
                    if (has_val(state$sce)) class(state$sce)[1] else "-"),
        stringsAsFactors = FALSE
      )

      items <- lapply(seq_len(nrow(df)), function(i) {
        ok <- isTRUE(df$Present[i])
        tags$li(
          style = "margin: 6px 0;",
          strong(df$Item[i]), " — ",
          tags$span(class = paste("badge", if (ok) "bg-green" else "bg-red"),
                    if (ok) "Present" else "Empty"),
          " ",
          tags$i(class = paste("fa", if (ok) "check-circle" else "times-circle"),
                 style = "margin-left:6px;"),
          tags$span(style = "margin-left:10px;color:#666;", df$Details[i])
        )
      })
      tags$ul(class = "list-unstyled", items)
    })
  })
}
