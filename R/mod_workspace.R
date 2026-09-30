# Workspace -------------------------------------------------------------------
# Save / load / clear the workspace, load the demo data and show what is loaded.

mod_workspace_ui <- function(id) {
  ns <- NS(id)
  tagList(
    card(
      card_header("Workspace"),
      tags$p(class = "small text-body-secondary",
             "A workspace file holds the data, thresholds and annotation tree of this session."),
      tags$div(class = "d-flex gap-2 align-items-start flex-wrap",
               downloadButton(ns("btnSaveWorkspace"), "Save workspace", class = "btn-primary"),
               actionButton(ns("btnClearWorkspace"), "Clear", icon = icon("eraser"),
                            class = "btn-outline-danger")),
      tags$div(id = ns("load_form"), class = "mt-3",
               fileInput(ns("btnLoadWorkspace"), "Load workspace",
                         placeholder = "cycadas-workspace_*.rds",
                         multiple = FALSE, accept = c(".rds", ".RDS")))
    ),
    card(card_header("Loaded data"), uiOutput(ns("status")))
  )
}

# Demo data buttons, shown on the import card
mod_workspace_demo_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$p(class = "text-body-secondary small",
           "1,600 clusters and 27 markers from a PBMC mass cytometry dataset."),
    tags$div(class = "d-flex gap-2",
             actionButton(ns("btnLoadDemoData"), "Load unannotated", class = "btn-outline-primary"),
             actionButton(ns("btnLoadAnnoData"), "Load annotated", class = "btn-primary")),
    tags$p(class = "text-body-secondary small mt-2 mb-0",
           "The annotated version includes thresholds, an annotation tree, metadata and counts.")
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
      dims <- function(x) if (has_df(x)) sprintf("%d \u00d7 %d", nrow(x), ncol(x)) else "empty"
      tags$ul(
        class = "status-list",
        status_item(has_df(state$expr), "Cluster expression", dims(state$expr)),
        status_item(has_val(state$cell_freq), "Cluster frequencies", dims(state$cell_freq)),
        status_item(has_df(state$th), "Thresholds", dims(state$th)),
        status_item(!is.null(state$graph), "Annotation tree",
                    if (!is.null(state$graph)) sprintf("%d phenotypes", nrow(state$graph$nodes)) else "empty"),
        status_item(has_df(state$md), "Sample metadata", dims(state$md)),
        status_item(has_df(state$counts), "Cluster counts", dims(state$counts)),
        status_item(has_val(state$sce), "CATALYST object",
                    if (has_val(state$sce)) class(state$sce)[1] else "empty")
      )
    })
  })
}
