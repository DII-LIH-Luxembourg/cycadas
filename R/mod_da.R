# Differential Abundance tab --------------------------------------------------
# Pairwise Wilcoxon tests of phenotype proportions between conditions.

mod_da_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 320, open = "always", title = "Differential abundance",
      tags$p(class = "small text-body-secondary",
             "Pairwise Wilcoxon tests of phenotype proportions between conditions."),
      uiOutput(ns("inputs")),
      selectInput(ns("correction_method"), "P-value adjustment",
                  choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
      actionButton(ns("doDA"), "Run tests", icon = icon("play"), class = "btn-primary w-100"),
      tags$hr(),
      tags$div(class = "section-label", "Export"),
      downloadButton(ns("exportDA"), "Test results (CSV)", class = "btn-outline-secondary btn-sm w-100 mb-2"),
      downloadButton(ns("exportProp"), "Proportion table (CSV)", class = "btn-outline-secondary btn-sm w-100")
    ),
    card(full_screen = TRUE, card_header("Results"), DTOutput(ns("DA_result_table"))),
    accordion(
      open = FALSE,
      accordion_panel("Metadata preview", tableOutput(ns("md_table"))),
      accordion_panel("Counts preview", tableOutput(ns("counts_table")))
    )
  )
}

mod_da_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    result <- reactiveVal(NULL)
    observeEvent(state$dataset_version, result(NULL), ignoreInit = TRUE)

    output$md_table <- renderTable(state$md[1:5, ])
    output$counts_table <- renderTable(state$counts[1:5, 1:5])
    output$inputs <- renderUI({
      tags$ul(class = "status-list mb-3",
              status_item(!is.null(state$expr), "Annotated clusters"),
              status_item(has_df(state$md), "Sample metadata"),
              status_item(has_df(state$counts), "Cluster counts"))
    })

    observe(shinyjs::toggleState("doDA", condition = has_df(state$md) && has_df(state$counts)))

    output$DA_result_table <- DT::renderDT({
      req(result())
      DT::datatable(result(), rownames = FALSE, filter = "top",
                    options = list(pageLength = 25, scrollX = TRUE)) %>%
        DT::formatSignif("p-value", digits = 3)
    })

    observeEvent(input$doDA, {
      req(state$counts, state$md)

      tryCatch(
        withCallingHandlers({
          result(run_da(state$counts, state$md, state$expr$cell, state$graph,
                        input$correction_method))
          showNotification("Differential abundance testing completed.", type = "message")
        }, warning = function(w) {
          showNotification(conditionMessage(w), type = "warning")
          invokeRestart("muffleWarning")
        }),
        error = function(e) showNotification(conditionMessage(e), type = "error")
      )
    })

    output$exportDA <- downloadHandler(
      filename = function() {
        paste("DA_Table_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(result(), file)
      }
    )

    output$exportProp <- downloadHandler(
      filename = function() {
        paste("Merged_Proportions_Table_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        req(state$counts)
        write.csv(merged_prop_table(state$counts, state$expr$cell, state$graph), file)
      }
    )
  })
}
