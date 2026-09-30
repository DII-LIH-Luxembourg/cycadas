# Differential Abundance tab --------------------------------------------------
# Pairwise Wilcoxon tests of phenotype proportions between conditions.

mod_da_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(width = 6,
           box(width = NULL,
               fluidRow(
                 column(width = 5,
                        box(width = NULL, title = "Metadata preview", tableOutput(ns("md_table")))),
                 column(width = 5,
                        box(width = NULL, title = "Counts Table preview", tableOutput(ns("counts_table"))))
               )
           ),
           box(width = NULL,
               fluidRow(
                 column(width = 6,
                        box(width = NULL,
                            selectInput(ns("correction_method"), "Select:",
                                        choices = c("holm", "hochberg", "hommel", "bonferroni",
                                                    "BH", "BY", "fdr", "none")))),
                 column(width = 6,
                        box(width = NULL, title = "Do Analysis", actionButton(ns("doDA"), "Calculate")))
               )
           ),
           box(width = NULL,
               fluidRow(
                 column(width = 6,
                        box(width = NULL, title = "Export DA Result",
                            downloadButton(ns("exportDA"), "Download"))),
                 column(width = 6,
                        box(width = NULL, title = "Export Proportion Table",
                            downloadButton(ns("exportProp"), "Download")))
               )
           )
    ),
    column(width = 6,
           box(width = NULL, title = "DA Result", tableOutput(ns("DA_result_table"))))
  )
}

mod_da_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    result <- reactiveVal(NULL)
    observeEvent(state$dataset_version, result(NULL), ignoreInit = TRUE)

    output$md_table <- renderTable(state$md[1:5, ])
    output$counts_table <- renderTable(state$counts[1:5, 1:5])
    output$DA_result_table <- renderTable(result())

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
