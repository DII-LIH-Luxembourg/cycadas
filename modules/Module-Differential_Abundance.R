
Differential_Abundance_UI <- function(id) {
  ns <- NS(id)

  fluidRow(column(width = 6,
                 box(width = NULL,
                      fluidRow(
                        column(width = 5,
                               box(width = NULL,title = "Metadata preview",tableOutput(ns("md_table")))),
                        column(width = 5,
                               box(width = NULL,title = "Counts Table preview",tableOutput(ns("counts_table"))))  
                      ),
                  ),
                  box(width = NULL,
                      fluidRow(
                        column(width = 6,
                               box(width = NULL,selectInput(ns("correction_method"), "Select:",
                                                            choices = c("holm", "hochberg", "hommel", "bonferroni",
                                                                        "BH", "BY","fdr", "none")))),
                        column(width = 6,
                               box(width = NULL,title = "Do Analysis",actionButton(ns("doDA"), "Calculate")))
                        )
                      ),
                  box(width = NULL,
                      fluidRow(
                        column(width = 6,
                               box(width = NULL,title = "Export DA Result",downloadButton(ns("exportDA"), "Download"))),
                        column(width = 6,
                               box(width = NULL,title = "Export Proportion Table",downloadButton(ns("exportProp"), "Download"))),
                        )
                      ),
                  ),
           column(width = 6,
                  box(width = NULL,title = "DA Result",tableOutput(ns("DA_result_table")))
                  )
           )
  }

Differential_Abundance_Server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("Run Differential Abundance Tab")

      
      # Server - Differential Abundance Tab ----
      output$md_table <-
        renderTable(
          session$userData$vars$md[1:5, ]
        )
      
      output$counts_table <-
        renderTable(
          session$userData$vars$counts_table[1:5, 1:5]
        )
      
       
          # Do the differential abundance ----
    observeEvent(input$doDA, {
      req(session$userData$vars$counts_table, session$userData$vars$md)

      countsTable <- session$userData$vars$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df01Tree$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL

      props_table <- t(t(countsTable) / colSums(countsTable)) * 100

      # mm <- match(md$sample_id, colnames(props_table))
      mm <- match(colnames(props_table), session$userData$vars$md$sample_id)

      tmp_cond <- session$userData$vars$md$condition[mm]

      DA_df <- data.frame()
      for (i in 1:nrow(props_table)) {

        ##### Does the p.adjust.method work??? #####
        foo <- pairwise.wilcox.test(as.numeric(props_table[i,]), tmp_cond, p.adjust.method=input$correction_method)
        df <- subset(melt(foo$p.value), value!=0)
        df$cell <- rownames(countsTable)[i]
        DA_df <- rbind(DA_df, as.data.frame(df))

      }

      # browser()
      # change the name of nodes:
      # if node has children add "_remaining"
      my_list <- lapply(DA_df$cell, function(x) {
            # get the id
            nid <- session$userData$vars$graph$nodes$id[session$userData$vars$graph$nodes$label == x]
            # check if the id has children
            if(nid %in% session$userData$vars$graph$edges$to) {
              return (x <- paste0(x, "_remaining"))
            }
            else {
              return(x)
            }
          })

      DA_df$list <- unlist(my_list)
      
      colnames(DA_df) <- c("Cond1", "Cond2", "p-value", "Cell", "Naming")

      session$userData$vars$DA_result_table <- DA_df

      output$DA_result_table <-
        renderTable(
          session$userData$vars$DA_result_table
        )
         
      
       })
      
      

      

get_merged_prop_table <- function() {
        
        req(session$userData$vars$counts_table)
        
        countsTable <- session$userData$vars$counts_table
        ## aggregate the clusters by name:
        countsTable['cell'] <- df01Tree$cell
        # merge and aggregate by cell
        countsTable <- aggregate(. ~ cell, countsTable, sum)
        
        rownames(countsTable) <- countsTable$cell
        countsTable$cell <- NULL
        
        props_table <- t(t(countsTable) / colSums(countsTable)) * 100
        
        return(props_table)
        
      }
      
      
      
# Save the Merged Proportion Table ----
      output$exportProp <- downloadHandler(
        
        filename = function() {
          paste("Merged_Proportions_Table_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(get_merged_prop_table(), file)
        }
      )
      
      
# Save the DA Table ----
      output$exportDA <- downloadHandler(
        filename = function() {
          paste("DA_Table_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(session$userData$vars$DA_result_table, file)
        }
      )
      
    })}