TreeAnnotation_annotationdownload_UI <- function(id) {
  ns <- NS(id)
  downloadButton(ns("exportAnnotationBtn"), "Export Annotation")
    
  }

TreeAnnotation_annotationdownload_Server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {

# Export Annotation Tree ----
      output$exportAnnotationBtn <- downloadHandler(

        filename = function(){
          paste("cycadas_annotation_data_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file) {
          temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
          dir.create(temp_directory)

          export_df_nodes <- session$userData$vars$graph$nodes
          # Use apply() to concatenate the pm column for each row of the nodes dataframe
          pm_concatenated <- apply(export_df_nodes, 1, function(row) {

            pm_vector <- unlist(row["pm"])
            pm_vector <- paste(pm_vector, collapse = "|")
          })

          nm_concatenated <- apply(export_df_nodes, 1, function(row) {
            nm_vector <- unlist(row["nm"])
            nm_vector <- paste(nm_vector, collapse = "|")
          })

          # Add the concatenated pm column as a new column to the nodes dataframe
          export_df_nodes$pm <- pm_concatenated
          export_df_nodes$nm <- nm_concatenated

          download_list <- list(annTable = df01Tree,
                                nodesTable = export_df_nodes,
                                edgesTable = session$userData$vars$graph$edges)
          # download_list <- list(annTable = df01Tree["cell"],
          #                       nodesTable = export_df_nodes,
          #                       edgesTable = session$userData$vars$graph$edges)

          download_list %>%
            imap(function(x,y){
              if(!is.null(x)){
                file_name <- glue("{y}_data.csv")
                readr::write_csv(x, file.path(temp_directory, file_name))
              }
            })

          zip::zip(
            zipfile = file,
            files = dir(temp_directory),
            root = temp_directory
          )
        },
        contentType = "application/zip"
      )

      
    })}




# 

