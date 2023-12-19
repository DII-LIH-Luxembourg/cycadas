TreeAnnotation_UI <- function(id) {
  ns <- NS(id)
      ))
  }

TreeAnnotation_Server <- function(id,reactVals,df01,tab) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
#       # req(exists("df01"))
#  
#       # TEST5 <<- reactVals
#       print("Run Tree Annotation")
#       
# # Function: Update cluster labels ----
#       updateClusterLabels <- function(mydf) {
#         
#         mysum <- sum(cell_freq[rownames(mydf),]$clustering_prop)
#         
#         output$progressBox <- renderText({
#           paste0(mysum, "%")
#         })
#         output$progressBox2 <- renderText({
#           paste0(dim(mydf)[1], "Cluster")
#         })
#       }
#       
# # Upload the marker expression file ----
#       # observeEvent(c(input$fMarkerExpr, input$cluster_freq), {
#         
#         session$userData$vars$graph <- initTree()
#         # browser()
#         
#         # req(input$fMarkerExpr)
#         # req(input$cluster_freq)
#         
#         # df <- read.csv(input$fMarkerExpr$datapath)
#         # df_global <<- df
#         
#         annotationlist <<- list("Unassigned")
#         
#         cell_freq <<- read.csv(input$cluster_freq$datapath)
#         
#         labels_row <-
#           paste0(rownames(df), " (", cell_freq$clustering_prop , "%)")
#         
#         marker_names <- rownames(df)
#         
#         set.seed(1234)
#         
#         my_umap <- umap(df)
#         dr_umap <<- data.frame(
#           u1 = my_umap$layout[, 1],
#           u2 = my_umap$layout[, 2],
#           my_umap$data,
#           cluster_number = 1:length(my_umap$layout[, 1]),
#           check.names = FALSE
#         )
#         
#         df01 <<- df %>% normalize01()
#         myDF <<- df %>% normalize01()
#         
#         df01Tree <<- df %>% normalize01()
#         allMarkers <<- colnames(df)
#         df01Tree$cell <<- "Unassigned"
#         
#         selectedMarkers <<- colnames(df)
#         posPickerList <<- colnames(df)
        
        
        #### IS THIS DOING ANYTHING
        # ## load only at start to fill the picker list
        # updatePickerInput(
        #   session,
        #   inputId = "myPickerPos",
        #   label = "Select Positive Markers",
        #   choices = colnames(df01),
        #   # choices = NULL,
        #   options = list(
        #     `actions-box` = TRUE,
        #     size = 10,
        #     `selected-text-format` = "count > 3"
        #   )
        # )
        # updatePickerInput(
        #   session,
        #   inputId = "myPickerNeg",
        #   label = "Select Negative Markers",
        #   choices = colnames(df01),
        #   options = list(
        #     `actions-box` = TRUE,
        #     size = 10,
        #     `selected-text-format` = "count > 3"
        #   )
        # )
        

        
      #   session$userData$vars$th <- kmeansTH(df01)
      #
      #   annotaionDF <<- data.frame("cell" = "unassigned",
      #                              clusterSize = cell_freq$clustering_prop)
      #   at <<- reactiveValues(data = annotaionDF, dr_umap = dr_umap)
      #
      # # })
      # # Observe MenutItems ----
      # observeEvent(input$tabs, {
      #
      #   if (exists("df01")) {
      #     if(input$tabs=="thresholds") { ## Thresholds ----
      #
      #       if (is.null(session$userData$vars$th) & !is.null(df01)) {
      #         session$userData$vars$th <- kmeansTH(df01)
      #       }
      #     }
      #     else if(input$tabs=="treeannotation") { # Tree annotation ----
      #
      #       updatePickerInput(
      #         session,
      #         inputId = "parentPicker",
      #         choices = annotationlist
      #       )
      #       updatePickerInput(
      #         session,
      #         inputId = "updateNodePicker",
      #         choices = annotationlist,
      #         selected = ""
      #       )
      #       plotTree()
      #     }
      #   }
      # })










# 
# 
#       # Observe node update picker ----
#       observeEvent(input$updateNodePicker, {
# 
#         browser()
# 
#         node <- session$userData$vars$graph$nodes %>% filter(label == input$updateNodePicker)
#         myid <- node$id
# 
#         filterPosMarkers <- unlist(node$pm)
#         filterNegMarkers <- unlist(node$nm)
# 
#         updatePickerInput(
#           session,
#           inputId = "updatePickerPos",
#           label = "Select Positive Markers",
#           choices = colnames(df01),
#           selected = filterPosMarkers,
#           options = list(
#             `actions-box` = TRUE,
#             size = 10,
#             `selected-text-format` = "count > 3"
#           )
#         )
#         updatePickerInput(
#           session,
#           inputId = "updatePickerNeg",
#           label = "Select Negative Markers",
#           choices = colnames(df01),
#           selected = filterNegMarkers,
#           options = list(
#             `actions-box` = TRUE,
#             size = 10,
#             `selected-text-format` = "count > 3"
#           )
#         )
#         updateTextInput(
#           session,
#           inputId = "renameNode",
#           label = NULL,
#           value = node$label)
#       })
# 
# 
# 
# 
# 
# 
#       importNodes <- function(id, label, pm, nm, parent_id, color) {
# 
#         if(id == parent_id) {return()}
# 
#         parent_row <- session$userData$vars$graph$nodes %>% filter(id == parent_id)
#         parent_name <- parent_row$label
# 
#         add_node(session$userData$vars$graph, parent_name, label, color)
#       }



      
    })}




# 

