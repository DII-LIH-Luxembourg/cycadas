
DeleteNode_UI <- function(id) {
  ns <- NS(id)
    }

DeleteNode_Server <- function(id,reactVals,filter) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("run Delete Node")
      # browser()
      
      # node <- session$userData$vars$graph$nodes %>% filter(label == input$updateNodePicker)
      node <- session$userData$vars$graph$nodes[session$userData$vars$graph$nodes$label == filter, ]
      
      # check and make sure that this node is a leaf node
      # browser()
      if(TRUE %in% (session$userData$vars$graph$edges$to == node$id)) {
        
        showModal(modalDialog(
          title = "Cannot delete Node",
          "The selected Node is not a Leaf Node!",
          easyClose = TRUE,
          footer = NULL
        ))
        
        return()
        
      }
      else {
        
        # before deleting the node, we must assign the cluster annotation
        # name for that node back to its parent
        parent_id <- session$userData$vars$graph$edges$to[session$userData$vars$graph$edges$from == node$id]
        parent_label <- session$userData$vars$graph$nodes$label[session$userData$vars$graph$nodes$id == parent_id]
        
        # In case there is a empty node with no the filtered clusters
        # we do not assign any labeling
        if( dim(df01Tree[df01Tree$cell == node$label,])[1] >0 ) {
          
          df01Tree[df01Tree$cell == node$label,]$cell <<- parent_label
          
        }
        
        annotationlist[annotationlist == node$label] <<- NULL
        
        updatePickerInput(
          session,
          inputId = "parentPicker",
          selected = NULL,
          choices = annotationlist
        )
        updatePickerInput(
          session,
          inputId = "updateNodePicker",
          selected = NULL,
          choices = annotationlist
        )
        updateTextInput(
          session,
          inputId = "renameNode",
          label = NULL,
          value = ""
        )
        updateTextInput(
          session,
          inputId = "newNode",
          label = NULL,
          value = ""
        )
        
        session$userData$vars$graph <- delete_leaf_node(session$userData$vars$graph, node$id)
        # plotTree()
      }
    })}
