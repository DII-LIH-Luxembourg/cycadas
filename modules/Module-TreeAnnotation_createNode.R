TreeAnnotation_createNode_UI <- function(id) {
  ns <- NS(id)
    }

TreeAnnotation_createNode_Server <- function(id,parent,newNodeName) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("run TreeAnnotation_createNode")
      # browser()
      if (is.null(session$userData$vars$treePickerNeg) & is.null(session$userData$vars$treePickerPos)) {
        showModal(modalDialog(title = "No Marker Selection","Select positive and / or negative markers!",easyClose = TRUE))
        } else if (newNodeName == "") {
          showModal(modalDialog(title = "Phenotype Name","Set a Name for this Phenotype!",easyClose = TRUE))
          } else if (newNodeName %in% session$userData$vars$graph$nodes$label) {
            showModal(modalDialog(title = "Naming Error!","The Name for this Phenotype is already taken!",easyClose = TRUE))
            } else {
              name <- newNodeName
              # receive the parent settings, resp. parent hm
              # filter hm by parent cell name
              tmp <- df01Tree[df01Tree$cell == parent, ]
              tmp <- filterHM(tmp,session$userData$vars$treePickerPos, session$userData$vars$treePickerNeg, session$userData$vars$th)
          # Make sure the selection is not empty!
          if(dim(tmp)[1] == 0) {
            showModal(modalDialog(title = "No Result","The Selection for this Phenotype is empty!",easyClose = TRUE))
            } else {
            session$userData$vars$annotationlist <- append(session$userData$vars$annotationlist, name)
            posmarker <- list(session$userData$vars$treePickerPos)
            negmarker <- list(session$userData$vars$treePickerNeg)
            session$userData$vars$graph <- add_node(session$userData$vars$graph,parent,name,posmarker,negmarker,"blue")
            df01Tree[rownames(tmp),]$cell <<- name
            #updateClusterLabels(tmp)
        } 
      print("finish create node")
      }
    })}
