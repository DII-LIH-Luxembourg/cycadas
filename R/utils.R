availMarkers <- function(sm) {

  m <- colnames(df)
  return(setdiff(m, sm))
}
subsetdf <- function(markers) {

  return(df[1:20, 1:8])
}


kmeansTH <- function(df) {
  th <- data.frame(cell = colnames(df), threshold = 0.0, color = "blue", bi_mod = 0)

  for (m in th$cell) {

    # check for bi-modal distribution, if not color = red to indicate
    # stated in Pfister et al., 2013
    bi_mod_value <- bimodality_coefficient(df[, m])
    th[th$cell == m, "bi_mod"] <- round(bi_mod_value, 3)
    if (bi_mod_value < 0.555) {
      th[th$cell == m, "color"] <- "red"
    }
    
    set.seed(42)

    z <- Ckmeans.1d.dp(df[, m], 2)
    th[th$cell == m, "threshold"] <- round(ahist(z, style="midpoints", data=df[, m], plot=FALSE)$breaks[2:2], 3)
  }

  rownames(th) <- th$cell
  return(th)
}


normalize01 <- function(hm) {

  eDR <- as.matrix(hm)
  rng <- colQuantiles(eDR, probs = c(0.01, 0.99))
  expr01 <- t((t(eDR) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0
  expr01[expr01 > 1] <- 1

  return(as.data.frame(expr01))
}


filterHM <- function(DF,posList, negList, th) {

  # browser()

  if (length(posList) == 0 & length(negList) == 0) {
    return(as.data.frame(DF))

  } else if (length(posList) == 0 & length(negList) != 0) {
    negTH <- th$threshold[th$cell %in% negList]
    for (i in 1:length(negList)) {
      marker <- negList[i]
      DF <- DF[DF[marker] < negTH[i], ]
    }
    return(as.data.frame(DF))
  } else if (length(posList) != 0 & length(negList) == 0) {
    posTH <- th$threshold[th$cell %in% posList]
    for (i in 1:length(posList)) {
      marker <- posList[i]
      DF <- DF[DF[marker] >= posTH[i], ]
    }
    return(as.data.frame(DF))

  } else {

    # browser()

    posTH <- th$threshold[th$cell %in% posList]
    negTH <- th$threshold[th$cell %in% negList]
    # first reduce by the positive markers

    for (i in 1:length(posList)) {
      marker <- posList[i]
      DF <- DF[DF[marker] >= posTH[i], ]
      # DF <- DF[DF[marker] > th[marker,]$threshold, ]
    }
    # next reduce by the negative markers

    for (i in 1:length(negList)) {
      marker <- negList[i]
      DF <- DF[DF[marker] < negTH[i], ]
      # DF <- DF[DF[marker] < th[marker,]$threshold, ]
    }
    return(as.data.frame(DF))
  }
}


filterColor <- function(DF,hm) {
  # browser()
  tmp <- replicate(nrow(DF), "other clusters")
  tmp[as.numeric(rownames(hm))] <- "selected phenotype"

  return(unname(tmp))

}
filterAnnotation <- function(hm_tmp,at_tmp) {

  rn <- rownames(hm_tmp)
  sub_annotation <- at_tmp$data[rn,]
  sub_annotation <- as.data.frame(sub_annotation)
  rownames(sub_annotation) <- rn

  return(sub_annotation)
}
set_ptname <- function(ptname, df) {
  df$phenotype <- ptname
}


setPhenotypeName <- function(markers, s, ph_name) {
  if (s == "pos") {
    ph_name <- paste0(ph_name, paste0(markers, collapse = "+"), "+")
  } else if (s == "neg") {
    ph_name <- paste0(ph_name, paste0(markers, collapse = "-"), "-")
  }
  return(ph_name)
}

# a function which generates the dataframe for plotting the distribution
getMarkerDistDF <- function(marker, myScale){

  marker_expr <- NULL

  if (myScale == "1") {
    marker_expr <- as.data.frame(df01[, c(marker)])
    marker_expr <- cbind(marker_expr, rnorm(1:dim(marker_expr)[1]))

  } else if (myScale == "2") {
    tmp <- log10(sinh(df_global[, c(marker)]) * 5)
    marker_expr <- as.data.frame(tmp)
    marker_expr <- cbind(marker_expr, rnorm(1:dim(marker_expr)[1]))
  }

  return(marker_expr)
}


# a function to add a new row for nodes and edges
add_node <- function(graph,parent,name,pos_m,neg_m,color) {

  next_id <- max(graph$nodes$id) + 1
  parent_row <- graph$nodes %>% filter(label == parent)
  # node_level <- paste0(parent_row$id,".",2)
  graph$nodes <- graph$nodes %>% add_row(id = next_id,
                                         label = name,
                                         pm = pos_m,
                                         nm = neg_m,
                                         color = color)

  graph$edges <- graph$edges %>% add_row(from = next_id, to = parent_row$id)

  return(graph)
}

# define recursive function to delete a node and all its children
delete_child_nodes <- function(graph_data, node_id) {
  # find the children nodes
  # browser()
  children <- graph_data$edges$from[graph_data$edges$to == node_id]

  # recursively delete children nodes
  if (length(children) > 0) {
    for (child in children) {
      graph_data <- delete_child_nodes(graph_data, child)
    }
  }

  # delete the node and its edges from the graph data
  graph_data$nodes <- graph_data$nodes[graph_data$nodes$id != node_id, ]
  graph_data$edges <- graph_data$edges[!(graph_data$edges$from == node_id | graph_data$edges$to == node_id), ]

  return(graph_data)
}

all_my_children <- function(graph_data, node_id) {

  # we don't want the root node
  if (node_id > 1) {
    children <- graph_data$edges$from[graph_data$edges$to == node_id]
    i <- 1
    while (i <= length(children)) {
      children <- append(children, graph_data$edges$from[graph_data$edges$to == children[i]])
      i <- i + 1
    }
    return (children)
  } else {
    return (c())
  }



}

# define recursive function to delete a node and all its children
delete_leaf_node <- function(graph_data, node_id) {
  # find the children nodes
  # browser()
  children <- graph_data$edges$from[graph_data$edges$to == node_id]

  # recursively delete children nodes
  if (length(children) > 0) {

    showModal(modalDialog(
      title = "Select Leaf Node",
      "Only Leaf Nodes can be deleted or modified!",
      easyClose = TRUE,
      footer = NULL
    ))
    return(graph_data)

  } else {
    # browser()
    # delete the node and its edges from the graph data
    graph_data$nodes <- graph_data$nodes[graph_data$nodes$id != node_id, ]
    graph_data$edges <- graph_data$edges[!(graph_data$edges$from == node_id | graph_data$edges$to == node_id), ]

    return(graph_data)
  }

}



# deleteChildNodes <- function(graph, node) {
#   children <- successors(graph, node)  # Get the immediate children of the node
#   if (length(children) > 0) {  # If the node has children, recursively delete them
#     for (child in children) {
#       graph <- delete.edges(graph, c(node, child))  # Remove the edge between the node and its child
#       graph <- deleteChildNodes(graph, child)  # Recursively delete the child's children
#     }
#   }
#   return(graph)  # Return the modified graph
# }
#
#
#
# getUpstreamMarkers <- function(label, nodes, edges) {
#
#   node <- nodes %>% filter(label == label)
#   myid <- node$id
#
#   filterPosMarkers <- unlist(node$pm)
#   filterNegMarkers <- unlist(node$nm)
#
#   # collect all positive and negative markers
#   # upwards from parent
#   # go through edges until we reach one level below master node
#   while (myid > 1) {
#     print(myid)
#
#     edge <- edges %>% filter(from == myid)
#     next_id <- edge$to
#
#     myid <- next_id
#     next_node <- nodes %>% filter(id == next_id)
#     filterPosMarkers <- append(filterPosMarkers, unlist(next_node$pm))
#     filterNegMarkers <- c(filterNegMarkers, unlist(next_node$nm))
#
#   }
#   print(filterPosMarkers)
#   print(filterNegMarkers)
#
#   # remove the empty strings in the markers
#   filterPosMarkers <- filterPosMarkers[nzchar(filterPosMarkers)]
#   filterNegMarkers <- filterNegMarkers[nzchar(filterNegMarkers)]
#
#   return(list(filterPosMarkers,filterNegMarkers))
# }

# plotTree <- function(mygraph) {
#
#   # browser()
#
#   # before plotting the tree, update its properties
#   # based on the annotated DF
#   if (!is.null(mygraph)) {
#     for(i in 1:nrow(mygraph$nodes)) {
#       l <- mygraph$nodes$label[i]
#
#       # get the number of clusters with that label
#       nLabel <- sum(df01Tree$cell == l)
#       if (nLabel == 0) {
#         browser()
#         mygraph$nodes$color[i] <- "red"
#       }
#     }
#   }
#
#   if(input$treeSwitch == T) {
#     output$mynetworkid <- renderVisNetwork({
#
#       visNetwork(mygraph$nodes, mygraph$edges, width = "100%") %>%
#         visEdges(arrows = "from")
#     })
#
#   } else {
#     output$mynetworkid <- renderVisNetwork({
#
#       visNetwork(mygraph$nodes, mygraph$edges, width = "100%") %>%
#         visEdges(arrows = "from") %>%
#         visHierarchicalLayout()
#     })
#   }
#
# }
#











