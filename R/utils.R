
# Estimate the threshold values -----------------------------------------------
kmeansTH <- function(df, th_mode="km") {
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

updateTH <- function(df, th, th_mode) {
  
  set.seed(42)
  
  if (th_mode == "km") {
    
    for (m in th$cell) {
      set.seed(42)
      z <- Ckmeans.1d.dp(df[, m], 2)
      th[th$cell == m, "threshold"] <- round(ahist(z, style="midpoints", data=df[, m], plot=FALSE)$breaks[2:2], 3)
    }
  }
  else if (th_mode == "gmm_mid") {
    
    # browser()
    for (m in th$cell) {
      set.seed(42)
      fit <- normalmixEM(df[, m], k = 2)
      # Get component means
      means <- fit$mu
      # Define high expression threshold (e.g., mean of high component + standard deviation)
      # th[th$cell == m, "threshold"] <- round(means[2] - fit$sigma[2], 3)
      th[th$cell == m, "threshold"] <- round((means[2] - means[1]) / 2, 3)
    }
  }
    else if (th_mode == "gmm_high") {
      
      # browser()
      for (m in th$cell) {
        set.seed(42)
        fit <- normalmixEM(df[, m], k = 2)
        # Get component means
        means <- fit$mu
        # Define high expression threshold (e.g., mean of high component + standard deviation)
        th[th$cell == m, "threshold"] <- round(means[2] - fit$sigma[2], 3)
        # th[th$cell == m, "threshold"] <- round((means[2] - means[1]) / 2, 3)
      }
    }
    else if (th_mode == "gmm_low") {
      
      # browser()
      for (m in th$cell) {
        set.seed(42)
        fit <- normalmixEM(df[, m], k = 2)
        # Get component means
        means <- fit$mu
        # Define high expression threshold (e.g., mean of high component + standard deviation)
        th[th$cell == m, "threshold"] <- round(means[1] + fit$sigma[1], 3)
        # th[th$cell == m, "threshold"] <- round((means[2] - means[1]) / 2, 3)
      }
    }
    else if (th_mode == "mclust") {
      
      # browser()
      for (m in th$cell) {
        set.seed(42)
        
        # Fit the GMM model using Mclust
        model <- Mclust(df[, m], G = 2)

        max_clus1 <- max(model$data[model$classification==1])
        min_clus2 <- min(model$data[model$classification==2])

        th[th$cell == m, "threshold"] <- mean(max_clus1, min_clus2) 

      }
    }
    else if (th_mode == "kde") {
      
      # browser()
      for (m in th$cell) {
        set.seed(42)
        
        dens <- density(asinh(df[, m]))
        
        # Fit Gaussian Mixture Model (GMM) to the peaks
        peaks <- data.frame(x = dens$x, y = dens$y)
        fit <- Mclust(peaks, G = 2) # Assuming 2 clusters
        
        # Get means of the Gaussian components (cluster centers)
        cluster_means <- fit$parameters$mean
        
        # Sort means to identify the two peaks
        sorted_means <- sort(cluster_means)
        
        # Calculate threshold as the midpoint between the two cluster means
        threshold <- mean(sorted_means)
        
        # th[th$cell == m, "threshold"] <- quantile(dens$x, probs = 0.75)
        th[th$cell == m, "threshold"] <- threshold
        
      }
    }
  return(th)
}


# Scale the expression values between 0 and 1 ---------------------------------
normalize01 <- function(hm) {

  eDR <- as.matrix(hm)
  rng <- colQuantiles(eDR, probs = c(0.01, 0.99))
  expr01 <- t((t(eDR) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0
  expr01[expr01 > 1] <- 1

  return(as.data.frame(expr01))
}

# Filter the Heatmap ----------------------------------------------------------
filterHM <- function(DF,posList, negList, th) {

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

    posTH <- th$threshold[th$cell %in% posList]
    negTH <- th$threshold[th$cell %in% negList]
    
    # first reduce by the positive markers
    for (i in 1:length(posList)) {
      marker <- posList[i]
      DF <- DF[DF[marker] >= posTH[i], ]
    }
    # next reduce by the negative markers
    for (i in 1:length(negList)) {
      marker <- negList[i]
      DF <- DF[DF[marker] < negTH[i], ]
    }
    return(as.data.frame(DF))
  }
}

# Filter the color for UMAP ---------------------------------------------------
filterColor <- function(DF,hm) {
  # browser()
  tmp <- replicate(nrow(DF), "other clusters")
  tmp[as.numeric(rownames(hm))] <- "selected phenotype"

  return(unname(tmp))

}

# TODO: check if needed!
filterAnnotation <- function(hm_tmp,at_tmp) {

  rn <- rownames(hm_tmp)
  sub_annotation <- at_tmp$data[rn,]
  sub_annotation <- as.data.frame(sub_annotation)
  rownames(sub_annotation) <- rn

  return(sub_annotation)
}

# TODO: check if needed!
set_ptname <- function(ptname, df) {
  df$phenotype <- ptname
}

# Set a phonetype name from marker selection ----------------------------------
# Currently not used!
setPhenotypeName <- function(markers, s, ph_name) {
  if (s == "pos") {
    ph_name <- paste0(ph_name, paste0(markers, collapse = "+"), "+")
  } else if (s == "neg") {
    ph_name <- paste0(ph_name, paste0(markers, collapse = "-"), "-")
  }
  return(ph_name)
}

# a function to add a new row for nodes and edges -----------------------------
add_node <- function(graph,parent,name,pos_m,neg_m,color) {

  next_id <- max(graph$nodes$id) + 1
  parent_row <- graph$nodes %>% dplyr::filter(label == parent)
  # node_level <- paste0(parent_row$id,".",2)
  graph$nodes <- graph$nodes %>% add_row(id = next_id,
                                         label = name,
                                         pm = pos_m,
                                         nm = neg_m,
                                         color = color)

  graph$edges <- graph$edges %>% add_row(from = next_id, to = parent_row$id)

  return(graph)
}

# define recursive function to delete a node and all its children -------------
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

# Get all Children from a node selection --------------------------------------
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

# Initialize the Tree at start ------------------------------------------------
initTree <- function() {
  
  # create initial master node of all Unassigned clusters
  nodes <- tibble(id = 1,
                  label = "Unassigned",
                  pm = list(""),
                  nm = list(""),
                  color = "blue"
  )
  
  edges <- data.frame(from = c(1), to = c(1))
  
  return (list(nodes = nodes, edges = edges))
}

# Build Graph from nodes and edges --------------------------------------------
getGraphFromLoad <- function(df_nodes, df_edges) {
  
  df_nodes$pm[is.na(df_nodes$pm)] <- ""
  df_nodes$pm <- strsplit(df_nodes$pm, "\\|")
  
  df_nodes$nm[is.na(df_nodes$nm)] <- ""
  df_nodes$nm <- strsplit(df_nodes$nm, "\\|")
  
  return (list(nodes = df_nodes, edges = df_edges))
}


# Build UMAP and return as DF -------------------------------------------------
buildUMAP <- function(df_expr) {
  
  my_umap <- umap(df_expr)
  
  df_umap <- data.frame(
    u1 = my_umap$layout[, 1],
    u2 = my_umap$layout[, 2],
    my_umap$data,
    cluster_number = 1:length(my_umap$layout[, 1]),
    check.names = FALSE
  )
  
  return(df_umap)
}

# define recursive function to delete a node and all its children -------------
delete_leaf_node <- function(graph_data, node_id) {
  # find the children nodes
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
    # delete the node and its edges from the graph data
    graph_data$nodes <- graph_data$nodes[graph_data$nodes$id != node_id, ]
    graph_data$edges <- graph_data$edges[!(graph_data$edges$from == node_id | graph_data$edges$to == node_id), ]

    return(graph_data)
  }

}

# Rebuild Tree ----------------------------------------------------------------
# rebuilt the annotation tree if threshold values are changed, or
# if a tree is loaded
rebuiltTree <- function(graph, df_expr, th) {

  df_expr$cell <- "Unassigned"
  
  # for all rows in annotation table
  # get all parents for a row:
  if (nrow(graph$nodes) > 1) {
    graph$nodes$to <- graph$edges$to
    
    # Iterate through the dataframe row by row
    for (i in 1:nrow(graph$nodes)) {
      
      posMarker <- list()
      negMarker <- list()
      
      nodeID <- graph$nodes$id[i]
      parentID <- graph$nodes$to[i]
      
      # now for that nodeID, get all parents and collect
      # the pm and nm markers
      posMarker <- c(posMarker, unlist(graph$nodes$pm[graph$nodes$id == nodeID]))
      negMarker <- c(negMarker, unlist(graph$nodes$nm[graph$nodes$id == nodeID]))
      
      # now we have all positive and negative marker for that type collected
      # and we can start filtering the df
      # remove the empty strings in the markers
      posMarker <- posMarker[nzchar(posMarker)]
      negMarker <- negMarker[nzchar(negMarker)]
      
      # receive the parent settings, resp. parent hm
      # filter hm by parent cell name
      parent_label <- graph$nodes$label[graph$nodes$id == parentID]
      
      tmp_parent <- df_expr[df_expr$cell == parent_label, ]
      tmp <- filterHM(tmp_parent[, lineage_marker],unique(unlist(posMarker)), unique(unlist(negMarker)), th)
      
      df_expr[rownames(tmp), 'cell'] <- graph$node$label[i]
    }
  }

  return(df_expr$cell)
}

# Function for the creation of the oveerall expression DF
createExpressionDF <- function(df_expr, cell_freq) {
  
  df01_expr <- df_expr %>% normalize01()
  colnames(df01_expr) <- lineage_marker
  colnames(df_expr) <- lineage_marker_raw
  
  df_expr <- cbind(df_expr, df01_expr)
  
  ## Add frequencies and annotation
  df_expr$freq <- cell_freq$clustering_prop
  df_expr$cell <- "Unassigned"
  
  return(df_expr)
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











