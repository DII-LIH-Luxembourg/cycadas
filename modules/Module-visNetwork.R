
visNetwork_UI <- function(id) {
  ns <- NS(id)
  visNetworkOutput(ns("mynetworkid"))
  }

visNetwork_Server <- function(id,reactVals,filter) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("run visual Network plot")
      nodes <- session$userData$vars$graph$nodes
      edges <- session$userData$vars$graph$edges
      
      print (head(nodes))
      # req(nrow(session$userData$vars$graph$nodes)>0)
      
      output$mynetworkid <- renderVisNetwork({
        
        # browser()
        
        
        for(i in 1:nrow(nodes)) {
          
          l <- nodes$label[i]
          
          # get the number of clusters with that label
          nLabel <- sum(df01Tree$cell == l)
          # mysum <- sum(cell_freq[rownames(df01Tree),]$clustering_prop)
          
          if (nLabel == 0) {
            nodes[i,]$color <- "grey"
          } else {
            nodes[i,]$color <- "blue"
          }
          
        }
        
        visNetwork(nodes, edges, width = "100%") %>%
          visEdges(arrows = "from") %>%
          visHierarchicalLayout() %>%
          visExport(type = "png", name = "export-network",
                    float = "left", label = "Save network", background = "white", style= "")
      })
      
    })}
