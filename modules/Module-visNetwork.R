
visNetwork_UI <- function(id) {
  ns <- NS(id)
  visNetworkOutput(ns("mynetworkid"))
  }

visNetwork_Server <- function(id,reactVals,filter) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("run visual Network plot")
      nodes <- reactVals$graph$nodes
      edges <- reactVals$graph$edges
      
      print (head(nodes))
      
      output$mynetworkid <- renderVisNetwork({
        
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
