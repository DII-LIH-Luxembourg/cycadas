
Umap_UI <- function(id) {
  ns <- NS(id)
  plotOutput(ns("umap_tree"))
  }

Umap_Server <- function(id,reactVals,filter) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("run umap Plot")
      
      node <- reactVals$graph$nodes %>% 
        filter(label == filter)
      
      myid <- node$id
       
      print("node")
      print(head(node))
      print(myid)

      if(length(myid)==0){
        output$umap_tree <-renderPlot(plot(1,1))
        }else{
      
      filterPosMarkers <- unlist(node$pm)
      filterNegMarkers <- unlist(node$nm)
      
      # collect all positive and negative markers
      # upwards from parent
      # go through edges until we reach one level below master node
      while (myid > 1) {
        print(myid)
        
        edge <- reactVals$graph$edges %>% filter(from == myid)
        next_id <- edge$to
        
        myid <- next_id
        next_node <- reactVals$graph$nodes %>% filter(id == next_id)
        filterPosMarkers <- append(filterPosMarkers, unlist(next_node$pm))
        filterNegMarkers <- c(filterNegMarkers, unlist(next_node$nm))
        
      }
      # remove the empty strings in the markers
      filterPosMarkers <- filterPosMarkers[nzchar(filterPosMarkers)]
      filterNegMarkers <- filterNegMarkers[nzchar(filterNegMarkers)]

      # tmp <- filterHM(df01Tree,filterPosMarkers, filterNegMarkers, reactVals$th)

      tmp <- filterHM(df01Tree,unique(unlist(filterPosMarkers)), unique(unlist(filterNegMarkers)), reactVals$th)
      
      tmp <- tmp[tmp$cell == node$label ,]
      
      # Function: Update cluster labels ----
      updateClusterLabels <- function(mydf) {
        
        mysum <- sum(cell_freq[rownames(mydf),]$clustering_prop)
        
        output$progressBox <- renderText({
          paste0(mysum, "%")
        })
        output$progressBox2 <- renderText({
          paste0(dim(mydf)[1], "Cluster")
        })
      }
      
      
      updateClusterLabels(tmp)
      #
      reactVals$hm <- tmp

      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = colnames(df01),
        selected = NULL,
        disabledChoices = filterPosMarkers
      )
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = colnames(df01),
        selected = NULL,
        disabledChoices = filterNegMarkers
      )
      
      # update umap plot for Tree ----
      
      filterColor <- function(DF,hm) {
        # browser()
        tmp <- replicate(nrow(DF), "other clusters")
        tmp[as.numeric(rownames(hm))] <- "selected phenotype"
        
        return(unname(tmp))
        
      }
      
      ClusterSelection <<- filterColor(df_global,tmp)

            output$umap_tree <-
        renderPlot(
          
          ggplot(dr_umap, aes(
            x = u1, y = u2, color = ClusterSelection)) +
            geom_point(size = 1.0) +
            theme_bw() +
            theme(legend.text = element_text(size = 8)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
        
          
          
          )
      
      }
      
       
    })}
