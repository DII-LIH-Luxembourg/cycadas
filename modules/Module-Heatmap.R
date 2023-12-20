
Heatmap_UI <- function(id) {
  ns <- NS(id)
    plotOutput(ns("hm_tree"))
  }

Heatmap_Server <- function(id,filter) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("run heatmaps plot")
      output$hm_tree <- renderPlot({
        
#TODO:  # expression
        # scaling
        # filter 
        
        
        
        heatmap_matrix <- session$userData$vars$hm %>%
          # filter(cell == filter) %>%
          select(-c("cell"))
        if(nrow(heatmap_matrix) > 0 & ncol(heatmap_matrix) > 0 ){
          pheatmap(heatmap_matrix,cluster_cols = F,cluster_rows = T)
        } else {plot(2,2)}
        })
    })}
