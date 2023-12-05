
Heatmap_UI <- function(id) {
  ns <- NS(id)
    plotOutput(ns("hm_tree"))
  }

Heatmap_Server <- function(id,df,filter) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
     
      print("run heatmaps plot")
      # req(df)
      # req(filter)
          print("run heatmaps plot")
            print(filter)
            # print(df)
            # print(Sys.getenv)


          output$hm_tree <- renderPlot({
            plot(1,1)
            
            
            # heatmap_matrix <- df$hm %>%
            #     filter(cell == filter) %>%
            #     select(-c("cell"))
            # 
            #   if(nrow(heatmap_matrix) > 0 & ncol(heatmap_matrix) > 0 ){
            #     pheatmap(heatmap_matrix,cluster_cols = F,cluster_rows = T)
            #   } else {plot(1,1)}
            # print(paste0(nrow(heatmap_matrix),ncol(heatmap_matrix)))
            })
    })}
