umap_interactive_UI <- function(id) {
  ns <- NS(id)
  c(
  fluidRow(column(width = 6,
                  box(
                    width = NULL,
                    plotOutput(ns("umap2"), brush = ns("umap2_brush"))
                  )),
           column(width = 6,
                  box(
                    width = NULL,
                    plotOutput(ns("hm2"))
                  ))),
  fluidRow(column(width = 12,
                  box(
                    width = NULL,
                    DTOutput(ns("umap_data"))
                  ))))
  }

umap_interactive_Server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("Run umap_interactive Tab")

      # Server - UMAP interactive Tab -------------------------------------------
      output$umap2 <- renderPlot(
        ggplot(dr_umap, aes(x = u1, y = u2)) +
          geom_point(size = 1.0) +
          theme_bw() +
          theme(legend.text=element_text(size=8)) +
          guides(color = guide_legend(override.aes = list(size = 4)))
      )
      
      output$hm2 <- renderPlot({
        my_new_hm <- round(brushedPoints(dr_umap, input$umap2_brush), 2)
        
        if (dim(my_new_hm)[1] > 0) {
          
          pheatmap(my_new_hm %>% select(-c("cluster_number", "u1", "u2")), cluster_cols = F)
        } else {
          ggplot() + theme_void() + ggtitle("Select area on Umap to plot Heatmap")
        }
      })
      
      output$umap_data <-
        DT::renderDT(server = FALSE, {
          DT::datatable(
            round(brushedPoints(dr_umap, input$umap2_brush), 2),
            filter = 'top',
            extensions = 'Buttons',
            options = list(
              scrollY = 600,
              scrollX = TRUE,
              dom = '<"float-left"l><"float-right"f>rt<"row"<"col-sm-4"B><"col-sm-4"i><"col-sm-4"p>>',
              lengthMenu =  list(c(10, 25, 50,-1),
                                 c('10', '25', '50', 'All')),
              scrollCollapse = TRUE,
              lengthChange = TRUE,
              widthChange = TRUE,
              rownames = TRUE
            )
          )
        })
    })}