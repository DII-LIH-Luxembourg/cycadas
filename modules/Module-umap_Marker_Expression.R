
UMAP_ME_UI <- function(id) {
  ns <- NS(id)
  
  fluidRow(column(width = 10,
                  box(width = NULL,plotOutput(ns("umap3"))
                  )),
           column(width = 2,
                  box(width = NULL,selectInput(ns("markerSelect"), "Select:", choices = NULL)
                  )))

  
  
}

UMAP_ME_Server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("Run Umap Marker expression Tab")

      observe({
        updateSelectInput(session, "markerSelect", "Select:", session$userData$vars$th$cell)
      })
      
      # Server - UMAP Marker Expression Tab ----
      output$umap3 <-
        renderPlot(
          ggplot(dr_umap, aes_string(
            x = "u1",
            y = "u2",
            color = input$markerSelect
          )) +
            geom_point(size = 1.0) +
            theme_bw() +
            scale_color_gradientn(input$markerSelect,
                                  colours = colorRampPalette(rev(brewer.pal(
                                    n = 11, name = "Spectral"
                                  )))(50))
        )
       
    })}