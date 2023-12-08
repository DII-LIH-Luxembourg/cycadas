
Threshold_UI <- function(id) {
  ns <- NS(id)
  fluidRow(column(
    width = 6,
    box(
      width = NULL,
      title = "Marker Expression:",
      plotOutput(ns("plot"), click = ns("plot_click"))
    ),
    box(width = NULL,
        title = "Histogram:",
        plotOutput(ns("plot2")))
  ),
  column(
    width = 6,
    # box(
    #   width = NULL,
    #   id = "transformationBox",
    #   title = "Data transformation",
    #   radioButtons(
    #     "radio",
    #     label = NULL,
    #     choices = list("0 to 1" = "1"),
    #     selected = "1"
    #   )
    # ),
    box(width = NULL,
        title = "Markers",
        DTOutput(ns('table')))
  ))
}


# ThresholdPlot_Server <- function(id,reactVals,selectedid) {
#   moduleServer(
#     id = id,
#     module = function(input, output, session) {
# 
# 
#       })}




threshold_Server <- function(id,reactVals) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      req(reactVals$th)
      print ("threshold_Server run")
      output$table = renderDataTable(
        reactVals$th[, c("cell", "threshold", "bi_mod")],
        editable = F,
        extensions = c('Buttons', 'Scroller'),
        selection = 'single',
        # selection = list(selection='single',selected = 5),
        options = list(
          dom = 'Bfrtip',
          server = FALSE,
          deferRender = TRUE,
          scrollY = 500,
          scroller = TRUE,
          buttons = list(
            list(extend = 'csv', filename = "MarkerThresholds")
          )
        )
      )
      
 
      
      observeEvent(c(input$plot_click,input$table_rows_selected), {
        
        req(input$table_rows_selected)
        req(reactVals$th)
        
        TEST <<- input$plot_click
        
        if(!is.null(input$plot_click)){
          print("Threshold click")
          reactVals$th[input$table_rows_selected, "threshold"] <- round(input$plot_click$x, 3)}
        
        selectedid <- input$table_rows_selected
        selRow <- reactVals$th[selectedid,]
        marker <- reactVals$th[selectedid, "cell"]
        # threshold value for vertical line
        myTH <- reactVals$th[selectedid, "threshold"]
        myCol <- reactVals$th[selectedid, "color"]
        marker_expr <- getMarkerDistDF(marker, "1")
        me<-marker_expr
        
        
        output$plot <- renderPlot({
          set.seed(1)
          
          ggplot(me, aes_string(x=me[,1], y=me[,2])) +
            geom_point(size=1) +
            theme(axis.title.y = element_blank(),
                  axis.ticks.y  = element_blank(),
                  axis.text.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank()) +
            labs(x = "Scale 0 to 1",title = marker ) +
            geom_vline(xintercept = myTH, linetype="dotted",
                       color = myCol, size=1.5)
        })
        
        output$plot2 <- renderPlot({
          
          ggplot(me, aes_string(x = me[, 1])) +
            geom_histogram(bins = 80) +
            labs(x = "Scale 0 to 1") +
            geom_vline(
              xintercept = myTH,
              linetype = "dotted",
              color = myCol,
              size = 1.5
            )
        })

        updateTreeAnnotation_server(id="updateTreeAnnotation")
        # updateTreeAnnotation(reactVals$th[input$table_rows_selected])
        
        
        
      })
    })}



updateTreeAnnotation_server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      df01Tree$cell <- "Unassigned"
      
      # for all rows in annotation table
      # get all parents for a row:
      if (nrow(reactVals$graph$nodes) > 1) {
        reactVals$graph$nodes$to <- reactVals$graph$edges$to
        
        # Iterate through the dataframe row by row
        for (i in 1:nrow(reactVals$graph$nodes)) {
          
          posMarker <- list()
          negMarker <- list()
          
          nodeID <- reactVals$graph$nodes$id[i]
          
          # now for that nodeID, get all parents and collect
          # the pm and nm markers
          while (nodeID > 1) {
            
            posMarker <- c(posMarker, unlist(reactVals$graph$nodes$pm[reactVals$graph$nodes$id == nodeID]))
            negMarker <- c(negMarker, unlist(reactVals$graph$nodes$nm[reactVals$graph$nodes$id == nodeID]))
            
            nodeID <- reactVals$graph$edges$to[reactVals$graph$edges$from == nodeID]
            
          }
          # now we have all positive and negative marker for that type collected
          # and we can start filtering the df
          # remove the empty strings in the markers
          posMarker <- posMarker[nzchar(posMarker)]
          negMarker <- negMarker[nzchar(negMarker)]
          
          tmp <- filterHM(df01Tree,unique(unlist(posMarker)), unique(unlist(negMarker)), reactVals$th)
          # tmp <- filterHM(df01Tree,unique(unlist(posMarker)), unique(unlist(negMarker)), my_th)
          
          df01Tree[rownames(tmp), 'cell'] <<- reactVals$graph$node$label[i]
          
          # updateClusterLabels(tmp)
        }
      }
      
      
      
    })}



