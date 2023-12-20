TreeAnnotation_valuebox_UI <- function(id) {
  ns <- NS(id)
  column(width = 12,
  valueBox(width = 12,uiOutput(ns("progressBox")),"Percentage", color = "purple"),
  valueBox(width = 12,uiOutput(ns("progressBox2")),"Clusters", color = "purple"))
  
  }

TreeAnnotation_valuebox_Server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      # update valueboxes ----
      # browser()
      tmp=filterHM(session$userData$vars$hm,
                   session$userData$vars$myPickerPos,
                   session$userData$vars$myPickerNeg,
                   session$userData$vars$th)
      
      mydf=annotaionDF[rownames(tmp),]
      rownames(mydf)=rownames(tmp)
      mysum <- sum(mydf$clusterSize)
      
      output$progressBox <- renderText({
        format(mysum, digits = 3)
      })
      
      output$progressBox2 <- renderText({
        prettyNum(dim(mydf)[1], big.mark=",")
      })
      })}




# 

