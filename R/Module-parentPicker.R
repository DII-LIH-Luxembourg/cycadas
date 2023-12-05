
parentPicker_UI <- function(id) {
  ns <- NS(id)
  pickerInput(inputId = ns("parentPicker"),label = "Select Parent:",choices = NULL,
              options = list(`actions-box` = TRUE,size = 10,`selected-text-format` = "count > 3"),
              multiple = F)
}

parentPicker_Server <- function(id,annotationlist) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
 
      print("Run Parent Picker")
 
      
      
      
      return(
        list(
          parentPicker = reactive({ input$parentPicker }))
      )
      
      
      
    })}