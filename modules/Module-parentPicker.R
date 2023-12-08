
parentPicker_UI <- function(id) {
  # ns <- NS(id)
  selectizeInput(inputId = ("parentPicker"),label = "Select Parent:",choices=NULL)
  
  # selectInput(inputId = ns("parentPicker"),label = "Select Parent:",choices = NULL)
  # ,
  #             options = list(`actions-box` = TRUE,size = 10,`selected-text-format` = "count > 3"),
  #             multiple = F)
}

parentPicker_Server <- function(id,options) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      # ns <- NS(id)
 
      print("Run Parent Picker")
      TEST<<-options
      # print(options)
      print ("update paraent picker")
      # updateSelectInput(session,inputId = ns("parentPicker"),choices = options) 
      

       
          # reactVals$ParentValue <- input$parentPicker
      # return(input$parentPicker)
      session$userData$vars$parentPicker <- input$parentPicker
      
    })}