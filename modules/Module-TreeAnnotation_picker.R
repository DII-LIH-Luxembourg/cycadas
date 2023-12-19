
TreeAnnotation_picker_UI <- function(id) {
  ns <- NS(id)
  splitLayout(cellWidths = c("50%", "50%"),
              checkboxGroupButtons(inputId = ns("treePickerPos"),label = "Positive:",choices = c("A"),direction = "vertical"),
              checkboxGroupButtons(inputId = ns("treePickerNeg"),label = "Negative:",choices = c("A"),direction = "vertical"))
  }

TreeAnnotation_picker_Server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("run TreeAnnotation_picker")

      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = session$userData$vars$th$cell,
        selected = NULL
      )
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = session$userData$vars$th$cell,
        selected = NULL
      )
      
      
observeEvent(input$treePickerPos, {
  session$userData$vars$treePickerPos <- input$treePickerPos
  session$userData$vars$treePickerNeg <- input$treePickerNeg   
  updateCheckboxGroupButtons(
          session,
          inputId = "treePickerNeg",
          choices = colnames(df01),
          selected = input$treePickerNeg,
          disabledChoices = input$treePickerPos
        )
  })
observeEvent(input$treePickerNeg, {
  session$userData$vars$treePickerPos <- input$treePickerPos
  session$userData$vars$treePickerNeg <- input$treePickerNeg   
  updateCheckboxGroupButtons(
          session,
          inputId = "treePickerPos",
          choices = colnames(df01),
          selected = input$treePickerPos,
          disabledChoices = input$treePickerNeg
        )
      })
   

    })}
