
TreeAnnotation_UI <- function(id) {
  ns <- NS(id)

  fluidRow(
    column(width = 4,
    box(width = NULL,title = "Create Node",
      textOutput(("progressBox")),
      textOutput(("progressBox2")),
      textInput("newNode", "Set Name..."),
      pickerInput(inputId = "parentPicker",label = "Select Parent:",choices = NULL,
                  options = list(`actions-box` = TRUE,size = 10,`selected-text-format` = "count > 3"),
                  multiple = F)
      ),
    box(width = NULL,
        column(width = 4,
               checkboxGroupButtons(inputId = "treePickerPos",label = "Positive:",choices = c("A"),direction = "vertical")
               ),
        column(width = 4,
               checkboxGroupButtons(inputId = "treePickerNeg",label = "Negative:",choices = c("A"),direction = "vertical")
               )
        ),
    box(width = NULL,title = "Create New Node",
      actionButton("createNodeBtn", "Create Node")
      ),
    box(width = NULL,title = "Delete Node",
      actionButton("deleteNodeBtn", "Delete Node")
      ),
    box(width = NULL,title = "Export Annotation",
      downloadButton("exportAnnotationBtn", "Export Annotation")
      ),
    box(width = NULL,title = "Export Tree as Image",
      actionButton("exportTreeGraphics", "Export Tree Image")
      )
    ),
    # end column
  column(width = 8,
         box(width = NULL,title = "Annotation Tree",
             visNetworkOutput("mynetworkid")
             ),
         box(width = NULL,title = "Heatmap",
             plotOutput("hm_tree")
             ),
         box(width = NULL,
             plotOutput("umap_tree"))
  ))}

TreeAnnotation_Server <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
 
      
      
      
      
    })}