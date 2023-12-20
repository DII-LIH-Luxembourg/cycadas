
# List of packages you want to check and install if needed
packages_to_install <- c("shiny", "DT", "ggplot2", "matrixStats", "tidyverse", "stats",
                         "pheatmap", "Ckmeans.1d.dp", "umap", "RColorBrewer", "shinydashboard",
                         "shinyWidgets", "visNetwork", "glue", "purrr", "reshape2", "mousetrap")

# Check if each package is already installed, and install if not
for (package in packages_to_install) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }

}


source("modules/utils.R")
source("modules/Module-Settings.R")
source("modules/Module-Threshold.R")
source("modules/Module-TreeAnnotation_valuebox.R")
source("modules/Module-TreeAnnotation_picker.R")
source("modules/Module-TreeAnnotation_createNode.R")
source("modules/Module-TreeAnnotation_annotationdownload.R")
source("modules/Module-visNetwork.R")
source("modules/Module-Heatmap.R")
source("modules/Module-umap.R")
source("modules/Module-DeleteNode.R")
source("modules/Module-umap_Marker_Expression.R")
source("modules/Module-Differential_Abundance.R")
source("modules/Module-DA_intreractive.R")
source("modules/Module-umap_interactive.R")
source("modules/ui.R")

# cycadas <- function() {

  # Server function ----
  server = function(input, output, session) {
    
    session$userData$vars <- reactiveValues(median_expr = NULL,
                                            th = NULL,
                                            myTH = NULL,
                                            md = NULL,
                                            counts_table = NULL,
                                            DA_result_table = NULL,
                                            DA_interactive_table = NULL,
                                            graph = NULL,
                                            hm = NULL,
                                            merged_prop_table = NULL,
                                            treePickerPos = NULL,
                                            treePickerNeg = NULL,
                                            annotationlist = NULL)

# Load data based on pressing the unannotated data button
    #TODO:remove the reactVals and check it runs still
    observe({Settings_Server1(id="Settings",reactVals=reactVals)}) %>% 
      bindEvent(input$btnLoadDemoData)

# Load annotated dataset ----
    
    observe({Settings_Server2(id="Settings",reactVals=reactVals)}) %>% 
      bindEvent(input$btnLoadAnnoData)
    
# Load data based on pressing the Import Tree button  ----
    
    observeEvent(input$btnImportTree, {Settings_Server3(id="Settings")})

# Load threshold Tab    ----
    
    observe({
      req(!is.null(session$userData$vars$th))
      threshold_Server(id="threshold")})
    

# Create nodes ----
    
    observe({
      req(!is.null(session$userData$vars$median_expr))
      print("Create Nodes")
      TreeAnnotation_createNode_Server(id="TreeAnnotation_createNode",parent = input$parentPicker,newNodeName=input$newNode)
      updateSelectizeInput(session, ("parentPicker"),selected = NULL, choices =      session$userData$vars$annotationlist, server = TRUE)}) %>% 
      bindEvent(c(input$createNodeBtn))  
    
    
# Delete nodes ----
    
    observe({
      req(exists("df01Tree"))
      DeleteNode_Server(id="DeleteNode",filter = input$parentPicker)
      updateSelectizeInput(session, ("parentPicker"),selected = NULL, choices = session$userData$vars$annotationlist, server = TRUE)
      updatePickerInput(session,inputId = "updateNodePicker",selected = NULL,choices = session$userData$vars$annotationlist)
      updateTextInput(session,inputId = "renameNode",label = NULL,value = "")
      updateTextInput(session,inputId = "newNode",label = NULL,value = "")})  %>%
      bindEvent(c(input$deleteNodeBtn))    
 
    
# Run the Tree Annotation Tab picker  ----
    observe({
      # browser()
      TreeAnnotation_picker_Server(id="TreeAnnotation_picker")}) %>%
      bindEvent(c(input$btnLoadDemoData,input$btnLoadAnnoData,input$createNodeBtn))   
    
# Update Parent picker based on annotation list once Annotated demo is pressed ----
    
    observe({
      print("update Parent Picker")
      updateSelectizeInput(session, ("parentPicker"), choices = session$userData$vars$annotationlist, server = TRUE)
      })  %>%
    bindEvent(c(input$btnLoadAnnoData,input$btnLoadDemoData,input$deleteNodeBtn,input$createNodeBtn))
    
# Update TreeAnnotation_annotationdownload  ----
    observe({
      TreeAnnotation_annotationdownload_Server(id="TreeAnnotation_annotationdownload")}) %>% 
      bindEvent(input$btnLoadAnnoData,input$btnLoadDemoData,input$deleteNodeBtn,input$createNodeBtn)   
    
    
    
# Update Tree annotation visuals based on parent picker ----
    observe({
      req(!is.null(session$userData$vars$hm))
      TreeAnnotation_valuebox_Server(id="TreeAnnotation_valuebox")
      Heatmap_Server(id="Heatmap",filter=input$parentPicker)
      visNetwork_Server(id="visNetwork",filter=input$parentPicker)
      Umap_Server(id="Umap",filter=input$parentPicker)}) %>% 
      bindEvent(input$parentPicker)
    
# Run Marker Expression Tab ----
    observe({
      UMAP_ME_Server(id="UMAP_ME")}) %>% 
      bindEvent(input$btnLoadDemoData,input$btnLoadAnnoData)
    
# Run Differential Abundance Tab ----
    observe({
      Differential_Abundance_Server(id="Differential_Abundance")}) %>% 
      bindEvent(input$btnLoadDemoData,input$btnLoadAnnoData)
  
# Run DA_interactive Tab ----
    observe({
      output$interactiveTree <- renderVisNetwork({
        visNetwork(session$userData$vars$graph$nodes, session$userData$vars$graph$edges, width = "100%") %>%
          visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
      })}) %>% 
      bindEvent(input$btnLoadDemoData,input$btnLoadAnnoData)
    
    observe({
      print(input$current_node_id)
      DA_interactive_Server(id="DA_interactive",current_node_id=input$current_node_id)}) %>% 
      bindEvent(input$current_node_id)
    
# Run umap interactive Tab ----
  observe({
    umap_interactive_Server(id="umap_interactive")}) %>% 
    bindEvent(input$btnLoadDemoData,input$btnLoadAnnoData)
  }
  
  
  shinyApp(
    ui = ui,
    server = server
  )
# }

