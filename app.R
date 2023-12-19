
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
# source("modules/Module-TreeAnnotation.R")
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

  # browser()



  # constructs a string of positive or negative markers
  ph_name <<- ""

  # Server function ----
  server = function(input, output, session) {
    
    session$userData$vars <- reactiveValues(th = NULL)
    session$userData$vars <- reactiveValues(myTH = NULL)
    session$userData$vars <- reactiveValues(md = NULL)
    session$userData$vars <- reactiveValues(counts_table = NULL)
    session$userData$vars <- reactiveValues(DA_result_table = NULL)
    session$userData$vars <- reactiveValues(DA_interactive_table = NULL)
    session$userData$vars <- reactiveValues(graph = NULL)
    session$userData$vars <- reactiveValues(hm = NULL)
    session$userData$vars <- reactiveValues(merged_prop_table = NULL)
    session$userData$vars <- reactiveValues(treePickerPos = NULL)
    session$userData$vars <- reactiveValues(treePickerNeg = NULL)
    session$userData$vars <- reactiveValues(annotationlist = NULL)

    
    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    
        
    # Load the marker threshold file ---------------------------------------
    observeEvent(input$fTH,{

      # browser()

      th <<- read.csv(input$fTH$datapath)
      th$X <- NULL
      th$color <- "blue"

      th[th$bi_mod < 0.555, "color"] <- "red"

      rownames(th) <- th$cell

      session$userData$vars$th<- th

    })

   
    ## create Phenotype name from picker e.g. CD4+CD8+CD19 --------------------
    observeEvent(input$myPickerPos, {
      output$r1 <- renderText(setPhenotypeName(input$myPickerPos, "pos", ph_name))
    })
    observeEvent(input$myPickerNeg, {
      output$r2 <- renderText(setPhenotypeName(input$myPickerNeg, "neg", ph_name))
    })

    ## Annotation Table -------------------------------------------------------
    observeEvent(c(input$myPickerPos, input$myPickerNeg, session$userData$vars$th),{
      req(input$fMarkerExpr)
      req(input$cluster_freq)

      tmp=filterHM(df01,input$myPickerPos, input$myPickerNeg, session$userData$vars$th)

      # update myDF to the filtered HM
      myDF <<- tmp

      mydf=annotaionDF[rownames(tmp),]
      rownames(mydf)=rownames(tmp)

      # update heatmap plot ----
      output$hm <- renderPlot({
        if (dim(tmp)[1] < 2) {
          pheatmap(tmp, cluster_cols = F, cluster_rows = F)
        } else {
          message(dim(tmp)[1])
          pheatmap(tmp, cluster_cols = F)
        }
      })

      # update umap plot ----
      myColor <- filterColor(df_global,tmp)
      Temp1<<-myColor
      Temp2<<-dr_umap

      output$umap <-
        renderPlot(
          ggplot(dr_umap, aes(
            x = u1, y = u2, color = myColor
          )) +
            geom_point(size = 1.0) +
            theme_bw() +
            theme(legend.text = element_text(size =
                                               8)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
        )
      # update table ----
      output$tableAnnotation <- DT::renderDT(
        at$data,
        extensions = c('Buttons', 'Scroller'),
        options = list(
          dom = 'Bfrtip',
          buttons = c('csv'),
          paging = FALSE
        )
      )

      # update valueboxes ----
      mysum <- sum(mydf$clusterSize)
      output$progressBox <- renderValueBox({
        valueBox(paste0(mysum, "%"), "Selected", icon = icon("list"),color = "purple")
      })
      output$progressBox2 <- renderValueBox({
        valueBox(dim(mydf)[1], "Cluster", icon = icon("list"),color = "purple")
      })
    })

    # Set Phenotype names -----------------------------------------------------
    observeEvent(input$btnSetType, {

      # if cluster has already been assigned, add the new name after
      for (n in rownames(myDF)) {
        if (at$data[n, "cell"] == "unassigned") {
          at$data[n, "cell"] <- input$phenotype
        } else(
          at$data[n, "cell"] <- paste0(at$data[n, "cell"], "_", input$phenotype)
        )

      }
    })

    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################

    
    
# Load data based on pressing the unannotated data button
    observe({Settings_Server1(id="Settings",reactVals=reactVals)}) %>% 
      bindEvent(input$btnLoadDemoData)

# Load data based on pressing the annotated data button ----
    
    observe({Settings_Server2(id="Settings",reactVals=reactVals)}) %>% 
      bindEvent(input$btnLoadAnnoData)
    
# Load data based on pressing the Import Tree button  ----
    
    observeEvent(input$btnImportTree, {Settings_Server3(id="Settings")})
    
# Lead threshold Tab    ----
    
    observe({
      req(!is.null(session$userData$vars$th))
      threshold_Server(id="threshold")})
    
# # Run the Tree Annotation Tab based on pressing the Annotated demo button  ----
#     observe({
#       # req(exists("cellfreq"))
#       # req(exists("df01"))
#       TreeAnnotation_Server(id="TreeAnnotation",df01 = df01)}) %>%
#       bindEvent(c(input$btnLoadAnnoData))

    
# Create nodes ----
    
    observe({
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
      bindEvent(c(input$btnLoadAnnoData,input$createNodeBtn))   
    
# Update Parent picker based on annotation list once Annotated demo is pressed ----
    
    observe({
      req(!is.null("annotationlist"))
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
      Heatmap_Server(id="Heatmap",filter=input$parentPicker)
      visNetwork_Server(id="visNetwork",filter=input$parentPicker)
      Umap_Server(id="Umap",filter=input$parentPicker)}) %>% 
      bindEvent(input$parentPicker)
    
# Run Marker Expression Tab ----
    observe({
      UMAP_ME_Server(id="UMAP_ME")}) %>% 
      bindEvent(input$btnLoadAnnoData)
    
# Run Differential Abundance Tab ----
    observe({
      Differential_Abundance_Server(id="Differential_Abundance")}) %>% 
      bindEvent(input$btnLoadAnnoData)
  
# Run DA_interactive Tab ----
    observe({
      output$interactiveTree <- renderVisNetwork({
        visNetwork(session$userData$vars$graph$nodes, session$userData$vars$graph$edges, width = "100%") %>%
          visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
      })}) %>% 
      bindEvent(input$btnLoadAnnoData)
    
    observe({
      print(input$current_node_id)
      DA_interactive_Server(id="DA_interactive",current_node_id=input$current_node_id)}) %>% 
      bindEvent(input$current_node_id)
    
# Run umap interactive Tab ----
  observe({
    umap_interactive_Server(id="umap_interactive")}) %>% 
    bindEvent(input$btnLoadAnnoData)
  }
  
  
  shinyApp(
    ui = ui,
    server = server
  )
# }

