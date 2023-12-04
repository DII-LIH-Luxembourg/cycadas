
TreeAnnotation_UI <- function(id) {
  ns <- NS(id)

  fluidRow(
    column(width = 4,
    box(width = NULL,title = "Create Node",
      textOutput(ns("progressBox")),
      textOutput(ns("progressBox2")),
      textInput(ns("newNode"), "Set Name..."),
      pickerInput(inputId = ns("parentPicker"),label = "Select Parent:",choices = NULL,
                  options = list(`actions-box` = TRUE,size = 10,`selected-text-format` = "count > 3"),
                  multiple = F)
      ),
    box(width = NULL,
        column(width = 4,
               checkboxGroupButtons(inputId = ns("treePickerPos"),label = "Positive:",choices = c("A"),direction = "vertical")
               ),
        column(width = 4,
               checkboxGroupButtons(inputId = ns("treePickerNeg"),label = "Negative:",choices = c("A"),direction = "vertical")
               )
        ),
    box(width = NULL,title = "Create New Node",
      actionButton(ns("createNodeBtn"), "Create Node")
      ),
    box(width = NULL,title = "Delete Node",
      actionButton(ns("deleteNodeBtn"), "Delete Node")
      ),
    box(width = NULL,title = "Export Annotation",
      downloadButton(ns("exportAnnotationBtn"), "Export Annotation")
      ),
    box(width = NULL,title = "Export Tree as Image",
      actionButton(ns("exportTreeGraphics"), "Export Tree Image")
      )
    ),
    # end column
  column(width = 8,
         box(width = NULL,title = "Annotation Tree",
             visNetworkOutput(ns("mynetworkid"))),
         box(width = NULL,title = "Heatmap",
             plotOutput(ns("hm_tree"))
             ),
         box(width = NULL,
             plotOutput(ns("umap_tree")))
  ))}

TreeAnnotation_Server <- function(id,reactVals,df01,tab) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      req(exists("df01"))
 
      TEST5 <<- reactVals
      print("Run Tree Annotation")
      
# Function: Update cluster labels ----
      updateClusterLabels <- function(mydf) {
        
        mysum <- sum(cell_freq[rownames(mydf),]$clustering_prop)
        
        output$progressBox <- renderText({
          paste0(mysum, "%")
        })
        output$progressBox2 <- renderText({
          paste0(dim(mydf)[1], "Cluster")
        })
      }
      
updatePickerInput(session,inputId = "parentPicker",choices = annotationlist)
      
# Function: Plot the annotation tree ----

        
      # TEST7 <<- reactVals$graph$nodes
      # TEST8 <<- reactVals$graph$edges
      # nodes <- TEST7
      # edges <- TEST8
      # 
      nodes <- reactVals$graph$nodes
      edges <- reactVals$graph$edges
      
      
        output$mynetworkid <- renderVisNetwork({
          
          for(i in 1:nrow(nodes)) {
            
            l <- nodes$label[i]
            
            # get the number of clusters with that label
            nLabel <- sum(df01Tree$cell == l)
            # mysum <- sum(cell_freq[rownames(df01Tree),]$clustering_prop)
            
            if (nLabel == 0) {
              nodes[i,]$color <- "grey"
            } else {
              nodes[i,]$color <- "blue"
            }
            
          }
          
          visNetwork(nodes, edges, width = "100%") %>%
            visEdges(arrows = "from") %>%
            visHierarchicalLayout() %>%
            visExport(type = "png", name = "export-network",
                      float = "left", label = "Save network", background = "white", style= "")
        })
        
        
      


        
        # heatmap <- reactVals$hm
        # heatmap_selected <- input$parentPicker

         
# Annotation heatmap plot ----

        
               
        
observeEvent(input$parentPicker, {Heatmap_module(id = "TreeAnnotation",df = reactVals,filter = input$parentPicker)})
 
  #       observeEvent(input$parentPicker, {
  # heatmap2 <<- reactVals$hm
  # heatmap2_selected <<- input$parentPicker
  # # heatmap2_selected <<- "CD4+ TN"
  # 
  #       output$hm_tree <- renderPlot({
  #         
  #         heatmap_matrix <- heatmap2 %>% 
  #           filter(cell == heatmap2_selected) %>% 
  #           select(-c("cell"))
  #           
  #           if(nrow(heatmap_matrix) > 0 & ncol(heatmap_matrix) > 0 ){
  #             pheatmap(heatmap_matrix,cluster_cols = F,cluster_rows = T)
  #             }
  # 
  #         
  #         # # req(!is.null(reactVals$hm))
  #         # # req(ncol(heatmap2)>0)
  #         # if (dim(heatmap2)[1] < 2) {
  #         #   # browser()
  #         #   # TEST1<<-reactVals$hm
  #         #   pheatmap(heatmap2 %>% select(-c("cell")), cluster_cols = F, cluster_rows = F)
  #         # } else {
  #         #   # message(dim(heatmap2)[1])
  #         #   # browser()
  #         #   # TEST1<<-reactVals$hm
  #           
  #       
  #             # }
  #       })
  #       })
        
        
        
        # unique(hm1$cell)
              
      
# Upload the marker expression file ----
      # observeEvent(c(input$fMarkerExpr, input$cluster_freq), {
        
        reactVals$graph <- initTree()
        # browser()
        
        # req(input$fMarkerExpr)
        # req(input$cluster_freq)
        
        # df <- read.csv(input$fMarkerExpr$datapath)
        # df_global <<- df
        # 
        # annotationlist <<- list("Unassigned")
        # 
        # cell_freq <<- read.csv(input$cluster_freq$datapath)
        # 
        # labels_row <-
        #   paste0(rownames(df), " (", cell_freq$clustering_prop , "%)")
        # 
        # marker_names <- rownames(df)
        # 
        # set.seed(1234)
        # 
        # my_umap <- umap(df)
        # dr_umap <<- data.frame(
        #   u1 = my_umap$layout[, 1],
        #   u2 = my_umap$layout[, 2],
        #   my_umap$data,
        #   cluster_number = 1:length(my_umap$layout[, 1]),
        #   check.names = FALSE
        # )
        # 
        # df01 <<- df %>% normalize01()
        # myDF <<- df %>% normalize01()
        # 
        # df01Tree <<- df %>% normalize01()
        # allMarkers <<- colnames(df)
        # df01Tree$cell <<- "Unassigned"
        # 
        # selectedMarkers <<- colnames(df)
        # posPickerList <<- colnames(df)

        
        TEST2 <<- colnames(df01)
        
        ## load only at start to fill the picker list
        updatePickerInput(
          session,
          inputId = "myPickerPos",
          label = "Select Positive Markers",
          choices = colnames(df01),
          # choices = NULL,
          options = list(
            `actions-box` = TRUE,
            size = 10,
            `selected-text-format` = "count > 3"
          )
        )
        updatePickerInput(
          session,
          inputId = "myPickerNeg",
          label = "Select Negative Markers",
          choices = colnames(df01),
          options = list(
            `actions-box` = TRUE,
            size = 10,
            `selected-text-format` = "count > 3"
          )
        )
        updateSelectInput(session, "markerSelect", "Select:", colnames(df01))
        
        updateCheckboxGroupButtons(
          session,
          inputId = "treePickerPos",
          choices = colnames(df01),
          selected = NULL
        )
        updateCheckboxGroupButtons(
          session,
          inputId = "treePickerNeg",
          choices = colnames(df01),
          selected = NULL
        )
        
        reactVals$th <- kmeansTH(df01)
        
        annotaionDF <<- data.frame("cell" = "unassigned",
                                   clusterSize = cell_freq$clustering_prop)
        at <<- reactiveValues(data = annotaionDF, dr_umap = dr_umap)
        
      # })
      
      
      
      
      
      # Observe MenutItems ----
      # eventReactive(tab, {
      #   
        # if (exists("df01")) {
        #   if(tab=="thresholds") { ## Thresholds ----
        #     
        #     if (is.null(reactVals$th) & !is.null(df01)) {
        #       reactVals$th <- kmeansTH(df01)
        #     }
        #   }
          # else if(tab=="treeannotation") { # Tree annotation ----
            

            # updatePickerInput(
            #   session,
            #   inputId = "updateNodePicker",
            #   choices = annotationlist,
            #   selected = ""
            # )
            # plotTree()
          # }
        # }
      # })
      
      observeEvent(input$treePickerPos, {
        
        updateCheckboxGroupButtons(
          session,
          inputId = "treePickerNeg",
          choices = colnames(df01),
          selected = input$treePickerNeg,
          disabledChoices = input$treePickerPos
        )
      })
      observeEvent(input$treePickerNeg, {
        
        updateCheckboxGroupButtons(
          session,
          inputId = "treePickerPos",
          choices = colnames(df01),
          selected = input$treePickerPos,
          disabledChoices = input$treePickerNeg
        )
      })
      
      # # Create new node ----
      # observeEvent(input$createNodeBtn, {
      #   
      #   if (is.null(input$treePickerPos) & is.null(input$treePickerNeg)) {
      #     print("no selection")
      #     showModal(modalDialog(
      #       title = "No Marker Selection",
      #       "Select positive and / or negative markers!",
      #       easyClose = TRUE,
      #       footer = NULL,
      #       return()
      #     ))
      #   } else if(input$newNode == "") {
      #     showModal(modalDialog(
      #       title = "Phenotype Name",
      #       "Set a Name for this Phenotype!",
      #       easyClose = TRUE,
      #       footer = NULL,
      #       return()
      #     ))
      #   }
      #   else {
      #     
      #     name <- input$newNode
      #     
      #     # make sure name is not yet taken
      #     if(name %in% reactVals$graph$nodes$label) {
      #       showModal(modalDialog(
      #         title = "Naming Error!",
      #         "The Name for this Phenotype is already taken!",
      #         easyClose = TRUE,
      #         footer = NULL
      #       ))
      #       
      #       return()
      #     }
      #     
      #     # browser()
      #     parent <- input$parentPicker
      #     
      #     # receive the parent settings, resp. parent hm
      #     # filter hm by parent cell name
      #     tmp <- df01Tree[df01Tree$cell == parent, ]
      #     tmp <- filterHM(tmp,input$treePickerPos, input$treePickerNeg, reactVals$th)
      #     
      #     # Make sure the selection is not empty!
      #     if(dim(tmp)[1] == 0) {
      #       
      #       showModal(modalDialog(
      #         title = "No Result",
      #         "The Selection for this Phenotype is empty!",
      #         easyClose = TRUE,
      #         footer = NULL
      #       ))
      #       return()
      #       
      #     } else {
      #       
      #       annotationlist <<- append(annotationlist, name)
      #       
      #       updatePickerInput(
      #         session,
      #         inputId = "parentPicker",
      #         choices = annotationlist
      #       )
      #       updatePickerInput(
      #         session,
      #         inputId = "updateNodePicker",
      #         choices = annotationlist,
      #         selected = ""
      #       )
      #       
      #       posmarker <- list(input$treePickerPos)
      #       negmarker <- list(input$treePickerNeg)
      #       
      #       reactVals$graph <- add_node(reactVals$graph,parent,name,posmarker,negmarker,"blue")
      #       
      #       df01Tree[rownames(tmp),]$cell <<- name
      #       
      #       updateCheckboxGroupButtons(
      #         session,
      #         inputId = "treePickerPos",
      #         choices = colnames(df01),
      #         selected = NULL,
      #         disabledChoices = input$treePickerPos
      #         
      #       )
      #       updateCheckboxGroupButtons(
      #         session,
      #         inputId = "treePickerNeg",
      #         choices = colnames(df01),
      #         selected = NULL,
      #         disabledChoices = input$treePickerNeg
      #         
      #       )
      #       updatePickerInput(
      #         session,
      #         inputId = "parentPicker",
      #         selected = name
      #       )
      #       
      #       updateClusterLabels(tmp)
      #       
      #       plotTree()
      #       
      #     }
      #     
      #   }
      # })
      
      
      

      
      
      
      
      
      
      # # Observe Parent Node selection ----
      # observeEvent(input$parentPicker, {
      #   
      #   # browser()
      #   
      #   node <- reactVals$graph$nodes %>% filter(label == input$parentPicker)
      #   myid <- node$id
      #   
      #   filterPosMarkers <- unlist(node$pm)
      #   filterNegMarkers <- unlist(node$nm)
      #   
      #   # collect all positive and negative markers
      #   # upwards from parent
      #   # go through edges until we reach one level below master node
      #   while (myid > 1) {
      #     print(myid)
      #     
      #     edge <- reactVals$graph$edges %>% filter(from == myid)
      #     next_id <- edge$to
      #     
      #     myid <- next_id
      #     next_node <- reactVals$graph$nodes %>% filter(id == next_id)
      #     filterPosMarkers <- append(filterPosMarkers, unlist(next_node$pm))
      #     filterNegMarkers <- c(filterNegMarkers, unlist(next_node$nm))
      #     
      #   }
      #   # remove the empty strings in the markers
      #   filterPosMarkers <- filterPosMarkers[nzchar(filterPosMarkers)]
      #   filterNegMarkers <- filterNegMarkers[nzchar(filterNegMarkers)]
      #   
      #   # tmp <- filterHM(df01Tree,filterPosMarkers, filterNegMarkers, reactVals$th)
      #   
      #   tmp <- filterHM(df01Tree,unique(unlist(filterPosMarkers)), unique(unlist(filterNegMarkers)), reactVals$th)
      #   
      #   tmp <- tmp[tmp$cell == node$label ,]
      #   
      #   updateClusterLabels(tmp)
      #   
      #   reactVals$hm <- tmp
      #   
      #   updateCheckboxGroupButtons(
      #     session,
      #     inputId = "treePickerPos",
      #     choices = colnames(df01),
      #     selected = NULL,
      #     disabledChoices = filterPosMarkers
      #   )
      #   updateCheckboxGroupButtons(
      #     session,
      #     inputId = "treePickerNeg",
      #     choices = colnames(df01),
      #     selected = NULL,
      #     disabledChoices = filterNegMarkers
      #   )
      #   
      #   # update umap plot for Tree ----
      #   ClusterSelection <- filterColor(df_global,tmp)
      #   
      #   output$umap_tree <-
      #     renderPlot(
      #       ggplot(dr_umap, aes(
      #         x = u1, y = u2, color = ClusterSelection
      #       )) +
      #         geom_point(size = 1.0) +
      #         theme_bw() +
      #         theme(legend.text = element_text(size = 8)) +
      #         guides(color = guide_legend(override.aes = list(size = 4)))
      #     )
      # })
      
      # # Observe node update picker ----
      # observeEvent(input$updateNodePicker, {
      #   
      #   browser()
      #   
      #   node <- reactVals$graph$nodes %>% filter(label == input$updateNodePicker)
      #   myid <- node$id
      #   
      #   filterPosMarkers <- unlist(node$pm)
      #   filterNegMarkers <- unlist(node$nm)
      #   
      #   updatePickerInput(
      #     session,
      #     inputId = "updatePickerPos",
      #     label = "Select Positive Markers",
      #     choices = colnames(df01),
      #     selected = filterPosMarkers,
      #     options = list(
      #       `actions-box` = TRUE,
      #       size = 10,
      #       `selected-text-format` = "count > 3"
      #     )
      #   )
      #   updatePickerInput(
      #     session,
      #     inputId = "updatePickerNeg",
      #     label = "Select Negative Markers",
      #     choices = colnames(df01),
      #     selected = filterNegMarkers,
      #     options = list(
      #       `actions-box` = TRUE,
      #       size = 10,
      #       `selected-text-format` = "count > 3"
      #     )
      #   )
      #   updateTextInput(
      #     session,
      #     inputId = "renameNode",
      #     label = NULL,
      #     value = node$label)
      # })
      
      # # Update Node ----
      # observeEvent(input$updateNodeBtn, {
      #   
      #   node <- reactVals$graph$nodes %>% filter(label == input$parentPicker)
      #   myid <- node$id
      #   
      #   reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker,]$pm <-
      #     list(c(reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker,]$pm[[1]], input$treePickerPos))
      #   
      #   reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker,]$nm <-
      #     list(c(reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker,]$nm[[1]], input$treePickerNeg))
      #   
      #   updateTreeAnnotation()
      #   
      #   reactVals$hm <- df01Tree[df01Tree$cell == input$parentPicker, ]
      #   
      #   # if "newNode" textfield is ot empty, update the node name
      #   if (input$newNode != "") {
      #     # browser()
      #     
      #     old_name <- input$parentPicker
      #     new_name <- input$newNode
      #     
      #     reactVals$graph$nodes[reactVals$graph$nodes$label == old_name,]$label <- new_name
      #     
      #     annotationlist[annotationlist == old_name] <<- new_name
      #     updatePickerInput(
      #       session,
      #       inputId = "parentPicker",
      #       selected = new_name,
      #       choices = annotationlist
      #     )
      #     
      #     # replace name in heatmap
      #     df01Tree$cell[df01Tree$cell == old_name] <<- new_name
      #     reactVals$hm <- df01Tree[df01Tree$cell == new_name, ]
      #     
      #     updateTextInput(
      #       session,
      #       inputId = "newNode",
      #       label = NULL,
      #       value = "",
      #       placeholder = NULL
      #     )
      #     
      #   }
      #   
      #   plotTree()
      # })
      
      # # Delete Node ----
      # observeEvent(input$deleteNodeBtn, {
      #   
      #   # browser()
      #   
      #   # node <- reactVals$graph$nodes %>% filter(label == input$updateNodePicker)
      #   node <- reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker, ]
      #   
      #   # check and make sure that this node is a leaf node
      #   if(TRUE %in% (reactVals$graph$edges$to == node$id)) {
      #     
      #     showModal(modalDialog(
      #       title = "Cannot delete Node",
      #       "The selected Node is not a Leaf Node!",
      #       easyClose = TRUE,
      #       footer = NULL
      #     ))
      #     
      #     return()
      #     
      #   }
      #   else {
      #     
      #     # before deleting the node, we must assign the cluster annotation
      #     # name for that node back to its parent
      #     parent_id <- reactVals$graph$edges$to[reactVals$graph$edges$from == node$id]
      #     parent_label <- reactVals$graph$nodes$label[reactVals$graph$nodes$id == parent_id]
      #     
      #     # In case there is a empty node with no the filtered clusters
      #     # we do not assign any labeling
      #     if( dim(df01Tree[df01Tree$cell == node$label,])[1] >0 ) {
      #       
      #       df01Tree[df01Tree$cell == node$label,]$cell <<- parent_label
      #       
      #     }
      #     
      #     annotationlist[annotationlist == node$label] <<- NULL
      #     
      #     updatePickerInput(
      #       session,
      #       inputId = "parentPicker",
      #       selected = NULL,
      #       choices = annotationlist
      #     )
      #     updatePickerInput(
      #       session,
      #       inputId = "updateNodePicker",
      #       selected = NULL,
      #       choices = annotationlist
      #     )
      #     updateTextInput(
      #       session,
      #       inputId = "renameNode",
      #       label = NULL,
      #       value = ""
      #     )
      #     updateTextInput(
      #       session,
      #       inputId = "newNode",
      #       label = NULL,
      #       value = ""
      #     )
      #     
      #     reactVals$graph <- delete_leaf_node(reactVals$graph, node$id)
      #     plotTree()
      #   }
      #   
      # })
      
      
      # importNodes <- function(id, label, pm, nm, parent_id, color) {
      #   
      #   if(id == parent_id) {return()}
      #   
      #   parent_row <- reactVals$graph$nodes %>% filter(id == parent_id)
      #   parent_name <- parent_row$label
      #   
      #   add_node(reactVals$graph, parent_name, label, color)
      # }
      
      
      # observeEvent(input$exportTreeGraphics, {
      #   
      #   # browser()
      #   
      #   # visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
      #   #   visEdges(arrows = "from") %>%
      #   #   visHierarchicalLayout() %>%
      #   #   visExport(type = "png", name = "export-network",
      #   #             float = "left", label = "Save network", background = "purple", style= "")
      #   
      #   
      # })
      
      # # Export Annotation Tree ----
      # output$exportAnnotationBtn <- downloadHandler(
      #   
      #   filename = function(){
      #     paste("cycadas_annotation_data_", Sys.Date(), ".zip", sep = "")
      #   },
      #   content = function(file) {
      #     temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
      #     dir.create(temp_directory)
      #     
      #     export_df_nodes <- reactVals$graph$nodes
      #     # Use apply() to concatenate the pm column for each row of the nodes dataframe
      #     pm_concatenated <- apply(export_df_nodes, 1, function(row) {
      #       
      #       pm_vector <- unlist(row["pm"])
      #       pm_vector <- paste(pm_vector, collapse = "|")
      #     })
      #     
      #     nm_concatenated <- apply(export_df_nodes, 1, function(row) {
      #       nm_vector <- unlist(row["nm"])
      #       nm_vector <- paste(nm_vector, collapse = "|")
      #     })
      #     
      #     # Add the concatenated pm column as a new column to the nodes dataframe
      #     export_df_nodes$pm <- pm_concatenated
      #     export_df_nodes$nm <- nm_concatenated
      #     
      #     download_list <- list(annTable = df01Tree,
      #                           nodesTable = export_df_nodes,
      #                           edgesTable = reactVals$graph$edges)
      #     # download_list <- list(annTable = df01Tree["cell"],
      #     #                       nodesTable = export_df_nodes,
      #     #                       edgesTable = reactVals$graph$edges)
      #     
      #     download_list %>%
      #       imap(function(x,y){
      #         if(!is.null(x)){
      #           file_name <- glue("{y}_data.csv")
      #           readr::write_csv(x, file.path(temp_directory, file_name))
      #         }
      #       })
      #     
      #     zip::zip(
      #       zipfile = file,
      #       files = dir(temp_directory),
      #       root = temp_directory
      #     )
      #   },
      #   contentType = "application/zip"
      # )
      
      
      # update Tree Annotation ----
      ## update the tree after a change in the thresholds
      ## walk through the list of nodes and adjust the annotation if
      ## needed based on the filtered DF
      
      
      ## Event on click scatter plot for setting the vertical line
      ##
      # observe click scatterplot -----
      ##
      
      ##
      # Thresholds plotting function ----
      ##
      
      ##
      # Server - Annotations Tab ------------------------------------------------
      ##
      # Load the annotation file ---------------------------------------------
      # observeEvent(input$annTable,{
      #   annData <<- read.csv(input$annTable$datapath)
      #   annData$X <- NULL
      #   at$data <- annData
      # })
      
      
    })}

Heatmap_module <- function(id,df,filter) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
    print("run heatmaps plot") 
      print(filter)
      print(df$hm)
      
      
    output$hm_tree <- renderPlot({
      heatmap_matrix <<- df$hm %>% 
          filter(cell == filter) %>% 
          select(-c("cell"))
        
        if(nrow(heatmap_matrix) > 0 & ncol(heatmap_matrix) > 0 ){
          pheatmap(heatmap_matrix,cluster_cols = F,cluster_rows = T)
        } else {plot(1,1)}
      print(paste0(nrow(heatmap_matrix),ncol(heatmap_matrix)))  
      })
      })}
