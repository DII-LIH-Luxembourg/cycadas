TreeAnnotation_UI <- function(id) {
  ns <- NS(id)
    c(box(width = NULL,
        column(width = 4,
               checkboxGroupButtons(inputId = ns("treePickerPos"),label = "Positive:",choices = c("A"),direction = "vertical")
               ),
        column(width = 4,
               checkboxGroupButtons(inputId = ns("treePickerNeg"),label = "Negative:",choices = c("A"),direction = "vertical")
               )
        ),
      box(width = NULL,title = "Create Node",
          textOutput(ns("progressBox")),
          textOutput(ns("progressBox2")),
          textInput(ns("newNode"), "Set Name..."),
          actionButton(ns("createNodeBtn"), "Create Node")),
      box(width = NULL,title = "Export Annotation",
          downloadButton(ns("exportAnnotationBtn"), "Export Annotation")
      ))
  }

TreeAnnotation_Server <- function(id,reactVals,df01,tab) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      # req(exists("df01"))
 
      # TEST5 <<- reactVals
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
      
# Upload the marker expression file ----
      # observeEvent(c(input$fMarkerExpr, input$cluster_freq), {
        
        session$userData$vars$graph <- initTree()
        # browser()
        
        # req(input$fMarkerExpr)
        # req(input$cluster_freq)
        
        df <- read.csv(input$fMarkerExpr$datapath)
        df_global <<- df
        
        annotationlist <<- list("Unassigned")
        
        cell_freq <<- read.csv(input$cluster_freq$datapath)
        
        labels_row <-
          paste0(rownames(df), " (", cell_freq$clustering_prop , "%)")
        
        marker_names <- rownames(df)
        
        set.seed(1234)
        
        my_umap <- umap(df)
        dr_umap <<- data.frame(
          u1 = my_umap$layout[, 1],
          u2 = my_umap$layout[, 2],
          my_umap$data,
          cluster_number = 1:length(my_umap$layout[, 1]),
          check.names = FALSE
        )
        
        df01 <<- df %>% normalize01()
        myDF <<- df %>% normalize01()
        
        df01Tree <<- df %>% normalize01()
        allMarkers <<- colnames(df)
        df01Tree$cell <<- "Unassigned"
        
        selectedMarkers <<- colnames(df)
        posPickerList <<- colnames(df)
        
        
        #### IS THIS DOING ANYTHING
        # ## load only at start to fill the picker list
        # updatePickerInput(
        #   session,
        #   inputId = "myPickerPos",
        #   label = "Select Positive Markers",
        #   choices = colnames(df01),
        #   # choices = NULL,
        #   options = list(
        #     `actions-box` = TRUE,
        #     size = 10,
        #     `selected-text-format` = "count > 3"
        #   )
        # )
        # updatePickerInput(
        #   session,
        #   inputId = "myPickerNeg",
        #   label = "Select Negative Markers",
        #   choices = colnames(df01),
        #   options = list(
        #     `actions-box` = TRUE,
        #     size = 10,
        #     `selected-text-format` = "count > 3"
        #   )
        # )
        
        print( getDefaultReactiveDomain())
        browser()
        
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
        
        session$userData$vars$th <- kmeansTH(df01)
        
        annotaionDF <<- data.frame("cell" = "unassigned",
                                   clusterSize = cell_freq$clustering_prop)
        at <<- reactiveValues(data = annotaionDF, dr_umap = dr_umap)
        
      # })
      # Observe MenutItems ----
      observeEvent(input$tabs, {
        
        if (exists("df01")) {
          if(input$tabs=="thresholds") { ## Thresholds ----
            
            if (is.null(session$userData$vars$th) & !is.null(df01)) {
              session$userData$vars$th <- kmeansTH(df01)
            }
          }
          else if(input$tabs=="treeannotation") { # Tree annotation ----
            
            updatePickerInput(
              session,
              inputId = "parentPicker",
              choices = annotationlist
            )
            updatePickerInput(
              session,
              inputId = "updateNodePicker",
              choices = annotationlist,
              selected = ""
            )
            plotTree()
          }
        }
      })
      
      
      
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
      
      
      
      # Create new node ----
      observeEvent(input$createNodeBtn, {

        if (is.null(input$treePickerPos) & is.null(input$treePickerNeg)) {
          print("no selection")
          showModal(modalDialog(
            title = "No Marker Selection",
            "Select positive and / or negative markers!",
            easyClose = TRUE,
            footer = NULL,
            return()
          ))
        } else if(input$newNode == "") {
          showModal(modalDialog(
            title = "Phenotype Name",
            "Set a Name for this Phenotype!",
            easyClose = TRUE,
            footer = NULL,
            return()
          ))
        }
        else {

          name <- input$newNode

          # make sure name is not yet taken
          if(name %in% session$userData$vars$graph$nodes$label) {
            showModal(modalDialog(
              title = "Naming Error!",
              "The Name for this Phenotype is already taken!",
              easyClose = TRUE,
              footer = NULL
            ))

            return()
          }

          # browser()
          parent <- input$parentPicker

          # receive the parent settings, resp. parent hm
          # filter hm by parent cell name
          tmp <- df01Tree[df01Tree$cell == parent, ]
          tmp <- filterHM(tmp,input$treePickerPos, input$treePickerNeg, session$userData$vars$th)

          # Make sure the selection is not empty!
          if(dim(tmp)[1] == 0) {

            showModal(modalDialog(
              title = "No Result",
              "The Selection for this Phenotype is empty!",
              easyClose = TRUE,
              footer = NULL
            ))
            return()

          } else {

            annotationlist <<- append(annotationlist, name)

            updatePickerInput(
              session,
              inputId = "parentPicker",
              choices = annotationlist
            )
            updatePickerInput(
              session,
              inputId = "updateNodePicker",
              choices = annotationlist,
              selected = ""
            )

            posmarker <- list(input$treePickerPos)
            negmarker <- list(input$treePickerNeg)

            session$userData$vars$graph <- add_node(session$userData$vars$graph,parent,name,posmarker,negmarker,"blue")

            df01Tree[rownames(tmp),]$cell <<- name

            
            
            updateCheckboxGroupButtons(
              session,
              inputId = "treePickerPos",
              choices = colnames(df01),
              selected = NULL,
              disabledChoices = input$treePickerPos

            )
            updateCheckboxGroupButtons(
              session,
              inputId = "treePickerNeg",
              choices = colnames(df01),
              selected = NULL,
              disabledChoices = input$treePickerNeg

            )
            updatePickerInput(
              session,
              inputId = "parentPicker",
              selected = name
            )

            updateClusterLabels(tmp)

            plotTree()

          }

        }
      })
      
      
      

      # Observe node update picker ----
      observeEvent(input$updateNodePicker, {

        browser()

        node <- session$userData$vars$graph$nodes %>% filter(label == input$updateNodePicker)
        myid <- node$id

        filterPosMarkers <- unlist(node$pm)
        filterNegMarkers <- unlist(node$nm)

        updatePickerInput(
          session,
          inputId = "updatePickerPos",
          label = "Select Positive Markers",
          choices = colnames(df01),
          selected = filterPosMarkers,
          options = list(
            `actions-box` = TRUE,
            size = 10,
            `selected-text-format` = "count > 3"
          )
        )
        updatePickerInput(
          session,
          inputId = "updatePickerNeg",
          label = "Select Negative Markers",
          choices = colnames(df01),
          selected = filterNegMarkers,
          options = list(
            `actions-box` = TRUE,
            size = 10,
            `selected-text-format` = "count > 3"
          )
        )
        updateTextInput(
          session,
          inputId = "renameNode",
          label = NULL,
          value = node$label)
      })
      
     
      

      
      
      importNodes <- function(id, label, pm, nm, parent_id, color) {

        if(id == parent_id) {return()}

        parent_row <- session$userData$vars$graph$nodes %>% filter(id == parent_id)
        parent_name <- parent_row$label

        add_node(session$userData$vars$graph, parent_name, label, color)
      }
      
      
      # Export Annotation Tree ----
      output$exportAnnotationBtn <- downloadHandler(

        filename = function(){
          paste("cycadas_annotation_data_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file) {
          temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
          dir.create(temp_directory)

          export_df_nodes <- session$userData$vars$graph$nodes
          # Use apply() to concatenate the pm column for each row of the nodes dataframe
          pm_concatenated <- apply(export_df_nodes, 1, function(row) {

            pm_vector <- unlist(row["pm"])
            pm_vector <- paste(pm_vector, collapse = "|")
          })

          nm_concatenated <- apply(export_df_nodes, 1, function(row) {
            nm_vector <- unlist(row["nm"])
            nm_vector <- paste(nm_vector, collapse = "|")
          })

          # Add the concatenated pm column as a new column to the nodes dataframe
          export_df_nodes$pm <- pm_concatenated
          export_df_nodes$nm <- nm_concatenated

          download_list <- list(annTable = df01Tree,
                                nodesTable = export_df_nodes,
                                edgesTable = session$userData$vars$graph$edges)
          # download_list <- list(annTable = df01Tree["cell"],
          #                       nodesTable = export_df_nodes,
          #                       edgesTable = session$userData$vars$graph$edges)

          download_list %>%
            imap(function(x,y){
              if(!is.null(x)){
                file_name <- glue("{y}_data.csv")
                readr::write_csv(x, file.path(temp_directory, file_name))
              }
            })

          zip::zip(
            zipfile = file,
            files = dir(temp_directory),
            root = temp_directory
          )
        },
        contentType = "application/zip"
      )
      
      
    })}




# 

