#' @export
cycadas <- function() {

  reactVals <- reactiveValues(th = NULL,
                              myTH = NULL, # user selected threshold
                              md = NULL,
                              counts_table = NULL,
                              DA_result_table = NULL,
                              DA_interactive_table = NULL,
                              myNode = "", # selected node for interactive DA
                              graph = NULL,
                              hm = NULL,
                              annotationlist = NULL,
                              sce = NULL,
                              metaClustLevel = NULL
                              )
  
  # Global static values
  lineage_marker <<- character(0)
  df_expr <<- NULL
  cell_freq <<- NULL
  # dr_umap <<- NULL
  dr_umap <- reactiveVal(NULL)   # starts empty
  rFreqs <- reactiveVal(NULL)
  init_app <<- TRUE
  # 2 column data.frame containing old_cluster and new_cluster IDs
  mergeTableCatalyst <<- NULL

  # Server function -------------------------
  server = function(input, output, session) {
    
    shinyjs::disable("metadiv")

    # Check if the package is available
    package_available <- requireNamespace("CATALYST", quietly = TRUE)
    
    # Render the conditional UI elements
    output$test <- renderUI({
      if (package_available) {
        library(CATALYST)
        library(SingleCellExperiment)
        fileInput("sce","Upload RDS",
                  placeholder = "Choose RDS File",
                  multiple = FALSE,
                  accept = c(".rds"))
      } else {
        p("The CATALYST package is not available.")
      }
    })
    
    # Load Expression Data and UMAP -------------------------------------------
    loadExprData <- function(path_expr, path_freq) {
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "loading Data...", value = 0)
      
      ## Load median expression and cell frequencies
      ## add the 0-1 scaled columns to the DF
      df_expr <<- read.csv(path_expr)
      cell_freq <<- read.csv(path_freq)
      
      lineage_marker <<- colnames(df_expr)
      lineage_marker_raw <<- paste0(lineage_marker, "_raw")
      
      df_expr <<- createExpressionDF(df_expr, cell_freq)
      reactVals$graph <- initTree()
      
      progress$set(message = "loading Data Cluster Expression Demo Data...", value = 0.2)
      
      reactVals$graph <- initTree()
      
      set.seed(1234)
      progress$set(message = "Building the UMAP...", value = 0.3)
      umap_result <- buildUMAP(df_expr[, lineage_marker_raw]) 
      dr_umap(umap_result)
      
    }
    
    
    # both files required
    gs_ready <- reactive({
      has_file(input$fMarkerExpr) && has_file(input$cluster_freq)
    })
    
    # keep button state in sync
    observe({
      if (gs_ready()) shinyjs::enable("btnImportExpr") else shinyjs::disable("btnImportExpr")
    })
    
    # live status lights
    output$gs_status <- renderUI({
      tagList(
        status_item(has_file(input$fMarkerExpr), "Marker Expressions (CSV)"),
        status_item(has_file(input$cluster_freq), "Cluster Frequencies (CSV)"),
        tags$div(style="margin-top:8px;",
                 if (gs_ready())
                   tags$span("Both files present — ready to import.", style="color:#28a745;")
                 else
                   tags$span("Please select both files to enable Import.", style="color:#dc3545;")
        )
      )
    })
    
    # Import click
    observeEvent(input$btnImportExpr, {
      req(gs_ready())  # safety
      
      ok <- try({
        loadExprData(
          path_expr = input$fMarkerExpr$datapath,
          path_freq = input$cluster_freq$datapath
        )
      }, silent = TRUE)
      
      if (isTRUE(ok)) {
        initExprData()
        showNotification("GigaSOM / FlowSOM data imported.", type = "message")
        # optional: disable until new files are chosen
        shinyjs::disable("btnImportExpr")
      } else {
        showNotification("Import failed. Please check your CSVs.", type = "error", duration = 8)
        # optional: clear the inputs so the user re-selects
        shinyjs::reset("gs_upload_form")
      }
    })
    
    # This initializes the checkboxes and thresholds at start of a new 
    # Project or for the non-annotated demo data
    # Init the checkboxes and generate the Thresholds -------------------------
    initExprData <- function() {
      
      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)
      
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = lineage_marker,
        selected = NULL
      )
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = lineage_marker,
        selected = NULL
      )
      
      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)
      
      reactVals$th <- kmeansTH(df_expr[, lineage_marker])
      reactVals$annotationlist <- c("Unassigned")
    }
    
    # Function: Update cluster labels -----------------------------------------
    updateClusterLabels <- function(mydf, my_node_id=0) {

      children <- all_my_children(reactVals$graph, my_node_id)
      
      if (is.null(children)) {
        children <- c()
      } 
      
      all_cluster_ids <- as.numeric(c(rownames(mydf), children))
      totalSum <- sum(cell_freq[all_cluster_ids,]$clustering_prop)

      mysum <- sum(cell_freq[rownames(mydf),]$clustering_prop)

      output$progressBox <- renderText({
        paste0(round(mysum, 3), "% in selection")
      })
      output$totalSelection <- renderText({
        paste0(round(totalSum, 3), "% in total")
      })
      output$progressBox2 <- renderText({
        paste0(dim(mydf)[1], " Cluster in selection")
      })
    }

    # Function: Plot the annotation tree --------------------------------------
    plotTree <- function() {
      # before plotting the tree, update its properties
      for(i in 1:nrow(reactVals$graph$nodes)) {
        
        l <- reactVals$graph$nodes$label[i]
        # get the number of clusters with that label
        nLabel <- sum(df_expr$cell == l)

        if (nLabel == 0) {
          reactVals$graph$nodes[i,]$color <- "grey"
        } else {
          reactVals$graph$nodes[i,]$color <- "blue"
        }
      }

      output$mynetworkid <- renderVisNetwork({
        visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
          visEdges(arrows = "from") %>%
          visHierarchicalLayout() %>%
          visEvents(select = "function(nodes) {
                Shiny.onInputChange('parent_node_id', nodes.nodes);
                ;}")
      })
    }
    
    
    # Load CATALYST data ----
    observeEvent(input$sce, {

      req(input$sce)
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "loading Data ...", value = 0.2)
      
      reactVals$sce <- readRDS(input$sce$datapath)
      
      reactVals$metaClustLevel <- colnames(reactVals$sce@metadata$cluster_codes)[1]

      # load only som codes as default
      df_expr_tmp <- as.data.frame(reactVals$sce@metadata$SOM_codes)
      
      rd <- rowData(reactVals$sce)
      
      lineage_marker <<- rd$marker_name[rd$marker_class == "type"]
      lineage_marker_raw <<- paste0(lineage_marker, "_raw")
      
      # cell_freq_tmp <- as.data.frame(cluster_ids(c_id = cluster_ids(reactVals$sce, "meta20")))
      cell_freq_tmp <- as.data.frame(table(c_id = cluster_ids(reactVals$sce)))
      
      cell_freq_tmp$clustering_prop <- cell_freq_tmp$Freq / sum(cell_freq_tmp$Freq)
      cell_freq_tmp$Freq <- NULL
      colnames(cell_freq_tmp) <- c("cluster", "clustering_prop")
      cell_freq <<- cell_freq_tmp
      
      
      df_expr <<- createExpressionDF(df_expr_tmp, cell_freq_tmp)
      reactVals$graph <- initTree()
      
      set.seed(1234)
      progress$set(message = "Building the UMAP...", value = 0.3)
      umap_result <- buildUMAP(df_expr[, lineage_marker_raw]) 
      dr_umap(umap_result)

      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)

      initExprData()
      
      mc_levels <- colnames(reactVals$sce@metadata$cluster_codes)
      #remove 2,3, and 4
      mc_levels <- mc_levels[-(2:5)]
      
      updatePickerInput(
        session,
        inputId = "metaLevel",
        choices = mc_levels,
        selected = reactVals$metaClustLevel
      )
      
      shinyjs::enable("metadiv")
    })
    
    # Observe Metalevel selection -------------------------------------------
    observeEvent(input$metaLevel, {
      
      req(input$sce)
      
      if (init_app == TRUE) {
        
        init_app <<- FALSE
        return()
      } else {

        reactVals$metaClustLevel=input$metaLevel
        
        # select the metalvel and merge the clusters 
        somExpr <- as.data.frame(reactVals$sce@metadata$SOM_codes)
        
        somExpr$metaLevels <- reactVals$sce@metadata$cluster_codes[reactVals$metaClustLevel]
        
        somExpr <- somExpr %>%
          group_by(metaLevels) %>%
          summarise_if(is.numeric, mean)
        
        somExpr <- as.data.frame(somExpr[, lineage_marker])

        cell_freq_tmp <- as.data.frame(table(c_id = cluster_ids(reactVals$sce, reactVals$metaClustLevel)))
        cell_freq_tmp$clustering_prop <- cell_freq_tmp$Freq / sum(cell_freq_tmp$Freq)
        cell_freq_tmp$Freq <- NULL
        colnames(cell_freq_tmp) <- c("cluster", "clustering_prop")
        cell_freq <<- cell_freq_tmp
        
        df_expr <<- createExpressionDF(somExpr, cell_freq_tmp)
        reactVals$graph <- initTree()
        
        set.seed(1234)

        if (nrow(df_expr) > 20) {
          umap_result <- buildUMAP(df_expr[, lineage_marker_raw])   
          dr_umap(umap_result)
        } else {
          dr_umap(NULL)
        }
        
        updateSelectInput(session, "markerSelect", "Select:", lineage_marker)
        
        initExprData()
        
        reactVals$hm <- df_expr[, lineage_marker]
      }
    })
    
    # Observe MenutItems ------------------------------------------------------
    observeEvent(input$tabs, {

      if (exists("df_expr")) {
        if(input$tabs=="thresholds") { 

          if (is.null(reactVals$th) & !is.null(df_expr)) {
            reactVals$th <- kmeansTH(df_expr)
          }
        }
        else if(input$tabs=="treeannotation") { 

          if (!is.null(reactVals$graph)) {
            updatePickerInput(
              session,
              inputId = "parentPicker",
              choices = reactVals$annotationlist
            )
          }
          if (!is.null(reactVals$graph))
            plotTree()  
        }
      }
    })

    # Observe Tree picker -----------------------------------------------------
    observeEvent(input$treePickerPos, {

      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = lineage_marker,
        selected = input$treePickerNeg,
        disabledChoices = input$treePickerPos
      )
    })
    observeEvent(input$treePickerNeg, {

      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = lineage_marker,
        selected = input$treePickerPos,
        disabledChoices = input$treePickerNeg
      )
    })

    # Create new node ---------------------------------------------------------
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
        if(name %in% reactVals$graph$nodes$label) {
          showModal(modalDialog(
            title = "Naming Error!",
            "The Name for this Phenotype is already taken!",
            easyClose = TRUE,
            footer = NULL
          ))
          return()
        }

        parent <- input$parentPicker
        # receive the parent settings, resp. parent hm
        # filter hm by parent cell name
        tmp <- df_expr[df_expr$cell == parent, ]
        tmp <- filterHM(tmp,input$treePickerPos, input$treePickerNeg, reactVals$th)

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

          reactVals$annotationlist <- append(reactVals$annotationlist, name)

          updatePickerInput(
            session,
            inputId = "parentPicker",
            choices = reactVals$annotationlist
          )
          updatePickerInput(
            session,
            inputId = "updateNodePicker",
            choices = reactVals$annotationlist,
            selected = ""
          )

          posmarker <- list(input$treePickerPos)
          negmarker <- list(input$treePickerNeg)

          reactVals$graph <- add_node(reactVals$graph,parent,name,posmarker,negmarker,"blue")

          df_expr[rownames(tmp),]$cell <<- name

          updateCheckboxGroupButtons(
            session,
            inputId = "treePickerPos",
            choices = lineage_marker,
            selected = NULL,
            disabledChoices = input$treePickerPos

          )
          updateCheckboxGroupButtons(
            session,
            inputId = "treePickerNeg",
            choices = lineage_marker,
            selected = NULL,
            disabledChoices = input$treePickerNeg

          )
          updatePickerInput(
            session,
            inputId = "parentPicker",
            selected = name
          )

          # browser()
          updateClusterLabels(tmp, parent)

          plotTree()
        }
      }
    })
    
    # Merge CATALYST function -------------------------------------------------
    observeEvent(input$mergeCatalystBtn, {

      # browser()
      mergeTableCatalyst <<- data.frame(old_cluster = as.numeric(rownames(df_expr)),
                                       new_cluster = df_expr$cell)

      reactVals$sce <- mergeClusters(reactVals$sce, k = reactVals$metaClustLevel, table = mergeTableCatalyst, 
                                     id = paste0("merging_", reactVals$metaClustLevel))

    })
    
    # Export CATALYST SingleCellObject ----------------------------------------
    output$exportCatalystBtn <- downloadHandler(
      filename = function() {
        paste("Annotated_sce_", Sys.Date(), ".rds", sep="")
      },
      content = function(file) {
        saveRDS(reactVals$sce, file)
      },
      contentType = "application/octet-stream"
    )
    
    # Annotation heatmap plot -------------------------------------------------
    output$hm_tree <- renderPlot({

      if (is.null(reactVals$hm)) {
        plot(1, type = "n", main = "No Data Available")
      }
      else if (nrow(reactVals$hm) == 0) {
        plot(1, type = "n", main = "Node is empty!")
      }
      else if (nrow(reactVals$hm) < 2) {
        pheatmap(reactVals$hm, cluster_cols = F, cluster_rows = F)
      } else {
        message(nrow(reactVals$hm))
        pheatmap(reactVals$hm, cluster_cols = F)
      }
    })

    # Observe Parent Node selection -------------------------------------------
    observeEvent(input$parentPicker, {
      
      nodeName <- input$parentPicker

      node <- reactVals$graph$nodes %>% dplyr::filter(label == nodeName)
      parent <- input$parentPicker
      
      nodeID <- reactVals$graph$nodes$id[reactVals$graph$nodes$label == nodeName]

      filterPosMarkers <- unlist(node$pm)
      filterNegMarkers <- unlist(node$nm)
      
      # receive the parent settings, resp. parent hm
      # filter hm by parent cell name
      tmp <- df_expr[df_expr$cell == parent, lineage_marker]
      tmp <- filterHM(tmp, input$treePickerPos, input$treePickerNeg, reactVals$th)

      updateClusterLabels(tmp, nodeID)
      
      reactVals$hm <- tmp
      
      updatePickerInput(
        session,
        inputId = "parentPicker",
        choices = reactVals$annotationlist,
        selected = parent
      )
      
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = lineage_marker,
        selected = NULL,
        disabledChoices = filterPosMarkers
      )
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = lineage_marker,
        selected = NULL,
        disabledChoices = filterNegMarkers
      )
      
      # update umap plot for Tree
      ClusterSelection <- filterColor(df_expr,tmp)
      
      #---- UMAP Annotation Tree
      output$umap_tree <-
        renderPlot(
          cbind(dr_umap(),ClusterSelection) %>%
            dplyr::mutate(ClusterSelection=as.factor(ClusterSelection)) %>%
            dplyr::mutate(ClusterSelection=fct_relevel(ClusterSelection,"other clusters","selected phenotype")) %>%
            arrange(desc(ClusterSelection)) %>%
            ggplot( aes(x = u1, y = u2, color = ClusterSelection)) +
            geom_point(size = 1.0) +
            # ggpubr::theme_pubr() +
            theme_classic() +
            theme(legend.text = element_text(size = 12),
                  legend.title = element_text(size = 20),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 20)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
        )
    })
    
    # Observe Interactive Parent Node selection -------------------------------
    observeEvent(input$parent_node_id, {

      node <- reactVals$graph$nodes %>% dplyr::filter(id == input$parent_node_id)
      parent <- node$label
      
      nodeID <- reactVals$graph$nodes$id[reactVals$graph$nodes$label == parent]
      
      filterPosMarkers <- unlist(node$pm)
      filterNegMarkers <- unlist(node$nm)
      
      # receive the parent settings, resp. parent hm
      # filter hm by parent cell name
      tmp <- df_expr[df_expr$cell == parent, lineage_marker]
      tmp <- filterHM(tmp, input$treePickerPos, input$treePickerNeg, reactVals$th)

      updateClusterLabels(tmp, nodeID)
      
      reactVals$hm <- tmp
      
      updatePickerInput(
        session,
        inputId = "parentPicker",
        choices = reactVals$annotationlist,
        selected = parent
      )
      
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = lineage_marker,
        selected = NULL,
        disabledChoices = filterPosMarkers
      )
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = lineage_marker,
        selected = NULL,
        disabledChoices = filterNegMarkers
      )
      
      # update umap plot for Tree
      ClusterSelection <- filterColor(df_expr,tmp)
      
      #---- UMAP Annotation Tree
      output$umap_tree <-
        renderPlot(
          cbind(dr_umap(),ClusterSelection) %>%
            dplyr::mutate(ClusterSelection=as.factor(ClusterSelection)) %>%
            dplyr::mutate(ClusterSelection=fct_relevel(ClusterSelection,"other clusters","selected phenotype")) %>%
            arrange(desc(ClusterSelection)) %>%
            ggplot( aes(x = u1, y = u2, color = ClusterSelection)) +
            geom_point(size = 1.0) +
            # ggpubr::theme_pubr() +
            theme_classic() +
            theme(legend.text = element_text(size = 12),
                  legend.title = element_text(size = 20),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 20)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
        )
    })

    # Observe node update picker ----------------------------------------------
    ## Function currently not ins use!
    observeEvent(input$updateNodePicker, {

      node <- reactVals$graph$nodes %>% dplyr::filter(label == input$updateNodePicker)
      myid <- node$id

      filterPosMarkers <- unlist(node$pm)
      filterNegMarkers <- unlist(node$nm)

      updatePickerInput(
        session,
        inputId = "updatePickerPos",
        label = "Select Positive Markers",
        choices = lineage_marker,
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
        choices = lineage_marker,
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

    # Update Node -------------------------------------------------------------
    ## Function currently not ins use!
    observeEvent(input$updateNodeBtn, {

      node <- reactVals$graph$nodes %>% dplyr::filter(label == input$parentPicker)
      myid <- node$id

      reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker,]$pm <-
              list(c(reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker,]$pm[[1]], input$treePickerPos))

      reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker,]$nm <-
        list(c(reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker,]$nm[[1]], input$treePickerNeg))

      reactVals <- rebuiltTree(reactVals)

      reactVals$hm <- reactVals[reactVals$cell == input$parentPicker, ]

      if (input$newNode != "") {

        old_name <- input$parentPicker
        new_name <- input$newNode

        reactVals$graph$nodes[reactVals$graph$nodes$label == old_name,]$label <- new_name

        reactVals$annotationlist[reactVals$annotationlist == old_name] <<- new_name
        updatePickerInput(
          session,
          inputId = "parentPicker",
          selected = new_name,
          choices = reactVals$annotationlist
        )

        # replace name in heatmap
        reactVals$cell[reactVals$cell == old_name] <<- new_name
        reactVals$hm <- reactVals[reactVals$cell == new_name, ]

        updateTextInput(
          session,
          inputId = "newNode",
          label = NULL,
          value = "",
          placeholder = NULL
        )
      }

      plotTree()
    })

    # Delete Node -------------------------------------------------------------
    observeEvent(input$deleteNodeBtn, {

      node <- reactVals$graph$nodes[reactVals$graph$nodes$label == input$parentPicker, ]
      # check and make sure that this node is a leaf node
      if(TRUE %in% (reactVals$graph$edges$to == node$id)) {

        showModal(modalDialog(
          title = "Cannot delete Node",
          "The selected Node is not a Leaf Node!",
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
      else {

        # before deleting the node, we must assign the cluster annotation
        # name for that node back to its parent
        parent_id <- reactVals$graph$edges$to[reactVals$graph$edges$from == node$id]
        parent_label <- reactVals$graph$nodes$label[reactVals$graph$nodes$id == parent_id]

        # In case there is a empty node with no the filtered clusters
        # we do not assign any labeling
        if( dim(df_expr[df_expr$cell == node$label,])[1] >0 ) {
          
          df_expr[df_expr$cell == node$label,]$cell <<- parent_label
        }

        # Remove label form annotation list
        tmplist <- reactVals$annotationlist 
        tmplist <- tmplist[!tmplist == node$label]
        reactVals$annotationlist <- tmplist

        updatePickerInput(
          session,
          inputId = "parentPicker",
          selected = NULL,
          choices = reactVals$annotationlist
        )
        updatePickerInput(
          session,
          inputId = "updateNodePicker",
          selected = NULL,
          choices = reactVals$annotationlist
        )
        updateTextInput(
          session,
          inputId = "renameNode",
          label = NULL,
          value = ""
        )
        updateTextInput(
          session,
          inputId = "newNode",
          label = NULL,
          value = ""
        )

        reactVals$graph <- delete_leaf_node(reactVals$graph, node$id)
        plotTree()
      }
    })

    # Load Annotation Tree ----------------------------------------------------
    observeEvent(input$btnImportTree, {

      req(input$fNodes)
      req(input$fEdges)
      req(input$fMarkerExpr)
      req(input$cluster_freq)

      df_nodes <- read.csv(input$fNodes$datapath)
      df_edges <- read.csv(input$fEdges$datapath)
      
      reactVals$graph <- getGraphFromLoad(df_nodes, df_edges)
      
      df_expr$cell <<- rebuiltTree(reactVals$graph, df_expr, reactVals$th)
      
      reactVals$annotationlist <- as.list(df_nodes$label)
    })


    # TODO: find a method to get a good resolution output
    observeEvent(input$exportTreeGraphics, {

    })

    # Export Annotation Tree --------------------------------------------------
    output$exportAnnotationBtn <- downloadHandler(

      filename = function(){
        paste("cycadas_annotation_data_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file) {

        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)

        export_df_nodes <- reactVals$graph$nodes
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

        download_list <- list(annTable = df_expr,
                              nodesTable = export_df_nodes,
                              edgesTable = reactVals$graph$edges)

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

    # Output Thresholds Tab ---------------------------------------------------

    output$table = renderDataTable(
      reactVals$th[, c("cell", "threshold", "bi_mod")],
      editable = F,
      extensions = c('Buttons', 'Scroller'),
      selection = 'single',
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

    # Observe marker selection in Thresholds Table ---------------------------
    # from marker, plot the expression in scatter-plot or histogram
    observeEvent(input$table_rows_selected, {

      selRow <- reactVals$th[input$table_rows_selected,]
      myTH <- reactVals$th[input$table_rows_selected, "threshold"]
      myColor <- reactVals$th[input$table_rows_selected, "color"]

      marker_expr <- as.data.frame(df_expr[, selRow$cell])
      marker_expr <- cbind(marker_expr, rnorm(1:nrow(marker_expr)))
      
      PlottingThreshold(marker_expr, myTH, myColor)
    })
    
    # Observe Threshold Method Selection --------------------------------------
    observeEvent(input$th_method, {

      if (!is.null(reactVals$th)) {
        # TODO: change naming: not only k-means
        reactVals$th <- updateTH(df_expr, reactVals$th, th_mode=input$th_method)
      }
    })
    
    # Observe click scatter-plot -----------------------------------------------
    observeEvent(input$plot_click, {

      req(input$table_rows_selected)
      selRow <- reactVals$th[input$table_rows_selected,]

      reactVals$th[input$table_rows_selected, "threshold"] <- round(input$plot_click$x, 3)
      
      myTH <- reactVals$th[input$table_rows_selected, "threshold"]
      myColor <- reactVals$th[input$table_rows_selected, "color"]
      marker_expr <- as.data.frame(df_expr[, selRow$cell])
      marker_expr <- cbind(marker_expr, rnorm(1:nrow(marker_expr)))
      
      PlottingThreshold(marker_expr, myTH, myColor)

      df_expr$cell <<- rebuiltTree(reactVals$graph, df_expr, reactVals$th)
    })

    # Plotting Threshold  function --------------------------------------------
    PlottingThreshold <- function(me, myTH, myCol){

      # Scatter-plot
      output$plot <- renderPlot({
        set.seed(1)

        ggplot(me, aes_string(x=me[,1], y=me[,2])) +
          geom_point(size=1) +
          theme_classic() +
          # ggpubr::theme_pubr() +
          theme(axis.title.y = element_blank(),
                axis.ticks.y  = element_blank(),
                axis.text.y = element_blank(),
                axis.text.x = element_text(size = 12),
                axis.title.x = element_text(size = 18),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank()) +
          labs(x = "Scale 0 to 1") +
          geom_vline(xintercept = myTH, 
                     linetype="dotted",
                     color = myCol, 
                     size=1.5)

      })
      
      # Histogram
      output$plot2 <- renderPlot({

        ggplot(me, aes_string(x = me[, 1])) +
          geom_histogram(bins = 80) +
          labs(x = "Scale 0 to 1") +
          # ggpubr::theme_pubr() +
          theme_classic() +
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 18)) +
          geom_vline(
            xintercept = myTH,
            linetype = "dotted",
            color = myCol,
            size = 1.5
          )
      })
    }

    # Save the DA Table -------------------------------------------------------
    output$exportDA <- downloadHandler(
      filename = function() {
        paste("DA_Table_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(reactVals$DA_result_table, file)
      }
    )
    
    # Get the merged phenotypes as proportion table  --------------------------
    get_merged_prop_table <- function() {
      
      req(reactVals$counts_table)
      
      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df_expr$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)
      
      # change the name of nodes:
      # if node has children add "_remaining"
      my_list <- lapply(countsTable$cell, function(x) {
        # get the id
        nid <- reactVals$graph$nodes$id[reactVals$graph$nodes$label == x]
        # check if the id has children
        if(nid %in% reactVals$graph$edges$to) {
          return (x <- paste0(x, "_remaining"))
        }
        else {
          return(x)
        }
      })
      
      countsTable$cell <- unlist(my_list)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL
      
      props_table <- t(t(countsTable) / colSums(countsTable)) * 100
      
      return(props_table)
    }
    
    # Save the Merged Proportion Table ----------------------------------------
    output$exportProp <- downloadHandler(
      
      filename = function() {
        paste("Merged_Proportions_Table_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(get_merged_prop_table(), file)
      }
    )
    
    # Load the marker threshold file ---------------------------------------
    observeEvent(input$fTH,{

      th <<- read.csv(input$fTH$datapath)
      th$X <- NULL
      th$color <- "blue"

      th[th$bi_mod < 0.555, "color"] <- "red"

      rownames(th) <- th$cell
      reactVals$th<- th
      
      # re-calculate the tree after threshold upload
      req(input$fMarkerExpr)
      df_expr$cell <<- rebuiltTree(reactVals$graph, df_expr, reactVals$th)
    })

    # Load the metadata ----------------------------------------------------
    observeEvent(input$metadata,{
      md <<- read.csv(input$metadata$datapath)
      md$X <- NULL
      reactVals$md<- md
    })

    # Load the counts talbe -------------------------------------------------
    observeEvent(input$counts_table,{
      ct <<- read.csv(input$counts_table$datapath)
      ct$X <- NULL
      reactVals$counts_table <- ct
    })

    # UMAP interactive Tab ----------------------------------------------------
    output$umap2 <- renderPlot({

      req(dr_umap())

      ggplot(dr_umap(), aes(x = u1, y = u2)) +
        geom_point(size = 1.0) +
        theme_classic() +
        theme(legend.text=element_text(size=8),
              legend.title = element_text(size = 20),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 20)) +
        guides(color = guide_legend(override.aes = list(size = 4)))        
    })

    # Heatmap from UMAP selection ---------------------------------------------
    output$hm2 <- renderPlot({
      
      req(dr_umap())
      my_new_hm <- round(brushedPoints(dr_umap(), input$umap2_brush), 2)

      if (dim(my_new_hm)[1] > 0) {

        pheatmap(my_new_hm %>% select(-c("cluster_number", "u1", "u2")), cluster_cols = F)
      } else {
        ggplot() + theme_void() + ggtitle("Select area on Umap to plot Heatmap")
      }
    })

    # Data Table from UMAP selection ------------------------------------------

    output$umap_data <- DT::renderDT(server = FALSE, {
      
      req(dr_umap)

      # 2) No selection (or empty selection)
      brushed <- brushedPoints(dr_umap(), input$umap2_brush)
      if (is.null(brushed) || nrow(brushed) == 0) {
        return(DT::datatable(
          data.frame(Message = "No points selected"),
          rownames = FALSE,
          options = list(dom = 't', paging = FALSE)
        ))
      }
      
      # 3) Show selected points
      DT::datatable(
        round(brushed, 2),
        filter = 'top',
        extensions = 'Buttons',
        options = list(
          scrollY = 600,
          scrollX = TRUE,
          dom = '<"float-left"l><"float-right"f>rt<"row"<"col-sm-4"B><"col-sm-4"i><"col-sm-4"p>>',
          lengthMenu = list(c(10, 25, 50, -1), c('10', '25', '50', 'All')),
          scrollCollapse = TRUE,
          lengthChange = TRUE,
          widthChange = TRUE,
          rownames = TRUE
        )
      )
    })
    

    # UMAP Marker Expression Tab ----------------------------------------------
    output$umap3 <-
      renderPlot({

      req(dr_umap())

        ggplot(dr_umap(), aes_string(
          x = "u1",
          y = "u2",
          color = paste0(input$markerSelect, "_raw")
        )) +
          geom_point(size = 1.0) +
          theme_bw() +
          theme(legend.text = element_text(size = 14),
                legend.title = element_text(size = 20),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 20)) +
          scale_color_gradientn(input$markerSelect,
                                colours = colorRampPalette(rev(brewer.pal(
                                  n = 11, name = "Spectral"
                                )))(50))
      })
    
    

    # Differential Abundance Tab ----------------------------------------------
    output$md_table <-
      renderTable(
        reactVals$md[1:5, ]
      )

    output$counts_table <-
      renderTable(
        reactVals$counts_table[1:5, 1:5]
      )

    output$DA_result_table <-
      renderTable(
        reactVals$DA_result_table
      )
    
    output$DA_interactive_table <-
      renderTable(
        reactVals$DA_interactive_table
      )
    
    # Do the differential abundance ------------------------------------------
    observeEvent(input$doDA, {
      req(reactVals$counts_table, reactVals$md)
      md <- reactVals$md
      
      # ---
      if (!all(c("sample_id", "condition") %in% colnames(md))) {
        showNotification("Metadata must contain 'sample_id' and 'condition'.", type = "error")
        return()
      }
      md$sample_id <- as.character(md$sample_id)
      
      # --- 1) Build counts table aggregated by cell label ----
      countsTable <- tryCatch({
        ct <- reactVals$counts_table
        # needs a 'cell' vector aligned with rows of df_expr
        if (is.null(df_expr) || is.null(df_expr$cell)) {
          showNotification("Missing 'cell' labels in expression table.", type = "error")
          return(NULL)
        }
        if (nrow(ct) != nrow(df_expr)) {
          # If your counts table is per-cluster (rows=clusters), this should match.
          # If not, adapt aggregation logic here.
          showNotification("Row mismatch between counts table and df_expr. Check inputs.", type = "error")
          return(NULL)
        }
        ct$cell <- df_expr$cell
        # aggregate by cell
        agg <- aggregate(. ~ cell, ct, sum, na.rm = TRUE)
        rownames(agg) <- agg$cell
        agg$cell <- NULL
        agg
      }, error = function(e) {
        showNotification(paste("Failed to prepare counts table:", e$message), type = "error")
        NULL
      })
      req(!is.null(countsTable))
      
      # --- 2) Convert to proportions (%) with zero-sum guard ----
      cs <- colSums(countsTable, na.rm = TRUE)
      if (any(cs == 0)) {
        zero_cols <- names(cs)[cs == 0]
        showNotification(sprintf("Some samples have zero total counts and will be dropped: %s",
                                 paste(zero_cols, collapse = ", ")), type = "warning")
      }
      keep_cols <- names(cs)[cs > 0]
      if (length(keep_cols) == 0) {
        showNotification("All samples have zero counts; cannot compute proportions.", type = "error")
        return()
      }
      countsTable <- countsTable[, keep_cols, drop = FALSE]
      props_table <- sweep(countsTable, 2, colSums(countsTable, na.rm = TRUE), "/") * 100
      
      # --- 3) Match props_table columns to metadata sample_id safely ----
      mm <- match(colnames(props_table), md$sample_id)
      if (anyNA(mm)) {
        missing <- colnames(props_table)[is.na(mm)]
        showNotification(sprintf("Samples missing in metadata and dropped: %s",
                                 paste(missing, collapse = ", ")), type = "warning")
      }
      keep <- which(!is.na(mm))
      if (length(keep) == 0) {
        showNotification("No overlap between props_table columns and metadata sample_id.", type = "error")
        return()
      }
      props_table <- props_table[, keep, drop = FALSE]
      tmp_cond <- droplevels(as.factor(md$condition[mm[keep]]))
      
      # --- 4) Validate groups for testing ----
      if (nlevels(tmp_cond) < 2) {
        showNotification("Need at least two conditions to run DA tests.", type = "error")
        return()
      }
      
      # --- 5) Run pairwise Wilcoxon per row (cluster/cell) ----
      DA_df <- data.frame(stringsAsFactors = FALSE)
      for (i in seq_len(nrow(props_table))) {
        y <- as.numeric(props_table[i, ])
        # Skip if all NA or constant
        if (all(!is.finite(y)) || length(unique(y[is.finite(y)])) < 2) next
        
        # Ensure each group has at least one finite value
        ok_groups <- tapply(y, tmp_cond, function(v) sum(is.finite(v)) > 0)
        if (any(!ok_groups)) next
        
        pw <- tryCatch(
          pairwise.wilcox.test(y, tmp_cond, p.adjust.method = input$correction_method),
          error = function(e) NULL
        )
        if (is.null(pw) || is.null(pw$p.value)) next
        
        # Melt p-value matrix (upper triangle); keep non-NA
        m <- reshape2::melt(pw$p.value, varnames = c("Cond1", "Cond2"), value.name = "p.value")
        m <- subset(m, is.finite(p.value))
        
        if (nrow(m) == 0) next
        
        m$Cell <- rownames(countsTable)[i]
        
        # --- 6) Optional: rename nodes with "_remaining" if node has children ----
        nm <- m$Cell
        if (!is.null(reactVals$graph) && !is.null(reactVals$graph$nodes) && !is.null(reactVals$graph$edges)) {
          nodes <- reactVals$graph$nodes
          edges <- reactVals$graph$edges
          if (all(c("id", "label") %in% names(nodes)) && all(c("from", "to") %in% names(edges))) {
            has_children <- function(label) {
              nid <- nodes$id[nodes$label == label]
              length(nid) > 0 && any(edges$to %in% nid)
            }
            nm <- vapply(m$Cell, function(x) if (has_children(x)) paste0(x, "_remaining") else x, character(1))
          }
        }
        m$Naming <- nm
        
        DA_df <- rbind(DA_df, m)
      }
      
      if (nrow(DA_df) == 0) {
        showNotification("No significant pairwise results computed (check inputs/conditions).", type = "warning")
      }
      
      # Standardize column order/names
      colnames(DA_df) <- c("Cond1", "Cond2", "p-value", "Cell", "Naming")
      reactVals$DA_result_table <- DA_df
      showNotification("Differential abundance testing completed.", type = "message")
    })
 
    # Render interactive Tree -------------------------------------------------
    output$interactiveTree <- renderVisNetwork({
      visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    })

    # Observe Interactive DA node selection -----------------------------------
    observeEvent(input$current_node_id, {
      # browser()
      reactVals$myNode <- input$current_node_id
      
      if (reactVals$myNode != 1) {
        doInteractiveDA()  
      }
      
    })

    output$selectedNode <- renderText({reactVals$myNode})
    
    doInteractiveDA <- function() {  

      # get selected "parent" first
      selected_labels <- reactVals$graph$nodes$label[reactVals$graph$nodes$id == reactVals$myNode]
      children <- all_my_children(reactVals$graph, reactVals$myNode)
      # attach children 
      if (!is.null(children) & !is_empty(children)) {
        
        row_ids <- reactVals$graph$nodes$id %in% children
        selected_labels <- c(reactVals$graph$nodes$label[row_ids], selected_labels)
      } 

      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df_expr$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL
      props_table <- t(t(countsTable) / colSums(countsTable)) * 100

      # check if any of the nodes is empty and therefore not in the props table
      # this cause the
      # make sure that sleecte labels are larger than on, otherwise it is
      # not a vector and treat it as single value
      if (length(selected_labels) > 1) {
        filter_vec <- selected_labels %in% rownames(props_table)
        props_table <- props_table[selected_labels[filter_vec],]
      }
      else {
        props_table <- t(as.data.frame(props_table[selected_labels,]))
      }

      ## subset the props table based on the selected node
      # props_table <- props_table[selected_labels, ]
      props_table <- t(as.data.frame(colSums(props_table)))
      mm <- match(colnames(props_table), md$sample_id)
      tmp_cond <- md$condition[mm]
      
      DA_df <- data.frame()

      foo <- pairwise.wilcox.test(props_table, tmp_cond, p.adjust.method="none")

      df <- subset(melt(foo$p.value), value!=0)
      df$cell <- reactVals$graph$nodes$label[reactVals$myNode]
      DA_df <- rbind(DA_df, as.data.frame(df))

      colnames(DA_df) <- c("Var1", "Var2", "p-value", "Cell")
      reactVals$DA_interactive_table <- DA_df

    }
    
    # Boxplot of interactive DA selection -------------------------------------
    output$boxplot <- renderPlot({

      # first check if a node is selected:
      if (reactVals$myNode == "") {
        
        plot(1, type = "n", main = "No Data Available")
      } else {

        # get selected "parent" first
        selected_labels <- reactVals$graph$nodes$label[reactVals$graph$nodes$id == reactVals$myNode]
        children <- all_my_children(reactVals$graph, reactVals$myNode)
        # attach children 
        if (!is.null(children) & !is_empty(children)) {
          
          row_ids <- reactVals$graph$nodes$id %in% children
          selected_labels <- c(reactVals$graph$nodes$label[row_ids], selected_labels)
        } 

        countsTable <- reactVals$counts_table
        ## aggregate the clusters by name:
        countsTable['cell'] <- df_expr$cell
        # merge and aggregate by cell
        countsTable <- aggregate(. ~ cell, countsTable, sum)
        
        rownames(countsTable) <- countsTable$cell
        countsTable$cell <- NULL

        props_table <- t(t(countsTable) / colSums(countsTable)) * 100
        
        # check if any of the nodes is empty and therefore not in the props table
        # this cause the
        # make sure that sleecte labels are larger than on, otherwise it is
        # not a vector and treat it as single value
        if (length(selected_labels) > 1) {
          filter_vec <- selected_labels %in% rownames(props_table)
          props_table <- props_table[selected_labels[filter_vec],]
        }
        else {
          props_table <- t(as.data.frame(props_table[selected_labels,]))
        }
        
        ## subset the props table based on the selected node
        props_table <- as.data.frame(colSums(props_table))
        mm <- match(rownames(props_table), md$sample_id)
        props_table$cond <- as.factor(md$condition[mm])
        
        combCond <- combn(levels(props_table$cond), 2)
        
        names(props_table) <- c("value", "cond")

        condPairList <- get_pairs(levels(props_table$cond))
        
        ggplot(props_table, aes(x = cond, y = value, fill=cond))+
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.2) +
          xlab("Condition") +
          ylab("Proportion") +
          geom_signif(
            # comparisons = list(c("HC", "PD")),
            comparisons = condPairList,
            map_signif_level = F, 
            textsize = 6,
            step_increase = 0.1
          ) +
          theme_classic() +
          theme(
            plot.title = element_text(size = 22),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            legend.position = "none") +
          ggtitle(reactVals$graph$nodes$label[reactVals$myNode])
      }
    })
    



    # Load Expression Demo Data ---------------------------------------------
    observeEvent(input$btnLoadDemoData, {
      
      df_expr <<- df_expr_demoData
      cell_freq <<- cluster_freq_demoData
      
      lineage_marker <<- colnames(df_expr)
      lineage_marker_raw <<- paste0(lineage_marker, "_raw")
      
      df_expr <<- createExpressionDF(df_expr, cell_freq)
      reactVals$graph <- initTree()
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "loading Data Cluster Expression Demo Data...", value = 0.2)
      
      reactVals$graph <- initTree()
      
      set.seed(1234)
      progress$set(message = "Building the UMAP...", value = 0.3)
      umap_results <- buildUMAP(df_expr[, lineage_marker_raw]) 
      dr_umap(umap_results)

      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)

      initExprData()
    })
    
    blank_plot <- function(msg = "Workspace cleared") {
      function() {
        plot.new()
        title(msg, line = -1)
        box(col = "#ddd")
      }
    }
    
    # Clear Workspace ----
    observeEvent(input$btnClearWorkspace, {

      # Global static values
      lineage_marker <<- character(0)
      lineage_marker_raw <<- NULL
      df_expr <<- NULL
      cell_freq <<- NULL
      # dr_umap <<- NULL
      dr_umap <- reactiveVal(NULL)   # starts empty
      rFreqs <- reactiveVal(NULL)
      init_app <<- TRUE
      # 2 column data.frame containing old_cluster and new_cluster IDs
      mergeTableCatalyst <<- NULL
      
      
      reactVals$th <- NULL
      reactVals$myTH <- NULL
      reactVals$md <- NULL
      reactVals$counts_table <- NULL
      reactVals$DA_result_table <- NULL
      reactVals$DA_interactive_table <- NULL
      reactVals$graph <- NULL
      reactVals$hm <- NULL
      reactVals$annotationlist <- NULL
      reactVals$sce <- NULL
      reactVals$metaClustLevel <- NULL
      
      # actively blank visible outputs
      output$umap_tree <- renderPlot(blank_plot()())
      # output$boxplot <- renderPlot(blank_plot()())
      
      # turn off any “ready” flags you use
      reactVals$ready <- FALSE
      
      # reset file inputs -> status lights go red automatically
      shinyjs::reset("gs_upload_form")
      shinyjs::reset("remotesom_form")
      
      # disable import buttons until new files are picked
      shinyjs::reset("fMarkerExpr")
      shinyjs::reset("cluster_freq")
      shinyjs::disable("btnImportExpr")
      shinyjs::disable("btnImportRemoteSom")
      shinyjs::reset("remotesom_features")
      shinyjs::reset("remotesom_medians")
      shinyjs::reset("remotesom_counts")
      
      # (optional) clear plots/tables you render
      # output$umap <- renderPlot({ plot.new(); title("Workspace cleared") })
      showNotification("Workspace cleared.", type = "message")
      
      # showNotification("Workspace cleared.", type = "message")
    })
    
    observe({
      if (has_file(input$fMarkerExpr) && has_file(input$cluster_freq))
        shinyjs::enable("btnImportExpr") else shinyjs::disable("btnImportExpr")
    })
    
    observe({
      if (has_file(input$remotesom_features) && has_file(input$remotesom_counts) && has_file(input$remotesom_medians))
        shinyjs::enable("btnImportRemoteSom") else shinyjs::disable("btnImportRemoteSom")
    })
    
    # Load Annotated Expr Demo Data -----------------------------------------
    observeEvent(input$btnLoadAnnoData, {
      
      df_expr <<- df_expr_demoData
      cell_freq <<- cluster_freq_demoData
      
      lineage_marker <<- colnames(df_expr)
      lineage_marker_raw <<- paste0(lineage_marker, "_raw")
      
      df_expr <<- createExpressionDF(df_expr, cell_freq)
      reactVals$graph <- initTree()
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "loading Data Cluster Expression Demo Data...", value = 0.2)
      
      reactVals$graph <- initTree()
      
      set.seed(1234)
      progress$set(message = "Building the UMAP...", value = 0.3)
      umap_result <- buildUMAP(df_expr[, lineage_marker_raw]) 
      dr_umap(umap_result)

      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)

      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = lineage_marker,
        selected = NULL
      )
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = lineage_marker,
        selected = NULL
      )
      
      ## Load demo thresholds
      th <- marker_th_demo_data

      th$X <- NULL
      th$color <- "blue"
      
      th[th$bi_mod < 0.555, "color"] <- "red"
      
      rownames(th) <- th$cell
      reactVals$th<- th

      ## Load metadata
      md <<- meta_demo_data
      md$X <- NULL
      reactVals$md<- md

      ## Load counts table
      ct <<- cluster_counts_demoData
      ct$X <- NULL
      reactVals$counts_table <- ct

      ## Load annotaiton Tree
      df_nodes <- nodes_demo_data
      df_edges <- edges_demo_data

      reactVals$graph <- getGraphFromLoad(df_nodes, df_edges)

      df_expr$cell <<- rebuiltTree(reactVals$graph, df_expr, reactVals$th)
      reactVals$annotationlist <- df_nodes$label

      reactVals$hm <- df_expr[, lineage_marker]
    })
    
    # Save Workspace ----------------------------------------------------------
    
    output$btnSaveWorkspace <- downloadHandler(
      filename = function() {
        paste0("cycadas-workspace_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".rds")
      },
      content = function(file) {
        
        # browser()
        
        # helpers ----
        has_df <- function(x) {
          is.data.frame(x) && nrow(x) > 0 && ncol(x) > 0
        }
        
        if (!is.null(reactVals$graph$nodes)) {
          export_df_nodes <- reactVals$graph$nodes
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
        } else {
          export_df_nodes <- NULL
        }

        state <- list(
          schema_version = "1.0",
          app_version    = as.character(utils::packageVersion("cycadas")),
          saved_at       = Sys.time(),
          median_expr    = if (has_df(df_expr)) df_expr else NULL,
          lineage_marker = if (!is.null(lineage_marker)) lineage_marker else NULL,
          lineage_marker_raw = if (!is.null(lineage_marker_raw)) lineage_marker_raw else NULL,
          cl_freq        = if (has_df(cell_freq)) cell_freq else NULL,
          metadata       = if (has_df(reactVals$md)) reactVals$md else NULL,
          thresholds     = if (has_df(reactVals$th)) reactVals$th else NULL,
          counts_table   = if (has_df(reactVals$counts_table)) reactVals$counts_table else NULL,
          nodes_df       = export_df_nodes,
          edges_df       = if (has_df(reactVals$graph$edges)) reactVals$graph$edges else NULL,
          annotation     = if (!is.null(reactVals$annotationlist)) reactVals$annotationlist else NULL,
          umap_coords    = if (has_df(dr_umap())) dr_umap() else NULL  
        )
        save_workspace(file, state)
      }
    )
    
    # Load Workspace ----------------------------------------------------------
    observeEvent(input$btnLoadWorkspace, {
      
      # browser()
      
      req(input$btnLoadWorkspace$datapath)
      ws <- load_workspace(input$btnLoadWorkspace$datapath)
      
      # restore into your app's reactives
      if (!is.null(ws$median_expr) && !is.null(ws$cl_freq)) {
        df_expr <<- ws$median_expr
        cell_freq <<- ws$cl_freq
        
        # Assign column names to variables
        lineage_marker <<- ws$lineage_marker
        lineage_marker_raw <<- ws$lineage_marker_raw
        
        reactVals$graph <- initTree()
      } 
      
      if (!is.null(ws$umap_coords)) dr_umap(ws$umap_coords)
      if (!is.null(ws$annotation)) reactVals$annotationlist <- ws$annotation
      if (!is.null(ws$counts_table)) reactVals$counts_table <- ws$counts_table
      if (!is.null(ws$metadata)) reactVals$md <- ws$metadata
      
      if (!is.null(ws$thresholds)) {
        ws$thresholds$X <- NULL
        reactVals$th  <- ws$thresholds
      }
      
      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)
      
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = lineage_marker,
        selected = NULL
      )
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = lineage_marker,
        selected = NULL
      )
      
      # reactVals$counts_table <- ws$counts_table
      # 
      if (!is.null(ws$nodes_df) && !is.null(ws$edges_df)) {
        reactVals$graph <- getGraphFromLoad(ws$nodes_df, ws$edges_df)
        # df_expr$cell <<- rebuiltTree(reactVals$graph, df_expr, reactVals$th)
        # browser()
        # reactVals$annotationlist <- df_expr$cell
      } else {
        # reactVals$annotationlist <- c("Unassigned")  
      }

      showNotification("Workspace loaded.", type = "message")
    })
    
    
    # Workspace status --------------------------------------------------------
    # 
    has_df <- function(x) is.data.frame(x) && nrow(x) > 0 && ncol(x) > 0
    has_val <- function(x) !is.null(x) && length(x) > 0
    
    status_df <- reactive({
      data.frame(
        Item    = c("Median expression",
                    "Metadata",
                    "Thresholds",
                    "Cell Frequencies",
                    "Cluster Counts table",
                    "CATALYST object"),
        Present = c(has_df(df_expr),
                    has_df(reactVals$md),
                    has_df(reactVals$th),
                    has_val(cell_freq),
                    has_df(reactVals$counts_table),
                    has_val(reactVals$sce)),
        Details = c(
          if (has_df(df_expr))      sprintf("%d×%d", nrow(df_expr), ncol(df_expr)) else "-",
          if (has_df(reactVals$md)) sprintf("%d×%d", nrow(reactVals$md), ncol(reactVals$md)) else "-",
          if (has_df(reactVals$th)) sprintf("%d×%d", nrow(reactVals$th), ncol(reactVals$th)) else "-",
          if (has_df(cell_freq)) sprintf("%d×%d", nrow(cell_freq), ncol(cell_freq)) else "-",
          if (has_df(reactVals$counts_table)) sprintf("%d×%d", nrow(reactVals$counts_table), ncol(reactVals$counts_table)) else "-",
          if (has_val(reactVals$sce)) class(reactVals$sce)[1] else "-"
        ),
        stringsAsFactors = FALSE
      )
    })
    
    output$ws_status <- renderUI({
      df <- status_df()
      # pretty HTML list with colored badges
      items <- lapply(seq_len(nrow(df)), function(i) {
        ok <- isTRUE(df$Present[i])
        badge_cl <- if (ok) "bg-green" else "bg-red"
        icon_cl  <- if (ok) "check-circle" else "times-circle"
        tags$li(
          style = "margin: 6px 0;",
          strong(df$Item[i]), " — ",
          tags$span(class = paste("badge", badge_cl), if (ok) "Present" else "Empty"),
          " ",
          tags$i(class = paste("fa", icon_cl), style = "margin-left:6px;"),
          tags$span(style="margin-left:10px;color:#666;", df$Details[i])
        )
      })
      tags$ul(class = "list-unstyled", items)
    })
    
    
    # start disabled (optional; ensures correct state on app load)
    shinyjs::disable("btnImportRemoteSom")
    # toggle button enabled/disabled as files are chosen
    observe({
      if (remote_ready()) shinyjs::enable("btnImportRemoteSom") else shinyjs::disable("btnImportRemoteSom")
    })
    
    # live status lights
    output$remotesom_status <- renderUI({
      tagList(
        status_item(has_file(input$remotesom_features), "Feature Names"),
        status_item(has_file(input$remotesom_counts),   "Cluster Counts"),
        status_item(has_file(input$remotesom_medians),  "Median Expression"),
        tags$div(style="margin-top:8px;",
                 if (remote_ready())
                   tags$span("All required files present — ready to import.", style="color:#28a745;")
                 else
                   tags$span("Select all three files to enable Import.", style="color:#dc3545;")
        )
      )
    })
    
    
    # all three required? ----
    remote_ready <- reactive({
      has_file(input$remotesom_features) &&
        has_file(input$remotesom_counts) &&
        has_file(input$remotesom_medians)
    })
    
    # Load Remotesom ----
    observeEvent(input$btnImportRemoteSom, {
      req(remote_ready())  # double safety
      
      # read files
      features_path <- input$remotesom_features$datapath
      counts_path   <- input$remotesom_counts$datapath
      medians_path  <- input$remotesom_medians$datapath
      
      # ---
      req(input$remotesom_features$datapath)
      
      feat <- tryCatch(read_json(input$remotesom_features$datapath,
                                 simplifyVector = TRUE),
                       error = function(e) NULL)
      
      if (length(feat)) {
        lineage_marker <<- feat
        lineage_marker_raw <<- paste0(lineage_marker, "_raw")
        
      } else {
        showNotification("feature-names.json could not be parsed.", type = "error")
      }
      # ---
      rMat <- tryCatch(read_json(input$remotesom_medians$datapath,
                                 simplifyVector = TRUE),
                       error = function(e) NULL)
      
      rMat <- as.data.frame(rMat)
      colnames(rMat) <- lineage_marker
      df_expr <<- createExpressionDF(as.data.frame(rMat), cell_freq)

      reactVals$graph <- initTree()
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "loading Data Cluster Expression Demo Data...", value = 0.2)
      
      reactVals$graph <- initTree()
      
      set.seed(1234)
      progress$set(message = "Building the UMAP...", value = 0.3)
      umap_results <- buildUMAP(df_expr[, lineage_marker_raw])
      dr_umap(umap_results)
      
      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)
      
      initExprData()
      # ---
      req(input$remotesom_counts$datapath)
      
      rCounts <- tryCatch(read_json(input$remotesom_counts$datapath,
                                    simplifyVector = TRUE),
                          error = function(e) NULL)
      
      tmp <- as.data.frame(rCounts / sum(rCounts))
      tmp$cluster <- 1:nrow(tmp)
      colnames(tmp) <- c("clustering_prop", "cluster")

      cell_freq <<- tmp
      
      showNotification("RemoteSOM data imported.", type = "message")
      
      # (optional) disable button until user changes files again
      shinyjs::disable("btnImportRemoteSom")
    })
    

  }

  shinyApp(
    ui = ui,
    server = server
  )
}


