# usethis::use_package("shiny")
# usethis::use_package("DT")
# usethis::use_package("ggplot2")
# usethis::use_package("matrixStats")
# usethis::use_package("dplyr")
# usethis::use_package("stats")
# usethis::use_package("pheatmap")
# usethis::use_package("Ckmeans.1d.dp")
# usethis::use_package("umap")
# usethis::use_package("RColorBrewer")
# usethis::use_package("shinydashboard")
# usethis::use_package("shinyWidgets")
# usethis::use_package("visNetwork")
# usethis::use_package("glue")
# usethis::use_package("purrr")
# usethis::use_package("reshape2")
# usethis::use_package("mousetrap")
# usethis::use_package("knitr")


# List of packages you want to check and install if needed
# packages_to_install <- c("shiny", "DT", "ggplot2", "matrixStats", "tidyverse", "stats", "knitr", "forcats",
#                          "pheatmap", "Ckmeans.1d.dp", "umap", "RColorBrewer", "shinydashboard", "mixtools",
#                          "shinyWidgets", "visNetwork", "glue", "purrr", "reshape2", "mousetrap", "ggpubr", 
#                          "SingleCellExperiment", "CATALYST", "shinyjs", "cluster", "mclust")
# 
# 
# # Check if each package is already installed, and install if not
# for (package in packages_to_install) {
#   if (!require(package, character.only = TRUE)) {
#     install.packages(package)
#     library(package, character.only = TRUE)
#   }
# }

shiny::addResourcePath("images", "./www")

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
  dr_umap <<- NULL
  init_app <<- TRUE
  # 2 column data.frame containing old_cluster and new_cluster IDs
  mergeTableCatalyst <<- NULL

  # Server function -------------------------
  server = function(input, output, session) {
    
    shinyjs::disable("metadiv")
    
    # browser()
    
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
    

    # Load Expression Data and UMAP -------------------------------------------------
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
      dr_umap <<- buildUMAP(df_expr[, lineage_marker_raw]) 

    }
    
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
      
      # posPickerList <<- lineage_marker
      
      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)
      
      reactVals$th <- kmeansTH(df_expr[, lineage_marker])
      reactVals$annotationlist <- c("Unassigned")
    }
    
    # Function: Update cluster labels -----------------------------------------
    updateClusterLabels <- function(mydf, my_node_id=0) {
      
      # browser()
      
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
    

    # Load the marker expression file ---------------------------------------
    observeEvent(c(input$fMarkerExpr, input$cluster_freq), {

      req(input$fMarkerExpr)
      req(input$cluster_freq)

      pathExpr <- input$fMarkerExpr$datapath
      pathFreq <- input$cluster_freq$datapath
      
      loadExprData(pathExpr, pathFreq)
      
      initExprData()
    })
    
    # Load CATALYST data ----
    observeEvent(input$sce, {
      # browser()
      req(input$sce)
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "loading Data ...", value = 0.2)
      
      reactVals$sce <- readRDS(input$sce$datapath)
      
      reactVals$metaClustLevel <- colnames(reactVals$sce@metadata$cluster_codes)[1]
      # reactVals$metaClustLevel <- mc_levels[5:length(mc_levels)]
      
      # my_sce <- readRDS(input$sce$datapath)
      # cluster_ids(sce)
      
      # print(dim(assay(reactVals$sce)))
      # print(colnames(colData(reactVals$sce)))
      # 
      # print(head(reactVals$sce@metadata$SOM_codes))
      
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
      dr_umap <<- buildUMAP(df_expr[, lineage_marker_raw]) 
      
      # df_expr
      # 
      # pathExpr <- "data/demo_data/median_expr_1600.csv"
      # pathFreq <- "data/demo_data/cluster_freq_1600.csv"
      # 
      # loadExprData(pathExpr, pathFreq)
      # posPickerList <<- lineage_marker
      
      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)
      
      # pathExpr <- "data/demo_data/median_expr_1600.csv"
      # pathFreq <- "data/demo_data/cluster_freq_1600.csv"
      # 
      # loadExprData(pathExpr, pathFreq)
      # 
      initExprData()
      
      mc_levels <- colnames(reactVals$sce@metadata$cluster_codes)
      #remove 2,3, and 4
      mc_levels <- mc_levels[-(2:5)]
      
      updatePickerInput(
        session,
        inputId = "metaLevel",
        # choices = colnames(reactVals$sce@metadata$cluster_codes),
        choices = mc_levels,
        selected = reactVals$metaClustLevel
        # selected = colnames(reactVals$sce@metadata$cluster_codes)[1]
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
        
        # browser()
        
        reactVals$metaClustLevel=input$metaLevel
        
        # select the metalvel and merge the clusters 
        somExpr <- as.data.frame(reactVals$sce@metadata$SOM_codes)
        
        somExpr$metaLevels <- reactVals$sce@metadata$cluster_codes[reactVals$metaClustLevel]
        
        somExpr <- somExpr %>%
          group_by(metaLevels) %>%
          summarise_if(is.numeric, mean)
        
        somExpr <- as.data.frame(somExpr[, lineage_marker])
        # cell_freq_tmp <- as.data.frame(cluster_ids(c_id = cluster_ids(reactVals$sce, "meta20")))
        cell_freq_tmp <- as.data.frame(table(c_id = cluster_ids(reactVals$sce, reactVals$metaClustLevel)))
        
        cell_freq_tmp$clustering_prop <- cell_freq_tmp$Freq / sum(cell_freq_tmp$Freq)
        cell_freq_tmp$Freq <- NULL
        colnames(cell_freq_tmp) <- c("cluster", "clustering_prop")
        cell_freq <<- cell_freq_tmp
        
        
        df_expr <<- createExpressionDF(somExpr, cell_freq_tmp)
        reactVals$graph <- initTree()
        
        set.seed(1234)
        # progress$set(message = "Building the UMAP...", value = 0.3)
        
        if (nrow(df_expr) > 20) {
          dr_umap <<- buildUMAP(df_expr[, lineage_marker_raw])   
        } else {
          dr_umap <<- NULL
        }
        
        
        # df_expr
        # 
        # pathExpr <- "data/demo_data/median_expr_1600.csv"
        # pathFreq <- "data/demo_data/cluster_freq_1600.csv"
        # 
        # loadExprData(pathExpr, pathFreq)
        # posPickerList <<- lineage_marker
        
        updateSelectInput(session, "markerSelect", "Select:", lineage_marker)
        
        # pathExpr <- "data/demo_data/median_expr_1600.csv"
        # pathFreq <- "data/demo_data/cluster_freq_1600.csv"
        # 
        # loadExprData(pathExpr, pathFreq)
        # 
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
      
      # x <- 1 # for debug
      
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
      
      # browser()

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
      
      # browser()
      
      nodeName <- input$parentPicker
      
      # if (selectMethod == "picker") {
      #   
      node <- reactVals$graph$nodes %>% dplyr::filter(label == nodeName)
      parent <- input$parentPicker
      
      nodeID <- reactVals$graph$nodes$id[reactVals$graph$nodes$label == nodeName]
      # } else {
      #   
      #   node <- reactVals$graph$nodes %>% dplyr::filter(id == nodeID)
      #   parent <- node$label
      # }
      # 
      filterPosMarkers <- unlist(node$pm)
      filterNegMarkers <- unlist(node$nm)
      
      # receive the parent settings, resp. parent hm
      # filter hm by parent cell name
      tmp <- df_expr[df_expr$cell == parent, lineage_marker]
      tmp <- filterHM(tmp, input$treePickerPos, input$treePickerNeg, reactVals$th)
      
      # browser()
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
          cbind(dr_umap,ClusterSelection) %>%
            dplyr::mutate(ClusterSelection=as.factor(ClusterSelection)) %>%
            dplyr::mutate(ClusterSelection=fct_relevel(ClusterSelection,"other clusters","selected phenotype")) %>%
            arrange(desc(ClusterSelection)) %>%
            ggplot( aes(x = u1, y = u2, color = ClusterSelection)) +
            geom_point(size = 1.0) +
            # ggpubr::theme_pubr() +
            theme(legend.text = element_text(size = 12),
                  legend.title = element_text(size = 20),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 20)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
          
          #     ggplot(dr_umap, aes(
          #       x = u1, y = u2, color = ClusterSelection
          #     )) +
          #       geom_point(size = 1.0) +
          #       theme_bw() +
          #       theme(legend.text = element_text(size = 12),
          #             legend.title = element_text(size = 20),
          #             axis.text = element_text(size = 12),
          #             axis.title = element_text(size = 20)) +
          #       guides(color = guide_legend(override.aes = list(size = 4)))
        )
    })
    
    # Observe Interactive Parent Node selection -------------------------------
    observeEvent(input$parent_node_id, {
      
      # browser()
      
      # parentNodeSelection("interactive", nodeID=input$parent_node_id)
      node <- reactVals$graph$nodes %>% dplyr::filter(id == input$parent_node_id)
      parent <- node$label
      
      nodeID <- reactVals$graph$nodes$id[reactVals$graph$nodes$label == parent]
      
      filterPosMarkers <- unlist(node$pm)
      filterNegMarkers <- unlist(node$nm)
      
      # receive the parent settings, resp. parent hm
      # filter hm by parent cell name
      tmp <- df_expr[df_expr$cell == parent, lineage_marker]
      tmp <- filterHM(tmp, input$treePickerPos, input$treePickerNeg, reactVals$th)
      
      # browser()
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
          cbind(dr_umap,ClusterSelection) %>%
            dplyr::mutate(ClusterSelection=as.factor(ClusterSelection)) %>%
            dplyr::mutate(ClusterSelection=fct_relevel(ClusterSelection,"other clusters","selected phenotype")) %>%
            arrange(desc(ClusterSelection)) %>%
            ggplot( aes(x = u1, y = u2, color = ClusterSelection)) +
            geom_point(size = 1.0) +
            # ggpubr::theme_pubr() +
            theme(legend.text = element_text(size = 12),
                  legend.title = element_text(size = 20),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 20)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
          
          #     ggplot(dr_umap, aes(
          #       x = u1, y = u2, color = ClusterSelection
          #     )) +
          #       geom_point(size = 1.0) +
          #       theme_bw() +
          #       theme(legend.text = element_text(size = 12),
          #             legend.title = element_text(size = 20),
          #             axis.text = element_text(size = 12),
          #             axis.title = element_text(size = 20)) +
          #       guides(color = guide_legend(override.aes = list(size = 4)))
        )
      
      
      # parentNodeSelection("picker", nodeName=input$parentPicker)
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

      # if "newNode" text-field is not empty, update the node name
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

      # visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
      #   visEdges(arrows = "from") %>%
      #   visHierarchicalLayout() %>%
      #   visExport(type = "png", name = "export-network",
      #             float = "left", label = "Save network", background = "purple", style= "")


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

        # export_freq <- data.frame(cluster=1:length(annotaionDF$cell),
        #                           clustering_prop = annotaionDF$clusterSize)

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
      
      # browser()
      
      if (!is.null(reactVals$th)) {
        # TODO: change naming: not only k-means
        reactVals$th <- updateTH(df_expr, reactVals$th, th_mode=input$th_method)
        # df_expr$cell <<- rebuiltTree(reactVals$graph, df_expr, reactVals$th)      
          
          # req(input$table_rows_selected)
          # selRow <- reactVals$th[input$table_rows_selected,]
          # 
          # reactVals$th[input$table_rows_selected, "threshold"] <- round(input$plot_click$x, 3)
        # req(input$table_rows_selected)  
        # myTH <- reactVals$th[input$table_rows_selected, "threshold"]
        # myColor <- reactVals$th[input$table_rows_selected, "color"]
        # marker_expr <- as.data.frame(df_expr[, selRow$cell])
        # marker_expr <- cbind(marker_expr, rnorm(1:nrow(marker_expr)))
        # 
        # PlottingThreshold(marker_expr, myTH, myColor)
        
        # df_expr$cell <<- rebuiltTree(reactVals$graph, df_expr, reactVals$th)        
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
      
      
      # rename the cells for the remaining definition
      # countsTable$cell <- reactVals$DA_result_table$Naming
      
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
    output$umap2 <- renderPlot(
      ggplot(dr_umap, aes(x = u1, y = u2)) +
        geom_point(size = 1.0) +
        theme_bw() +
        theme(legend.text=element_text(size=8),
              legend.title = element_text(size = 20),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 20)) +
        guides(color = guide_legend(override.aes = list(size = 4)))
    )

    # Heatmap from UMAP selection ---------------------------------------------
    output$hm2 <- renderPlot({
      
      my_new_hm <- round(brushedPoints(dr_umap, input$umap2_brush), 2)

      if (dim(my_new_hm)[1] > 0) {

        pheatmap(my_new_hm %>% select(-c("cluster_number", "u1", "u2")), cluster_cols = F)
      } else {
        ggplot() + theme_void() + ggtitle("Select area on Umap to plot Heatmap")
      }
    })

    # Data Table from UMAP selection ------------------------------------------
    output$umap_data <-
      DT::renderDT(server = FALSE, {
        DT::datatable(
          round(brushedPoints(dr_umap, input$umap2_brush), 2),
          filter = 'top',
          extensions = 'Buttons',
          options = list(
            scrollY = 600,
            scrollX = TRUE,
            dom = '<"float-left"l><"float-right"f>rt<"row"<"col-sm-4"B><"col-sm-4"i><"col-sm-4"p>>',
            lengthMenu =  list(c(10, 25, 50,-1),
                               c('10', '25', '50', 'All')),
            scrollCollapse = TRUE,
            lengthChange = TRUE,
            widthChange = TRUE,
            rownames = TRUE
          )
        )
      })

    # UMAP Marker Expression Tab ----------------------------------------------
    output$umap3 <-
      renderPlot(
        
        ggplot(dr_umap, aes_string(
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
      )

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

      countsTable <- reactVals$counts_table
      # aggregate the clusters by name:
      countsTable['cell'] <- df_expr$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL

      props_table <- t(t(countsTable) / colSums(countsTable)) * 100

      mm <- match(colnames(props_table), md$sample_id)
      tmp_cond <- md$condition[mm]

      DA_df <- data.frame()
      
      for (i in 1:nrow(props_table)) {
        foo <- pairwise.wilcox.test(as.numeric(props_table[i,]), tmp_cond, p.adjust.method=input$correction_method)

        df <- subset(melt(foo$p.value), value!=0)
        df$cell <- rownames(countsTable)[i]
        DA_df <- rbind(DA_df, as.data.frame(df))
      }

      # change the name of nodes:
      # if node has children add "_remaining"
      my_list <- lapply(DA_df$cell, function(x) {
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

      DA_df$list <- unlist(my_list)
      colnames(DA_df) <- c("Cond1", "Cond2", "p-value", "Cell", "Naming")
      reactVals$DA_result_table <- DA_df
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
      reactVals$myNode <- input$current_node_id
      
      if (reactVals$myNode != 1) {
        doInteractiveDA()  
      }
      
    })

    output$selectedNode <- renderText({reactVals$myNode})
    
    doInteractiveDA <- function() {  

      children <- all_my_children(reactVals$graph, reactVals$myNode)

      if (!is.null(children)) {
        selected_labels <- c(reactVals$graph$nodes$label[children], reactVals$graph$nodes$label[reactVals$myNode])
      } else {
        selected_labels <- reactVals$graph$nodes$label[reactVals$myNode]
        # TODO: find a better method 
        # selected_labels <- "Unassigned"
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

        children <- all_my_children(reactVals$graph, reactVals$myNode)
        
        if (!is.null(children)) {
          selected_labels <- c(reactVals$graph$nodes$label[children], reactVals$graph$nodes$label[reactVals$myNode])
        } else {
          selected_labels <- reactVals$graph$nodes$label[reactVals$myNode]
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
        
        # browser()
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
          # geom_pwc(aes(group = cond),method = "wilcox_test", label = "Wilcoxon, italic(p)= {p}")+
          # ggpubr::theme_pubr() +
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
    
    # Load Expression Data and UMAP -------------------------------------------------
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
      dr_umap <<- buildUMAP(df_expr[, lineage_marker_raw]) 
      
    }

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
      dr_umap <<- buildUMAP(df_expr[, lineage_marker_raw]) 
      
      # df_expr
      # 
      # pathExpr <- "data/demo_data/median_expr_1600.csv"
      # pathFreq <- "data/demo_data/cluster_freq_1600.csv"
      # 
      # loadExprData(pathExpr, pathFreq)
      # posPickerList <<- lineage_marker
      
      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)

      # pathExpr <- "data/demo_data/median_expr_1600.csv"
      # pathFreq <- "data/demo_data/cluster_freq_1600.csv"
      # 
      # loadExprData(pathExpr, pathFreq)
      # 
      initExprData()
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
      dr_umap <<- buildUMAP(df_expr[, lineage_marker_raw]) 
      
      # df_expr
      # 
      # pathExpr <- "data/demo_data/median_expr_1600.csv"
      # pathFreq <- "data/demo_data/cluster_freq_1600.csv"
      # 
      # loadExprData(pathExpr, pathFreq)
      # posPickerList <<- lineage_marker

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
    
  }

  shinyApp(
    ui = ui,
    server = server
  )
}


