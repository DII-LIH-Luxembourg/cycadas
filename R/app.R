
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


cycadas <- function() {

  # browser()

  reactVals <- reactiveValues(th = NULL,
                              myTH = NULL, # user selected threshold
                              md = NULL,
                              counts_table = NULL,
                              DA_result_table = NULL,
                              DA_interactive_table = NULL,
                              graph = NULL,
                              hm = NULL,
                              merged_prop_table = NULL,
                              annotationlist = NULL
                              )
  
  # Global static values
  lineage_marker <<- character(0)
  df_expr <<- NULL
  dr_umap <<- NULL

  # Server function -------------------------
  server = function(input, output, session) {
    
    # Function: Update cluster labels -----------------------------------------
    updateClusterLabels <- function(mydf) {

      mysum <- sum(cell_freq[rownames(mydf),]$clustering_prop)

      output$progressBox <- renderText({
        paste0(mysum, "%")
      })
      output$progressBox2 <- renderText({
        paste0(dim(mydf)[1], "Cluster")
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
      
      # browser()

      output$mynetworkid <- renderVisNetwork({
        visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
          visEdges(arrows = "from") %>%
          visHierarchicalLayout()
      })
    }

    # Upload the marker expression file ---------------------------------------
    observeEvent(c(input$fMarkerExpr, input$cluster_freq), {

      reactVals$graph <- initTree()
      # browser()

      req(input$fMarkerExpr)
      req(input$cluster_freq)

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

      ## load only at start to fill the picker list
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

    })
    
    # Observe MenutItems ------------------------------------------------------
    observeEvent(input$tabs, {
      
      # browser()    

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

          updateClusterLabels(tmp)

          plotTree()

        }
      }
    })

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

      node <- reactVals$graph$nodes %>% filter(label == input$parentPicker)
      myid <- node$id

      filterPosMarkers <- unlist(node$pm)
      filterNegMarkers <- unlist(node$nm)
      
      parent <- input$parentPicker
      
      # receive the parent settings, resp. parent hm
      # filter hm by parent cell name
      tmp <- df_expr[df_expr$cell == parent, lineage_marker]
      tmp <- filterHM(tmp, input$treePickerPos, input$treePickerNeg, reactVals$th)

      updateClusterLabels(tmp)

      reactVals$hm <- tmp

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

      output$umap_tree <-
        renderPlot(
          ggplot(dr_umap, aes(
            x = u1, y = u2, color = ClusterSelection
          )) +
            geom_point(size = 1.0) +
            theme_bw() +
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

      node <- reactVals$graph$nodes %>% filter(label == input$updateNodePicker)
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
    observeEvent(input$updateNodeBtn, {

      node <- reactVals$graph$nodes %>% filter(label == input$parentPicker)
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


    importNodes <- function(id, label, pm, nm, parent_id, color) {

      if(id == parent_id) {return()}

      parent_row <- reactVals$graph$nodes %>% filter(id == parent_id)
      parent_name <- parent_row$label

      add_node(reactVals$graph, parent_name, label, color)
    }

    # Load Annotation Tree ----
    observeEvent(input$btnImportTree, {

      # browser()

      req(input$fNodes)
      req(input$fEdges)
      # req(input$fAnno)

      req(input$fMarkerExpr)
      req(input$cluster_freq)

      df_nodes <- read.csv(input$fNodes$datapath)
      df_edges <- read.csv(input$fEdges$datapath)
      # df_anno <- read.csv(input$fAnno$datapath)

      # df_nodes$pm[is.na(df_nodes$pm)] <- ""
      # df_nodes$pm <- strsplit(df_nodes$pm, "\\|")
      # 
      # df_nodes$nm[is.na(df_nodes$nm)] <- ""
      # df_nodes$nm <- strsplit(df_nodes$nm, "\\|")
      # 
      # reactVals$graph$nodes <- df_nodes
      # reactVals$graph$edges <- df_edges
      
      reactVals$graph <- getGraphFromLoad(df_nodes, df_edges)
      
      reactVals <- rebuiltTree(reactVals)

      # df01Tree <<- df_anno

      annotationlist <<- as.list(df_nodes$label)

    })


    observeEvent(input$exportTreeGraphics, {

      # browser()

      # visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
      #   visEdges(arrows = "from") %>%
      #   visHierarchicalLayout() %>%
      #   visExport(type = "png", name = "export-network",
      #             float = "left", label = "Save network", background = "purple", style= "")


    })

    # Export Annotation Tree ----
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
        
        # browser()
        
        export_freq <- data.frame(cluster=1:length(annotaionDF$cell),
                                  clustering_prop = annotaionDF$clusterSize)

        download_list <- list(annTable = reactVals,
                              nodesTable = export_df_nodes,
                              edgesTable = reactVals$graph$edges,
                              freq = export_freq)
        # download_list <- list(annTable = df01Tree["cell"],
        #                       nodesTable = export_df_nodes,
        #                       edgesTable = reactVals$graph$edges)

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

    ##
    # Server - Thresholds Tab -------------------------------------------------
    ##
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

    ## react to selected rows in thresholds list: -----------------------------
    ## from marker, plot the expression in scatterplot or histogram
    observeEvent(input$table_rows_selected, {

      # browser()

      selRow <- reactVals$th[input$table_rows_selected,]
      myTH <- reactVals$th[input$table_rows_selected, "threshold"]
      myColor <- reactVals$th[input$table_rows_selected, "color"]

      marker_expr <- as.data.frame(df_expr[, selRow$cell])
      marker_expr <- cbind(marker_expr, rnorm(1:nrow(marker_expr)))
      
      PlottingThreshold(marker_expr, myTH, myColor)

      # myRenderFunction(marker_expr, selRow$threshold, selRow$color)
    })
    
    ##
    # Observe click scatterplot -----------------------------------------------
    ##
    observeEvent(input$plot_click, {

      # browser()

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

    ##
    # PlottingThreshold  function ---------------------------------------------
    ##
    PlottingThreshold <- function(me, myTH, myCol){

      output$plot <- renderPlot({
        set.seed(1)

        ggplot(me, aes_string(x=me[,1], y=me[,2])) +
          geom_point(size=1) +
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
      #
      output$plot2 <- renderPlot({

        ggplot(me, aes_string(x = me[, 1])) +
          geom_histogram(bins = 80) +
          labs(x = "Scale 0 to 1") +
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


    # Load the annotation file ------------------------------------------------
    # observeEvent(input$annTable,{
    #   annData <<- read.csv(input$annTable$datapath)
    #   annData$X <- NULL
    #   at$data <- annData
    # })
    
    # Save the DA Table -------------------------------------------------------
    output$exportDA <- downloadHandler(
      filename = function() {
        paste("DA_Table_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(reactVals$DA_result_table, file)
      }
    )
    
    # get_merged_prop_table ---------------------------------------------------
    get_merged_prop_table <- function() {
      
      req(reactVals$counts_table)
      
      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df_expr$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)
      
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

      # browser()

      th <<- read.csv(input$fTH$datapath)
      th$X <- NULL
      th$color <- "blue"

      th[th$bi_mod < 0.555, "color"] <- "red"

      rownames(th) <- th$cell

      reactVals$th<- th
      
      # re-calculate the tree after threshold upload
      req(input$fMarkerExpr)
      
      reactVals <- rebuiltTree(reactVals)
      

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

    ## create Phenotype name from picker e.g. CD4+CD8+CD19 --------------------
    # observeEvent(input$myPickerPos, {
    #   output$r1 <- renderText(setPhenotypeName(input$myPickerPos, "pos", ph_name))
    # })
    # observeEvent(input$myPickerNeg, {
    #   output$r2 <- renderText(setPhenotypeName(input$myPickerNeg, "neg", ph_name))
    # })

    ## Annotation Table -------------------------------------------------------
    # observeEvent(c(input$myPickerPos, input$myPickerNeg, reactVals$th),{
    #   req(input$fMarkerExpr)
    #   req(input$cluster_freq)
    # 
    #   tmp=filterHM(df01,input$myPickerPos, input$myPickerNeg, reactVals$th)
    # 
    #   # update myDF to the filtered HM
    #   myDF <<- tmp
    # 
    #   mydf=annotaionDF[rownames(tmp),]
    #   rownames(mydf)=rownames(tmp)
    # 
    #   # update heatmap plot ----
    #   output$hm <- renderPlot({
    #     if (dim(tmp)[1] < 2) {
    #       pheatmap(tmp, cluster_cols = F, cluster_rows = F)
    #     } else {
    #       message(dim(tmp)[1])
    #       pheatmap(tmp, cluster_cols = F)
    #     }
    #   })
    # 
    #   # update umap plot ----
    #   myColor <- filterColor(df_global,tmp)
    #   Temp1<<-myColor
    #   Temp2<<-dr_umap
    # 
    #   output$umap <-
    #     renderPlot(
    #       ggplot(dr_umap, aes(
    #         x = u1, y = u2, color = myColor
    #       )) +
    #         geom_point(size = 1.0) +
    #         theme_bw() +
    #         theme(legend.text = element_text(size =
    #                                            8)) +
    #         guides(color = guide_legend(override.aes = list(size = 4)))
    #     )
    #   # update table ----
    #   output$tableAnnotation <- DT::renderDT(
    #     at$data,
    #     extensions = c('Buttons', 'Scroller'),
    #     options = list(
    #       dom = 'Bfrtip',
    #       buttons = c('csv'),
    #       paging = FALSE
    #     )
    #   )
    # 
    #   # update valueboxes ----
    #   mysum <- sum(mydf$clusterSize)
    #   output$progressBox <- renderValueBox({
    #     valueBox(paste0(mysum, "%"), "Selected", icon = icon("list"),color = "purple")
    #   })
    #   output$progressBox2 <- renderValueBox({
    #     valueBox(dim(mydf)[1], "Cluster", icon = icon("list"),color = "purple")
    #   })
    # })

    # Set Phenotype names -----------------------------------------------------
    # observeEvent(input$btnSetType, {
    # 
    #   # if cluster has already been assigned, add the new name after
    #   for (n in rownames(myDF)) {
    #     if (at$data[n, "cell"] == "unassigned") {
    #       at$data[n, "cell"] <- input$phenotype
    #     } else(
    #       at$data[n, "cell"] <- paste0(at$data[n, "cell"], "_", input$phenotype)
    #     )
    # 
    #   }
    # })

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
          color = input$markerSelect
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

    ## Do the differential abundance ------------------------------------------
    observeEvent(input$doDA, {
      # browser()
      req(reactVals$counts_table, reactVals$md)

      countsTable <- reactVals$counts_table
      # aggregate the clusters by name:
      countsTable['cell'] <- df_expr$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL

      props_table <- t(t(countsTable) / colSums(countsTable)) * 100

      # mm <- match(md$sample_id, colnames(props_table))
      mm <- match(colnames(props_table), md$sample_id)

      tmp_cond <- md$condition[mm]

      DA_df <- data.frame()
      for (i in 1:nrow(props_table)) {
        foo <- pairwise.wilcox.test(as.numeric(props_table[i,]), tmp_cond, p.adjust.method=input$correction_method)

        df <- subset(melt(foo$p.value), value!=0)
        df$cell <- rownames(countsTable)[i]
        DA_df <- rbind(DA_df, as.data.frame(df))

      }

      # browser()
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

    # Render interactive Tree
    output$interactiveTree <- renderVisNetwork({
      visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    })

    myNode <- reactiveValues(selected = '')

    observeEvent(input$current_node_id, {
      myNode$selected <<- input$current_node_id
    })

    output$selectedNode <- renderText({myNode$selected})

    output$DA_interactive_table <- renderTable({

      # browser()
      children <- all_my_children(reactVals$graph, myNode$selected)

      if (!is.null(children)) {
        selected_labels <- c(reactVals$graph$nodes$label[children], reactVals$graph$nodes$label[myNode$selected])
      } else {
        selected_labels <- reactVals$graph$nodes$label[myNode$selected]
      }

      # browser()
      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df_expr$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL

      # countsTable <- countsTable[selected_labels,]

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

      ## subest the prps tabel based on the selected node
      # props_table <- props_table[selected_labels, ]
      props_table <- t(as.data.frame(colSums(props_table)))

      # mm <- match(md$sample_id, colnames(props_table))
      mm <- match(colnames(props_table), md$sample_id)

      tmp_cond <- md$condition[mm]

      DA_df <- data.frame()

      foo <- pairwise.wilcox.test(props_table, tmp_cond, p.adjust.method="none")

      df <- subset(melt(foo$p.value), value!=0)
      df$cell <- reactVals$graph$nodes$label[myNode$selected]
      DA_df <- rbind(DA_df, as.data.frame(df))

      colnames(DA_df) <- c("Var1", "Var2", "p-value", "Cell")
      reactVals$DA_interactive_table <- DA_df

    })

    # Boxplot of interactive DA selection -------------------------------------
    output$boxplot <- renderPlot({

      # browser()
      children <- all_my_children(reactVals$graph, myNode$selected)

      if (!is.null(children)) {
        selected_labels <- c(reactVals$graph$nodes$label[children], reactVals$graph$nodes$label[myNode$selected])
      } else {
        selected_labels <- reactVals$graph$nodes$label[myNode$selected]
      }

      # browser()
      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df_expr$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL

      # countsTable <- countsTable[selected_labels,]

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

      ## subest the prps tabel based on the selected node
      # props_table <- props_table[selected_labels, ]

      props_table <- as.data.frame(colSums(props_table))


      # mm <- match(md$sample_id, rownames(props_table))
      mm <- match(rownames(props_table), md$sample_id)

      props_table$cond <- md$condition[mm]

      names(props_table) <- c("value", "cond")

      ggplot(props_table, aes(x = cond, y = value, fill=cond))+
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2) +
        xlab("Condition") +
        ylab("Proportion") +
        theme(plot.title = element_text(size = 22),
              axis.text = element_text(size = 12),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 20),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20)) +
        ggtitle(reactVals$graph$nodes$label[myNode$selected])


    })
    
    loadDemoData <- function(path_expr, path_freq) {
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "loading Data...", value = 0)
      
      ## Load median expression and cell frequencies
      ## add the 0-1 scaled columns to the DF
      df_expr <<- read.csv("data/demo_data/median_expr_1600.csv")
      cell_freq <<- read.csv("data/demo_data/cluster_freq_1600.csv")
      
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

    # Upload expression Demo Data ---------------------------------------------
    observeEvent(input$btnLoadDemoData, {

      pathExpr <- "data/demo_data/median_expr_1600.csv"
      pathFreq <- "data/demo_data/cluster_freq_1600.csv"
      
      loadDemoData(pathExpr, pathFreq)
      
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
      
      posPickerList <<- lineage_marker
      
      updateSelectInput(session, "markerSelect", "Select:", lineage_marker)

      reactVals$th <- kmeansTH(df_expr[, lineage_marker])
      reactVals$annotationlist <- c("Unassigned")

    })
    
    # Upload Annotated Expr Demo Data -----------------------------------------
    observeEvent(input$btnLoadAnnoData, {
      
      pathExpr <- "data/demo_data/median_expr_1600.csv"
      pathFreq <- "data/demo_data/cluster_freq_1600.csv"
      
      loadDemoData(pathExpr, pathFreq)
      posPickerList <<- lineage_marker

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
      th <- read.csv("data/demo_data/MarkerThresholds.csv")

      th$X <- NULL
      th$color <- "blue"
      
      th[th$bi_mod < 0.555, "color"] <- "red"
      
      rownames(th) <- th$cell
      reactVals$th<- th

      ## Load metadata
      md <<- read.csv("data/demo_data/metadata.csv")
      md$X <- NULL
      reactVals$md<- md

      ## Load counts table
      ct <<- read.csv("data/demo_data/cluster_counts_1600.csv")
      ct$X <- NULL
      reactVals$counts_table <- ct

      ## Load annotaiton Tree
      df_nodes <- read.csv("data/demo_data/nodesTable_data.csv")
      df_edges <- read.csv("data/demo_data/edgesTable_data.csv")

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

