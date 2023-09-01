library(shiny)
library(DT)
library(ggplot2)
library(matrixStats)
library(tidyverse)
library(stats)
library(pheatmap)
library(Ckmeans.1d.dp)
library(umap)
library(RColorBrewer)
library(shinydashboard)
library(shinyWidgets)
library(visNetwork)
library(glue)
library(purrr)
library(reshape2)

cycadas <- function() {

  reactVals <- reactiveValues(th = NULL,
                              myTH = NULL,
                              md = NULL,
                              counts_table = NULL,
                              DA_result_table = NULL,
                              DA_interactive_table = NULL,
                              graph = NULL)
  # constructs a string of positive or negative markers
  ph_name <<- ""
  # initDATreePage <<- T

  # Server function ----
  server = function(input, output, session) {

    # Function: Update cluster labels ----
    updateClusterLabels <- function(mydf) {
      # browser()
      mysum <- sum(cell_freq[rownames(mydf),]$clustering_prop)

      output$progressBox <- renderText({
        paste0(mysum, "%")
      })
      output$progressBox2 <- renderText({
        paste0(dim(mydf)[1], "Cluster")
      })
    }

    # Function: Plot the annotation tree ----
    plotTree <- function() {

      # before plotting the tree, update its properties
      # based on the annotated DF
      for(i in 1:nrow(reactVals$graph$nodes)) {
        l <- reactVals$graph$nodes$label[i]

        # get the number of clusters with that label
        nLabel <- sum(df01Tree$cell == l)
        mysum <- sum(cell_freq[rownames(df01Tree),]$frequency)

        if (nLabel == 0) reactVals$graph$nodes[i,]$color <- "red"

      }

      if(input$treeSwitch == T) {
        output$mynetworkid <- renderVisNetwork({
          visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
            visEdges(arrows = "from")
        })

      } else {
        output$mynetworkid <- renderVisNetwork({
          visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
            visEdges(arrows = "from") %>%
            visHierarchicalLayout()
        })
      }

    }

    # Upload the marker expression file ----
    observeEvent(c(input$fMarkerExpr, input$cluster_freq), {

      # browser()

      req(input$fMarkerExpr)
      req(input$cluster_freq)

      df <- read.csv(input$fMarkerExpr$datapath)
      df_global <<- df

      # create initial master node of all Unassigned clusters
      nodes <- tibble(id = 1,
                      label = "Unassigned",
                      pm = list(""),
                      nm = list(""),
                      color = "blue"
      )

      edges <- data.frame(from = numeric(), to = numeric())
      reactVals$graph <- list(nodes = nodes, edges = edges)

      annotationlist <<- list("Unassigned")

      cell_freq <<- read.csv(input$cluster_freq$datapath) %>%
        mutate(total = sum(column_name)) %>%
        mutate(frequency = round(column_name / total * 100, 2))

      labels_row <-
        paste0(rownames(df), " (", cell_freq$frequency , "%)")

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
                                 clusterSize = cell_freq$frequency)
      at <<- reactiveValues(data = annotaionDF, dr_umap = dr_umap)

    })
    # Observe MenutItems ----
    observeEvent(input$tabs, {
      # req(input$fMarkerExpr)
      # req(input$cluster_freq)
      # req(df01)
      # req(graph)

      if (exists("df01")) {
        if(input$tabs=="thresholds") { ## Thresholds ----

          if (is.null(reactVals$th) & !is.null(df01)) {
            reactVals$th <- kmeansTH(df01)
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
        if(name %in% reactVals$graph$nodes$label) {
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

          reactVals$graph <- add_node(reactVals$graph,parent,name,posmarker,negmarker,"blue")

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

    # Observe Parent Node selection ----
    observeEvent(input$parentPicker, {

      # browser()

      node <- reactVals$graph$nodes %>% filter(label == input$parentPicker)
      myid <- node$id

      filterPosMarkers <- unlist(node$pm)
      filterNegMarkers <- unlist(node$nm)

      # collect all positive and negative markers
      # upwards from parent
      # go through edges until we reach one level below master node
      while (myid > 1) {
        print(myid)

        edge <- reactVals$graph$edges %>% filter(from == myid)
        next_id <- edge$to

        myid <- next_id
        next_node <- reactVals$graph$nodes %>% filter(id == next_id)
        filterPosMarkers <- append(filterPosMarkers, unlist(next_node$pm))
        filterNegMarkers <- c(filterNegMarkers, unlist(next_node$nm))

      }
      # print(filterPosMarkers)
      # print(filterNegMarkers)

      # remove the empty strings in the markers
      filterPosMarkers <- filterPosMarkers[nzchar(filterPosMarkers)]
      filterNegMarkers <- filterNegMarkers[nzchar(filterNegMarkers)]

      tmp <- filterHM(df01Tree,filterPosMarkers, filterNegMarkers, reactVals$th)
      tmp <- tmp[tmp$cell == node$label ,]

      updateClusterLabels(tmp)

      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerPos",
        choices = colnames(df01),
        selected = NULL,
        disabledChoices = filterPosMarkers
      )
      updateCheckboxGroupButtons(
        session,
        inputId = "treePickerNeg",
        choices = colnames(df01),
        selected = NULL,
        disabledChoices = filterNegMarkers
      )

      # update picker heatmap plot ----
      output$hm_tree <- renderPlot({
        if (dim(tmp)[1] < 2) {
          pheatmap(tmp[allMarkers], cluster_cols = F, cluster_rows = F)
        } else {
          message(dim(tmp)[1])
          pheatmap(tmp[allMarkers], cluster_cols = F)
        }
      })

      # update umap plot for Tree ----
      ClusterSelection <- filterColor(df_global,tmp)

      output$umap_tree <-
        renderPlot(
          ggplot(dr_umap, aes(
            x = u1, y = u2, color = ClusterSelection
          )) +
            geom_point(size = 1.0) +
            theme_bw() +
            # scale_color_manual(name='Cluster selection',
            #                    breaks=c('selected phenotype', 'other clusters'),
            #                    values=c('selected phenotype'='green', 'other clusters'='grey'))
            theme(legend.text = element_text(size =
                                               8)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
        )
    })

    # Observe node update picker ----
    observeEvent(input$updateNodePicker, {

      node <- reactVals$graph$nodes %>% filter(label == input$updateNodePicker)
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

    # Delete Node ----
    observeEvent(input$deleteNodeBtn, {

      node <- reactVals$graph$nodes %>% filter(label == input$updateNodePicker)

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
        if( dim(df01Tree[df01Tree$cell == node$label,])[1] >0 ) {

          df01Tree[df01Tree$cell == node$label,]$cell <<- parent_label

        }

        annotationlist[annotationlist == node$label] <<- NULL

        updatePickerInput(
          session,
          inputId = "parentPicker",
          selected = NULL,
          choices = annotationlist
        )
        updatePickerInput(
          session,
          inputId = "updateNodePicker",
          selected = NULL,
          choices = annotationlist
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


    importNodes <- function(id, label, pm, nm, parent_id) {

      if(id == parent_id) {return()}

      parent_row <- reactVals$graph$nodes %>% filter(id == parent_id)
      parent_name <- parent_row$label

      add_node(reactVals$graph, parent_name, label,)
    }

    # Load Annotation Tree ----
    observeEvent(input$btnImportTree, {

      req(input$fNodes)
      req(input$fEdges)
      req(input$fAnno)

      req(input$fMarkerExpr)
      req(input$cluster_freq)

      df_nodes <- read.csv(input$fNodes$datapath)
      df_edges <- read.csv(input$fEdges$datapath)
      df_anno <- read.csv(input$fAnno$datapath)

      df_nodes$pm[is.na(df_nodes$pm)] <- ""
      df_nodes$pm <- strsplit(df_nodes$pm, "\\|")

      df_nodes$nm[is.na(df_nodes$nm)] <- ""
      df_nodes$nm <- strsplit(df_nodes$nm, "\\|")

      reactVals$graph$nodes <- df_nodes
      reactVals$graph$edges <<- df_edges

      df01Tree <<- df_anno

      annotationlist <<- as.list(df_nodes$label)

    })

    # Export Annotation Tree ----
    output$exportAnnotationBtn <- downloadHandler(

      filename = function(){
        paste("cellcat_annotation_data_", Sys.Date(), ".zip", sep = "")
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

        download_list <- list(annTable = df01Tree["cell"],
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

    ##
    # Server - Thresholds Tab -------------------------------------------------
    ##
    output$table = renderDataTable(
      reactVals$th,
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

      selRow <- reactVals$th[input$table_rows_selected,]
      marker <- reactVals$th[input$table_rows_selected, 1]
      # threshold value for vertical line
      myTH <- reactVals$th[input$table_rows_selected, 2]

      marker_expr <- getMarkerDistDF(marker, input$radio)
      myRenderFunction(marker_expr, myTH)
    })

    # update Tree Annotation ----
    ## update the tree after a change in the thresholds
    ## walk through the list of nodes and adjust the annotation if
    ## needed based on the filtered DF
    updateTreeAnnotation <- function(df_row) {

      df01Tree$cell <- "Unassigned"

      # for all rows in annotation table
      # get all parents for a row:
      reactVals$graph$nodes$to <- reactVals$graph$edges$to

      # Iterate through the dataframe row by row
      for (i in 1:nrow(reactVals$graph$nodes)) {

        posMarker <- list()
        negMarker <- list()

        nodeID <- reactVals$graph$nodes$id[i]

        # now for that nodeID, get all parents and collect
        # the pm and nm markers
        while (nodeID > 1) {

          posMarker <- c(posMarker, unlist(reactVals$graph$nodes$pm[reactVals$graph$nodes$id == nodeID]))
          negMarker <- c(negMarker, unlist(reactVals$graph$nodes$nm[reactVals$graph$nodes$id == nodeID]))

          nodeID <- reactVals$graph$edges$to[reactVals$graph$edges$from == nodeID]

        }
        # now we have all positive and negative marker for that type collected
        # and we can start filtering the df
        # remove the empty strings in the markers
        posMarker <- posMarker[nzchar(posMarker)]
        negMarker <- negMarker[nzchar(negMarker)]

        tmp <- filterHM(df01Tree,unlist(posMarker), unlist(negMarker), reactVals$th)

        df01Tree[rownames(tmp), 'cell'] <<- reactVals$graph$node$label[i]

        updateClusterLabels(tmp)
      }
    }

    ## Event on click scatter plot for setting the vertical line
    ##
    # observe click scatterplot -----
    ##
    observeEvent(input$plot_click, {

      req(input$table_rows_selected)
      selRow <- reactVals$th[input$table_rows_selected,]
      marker <- reactVals$th[input$table_rows_selected, 1]

      if (input$radio == "1") {
        reactVals$th[input$table_rows_selected, 2] <- round(input$plot_click$x, 3)
        myTH <- reactVals$th[input$table_rows_selected, 2]
        marker_expr <- getMarkerDistDF(marker, input$radio)
        myRenderFunction(marker_expr, myTH)
      }
      if (input$radio == "2") {
        reactVals$th[input$table_rows_selected, 2] <- round(10^(input$plot_click$x), 6)
        myTH <- reactVals$th[input$table_rows_selected, 2]
        marker_expr <- getMarkerDistDF(marker, input$radio)
        myRenderFunction(marker_expr, myTH)
      }

      updateTreeAnnotation(reactVals$th[input$table_rows_selected])
    })

    ##
    ## transform according the selection event --------------------------------
    ##
    observeEvent(input$radio, {
      req(input$table_rows_selected)
      selRow <- reactVals$th[input$table_rows_selected,]
      marker <- reactVals$th[input$table_rows_selected, 1]
      myTH <- reactVals$th[input$table_rows_selected, 2]
      if ( input$radio == "1") {
        marker_expr <- getMarkerDistDF(marker, input$radio)
        myRenderFunction(marker_expr, myTH)
      } else if (input$radio == "2") {
        marker_expr <- getMarkerDistDF(marker, input$radio)
        myRenderFunction(marker_expr, myTH)
      }
    })

    ##
    # Thresholds plotting function ----
    ##
    myRenderFunction <- function(me, myTH){

      output$plot <- renderPlot({
        set.seed(1)
        if (input$radio == "1") {
          ggplot(me, aes_string(x=me[,1], y=me[,2])) +
            geom_point(size=1) +
            theme(axis.title.y = element_blank(),
                  axis.ticks.y  = element_blank(),
                  axis.text.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank()) +
            labs(x = "Scale 0 to 1") +
            geom_vline(xintercept = myTH, linetype="dotted",
                       color = "blue", size=1.5)

        } else if (input$radio == "2") {
          ggplot(me, aes_string(x=me[,1], y=me[,2])) +
            geom_point(size=1) +
            theme(axis.title.y = element_blank(),
                  axis.ticks.y  = element_blank(),
                  axis.text.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank()) +
            labs(x = "log 10") +
            geom_vline(xintercept = log10(myTH), linetype="dotted",
                       color = "blue", size=1.5)
        }
      })
      #
      output$plot2 <- renderPlot({
        if (input$radio == "1") {
          ggplot(me, aes_string(x = me[, 1])) +
            geom_histogram(bins = 80) +
            labs(x = "Scale 0 to 1") +
            geom_vline(
              xintercept = myTH,
              linetype = "dotted",
              color = "blue",
              size = 1.5
            )

        } else if (input$radio == "2") {
          ggplot(me, aes_string(x = me[, 1])) +
            geom_histogram(bins = 80) +
            labs(x = "log 10") +
            geom_vline(
              xintercept = log10(myTH),
              linetype = "dotted",
              color = "blue",
              size = 1.5
            )

        } else {
          ggplot(me, aes_string(x = me[, 1])) +
            geom_histogram(bins = 80) +
            geom_vline(
              xintercept = log2(myTH),
              linetype = "dotted",
              color = "blue",
              size = 1.5
            )
        }
      })
    }

    ##
    # Server - Annotations Tab ------------------------------------------------
    ##
    ## upload the annotation file ---------------------------------------------
    observeEvent(input$annTable,{
      annData <<- read.csv(input$annTable$datapath)
      annData$X <- NULL
      at$data <- annData
    })

    ## upload the marker threshold file ---------------------------------------
    observeEvent(input$fTH,{
      th <<- read.csv(input$fTH$datapath)
      th$X <- NULL
      reactVals$th<- th
    })

    ## upload the metadata ----------------------------------------------------
    observeEvent(input$metadata,{
      md <<- read.csv(input$metadata$datapath)
      md$X <- NULL
      reactVals$md<- md
    })

    ## upload the counts talbe -------------------------------------------------
    observeEvent(input$counts_table,{
      ct <<- read.csv(input$counts_table$datapath)
      ct$X <- NULL
      reactVals$counts_table <- ct
    })

    ## create Phenotype name from picker e.g. CD4+CD8+CD19 --------------------
    observeEvent(input$myPickerPos, {
      output$r1 <- renderText(setPhenotypeName(input$myPickerPos, "pos", ph_name))
    })
    observeEvent(input$myPickerNeg, {
      output$r2 <- renderText(setPhenotypeName(input$myPickerNeg, "neg", ph_name))
    })

    ## Annotation Table -------------------------------------------------------
    observeEvent(c(input$myPickerPos, input$myPickerNeg, reactVals$th),{
      req(input$fMarkerExpr)
      req(input$cluster_freq)

      tmp=filterHM(df01,input$myPickerPos, input$myPickerNeg, reactVals$th)

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

    ##
    # Annotation Tree ---------------------------------------------------------
    ##
    output$hm_tree <- renderPlot({
      my_new_hm <- dr_umap
      # my_new_hm <- my_new_hm %>% normalize01()
      if (dim(my_new_hm)[1] > 0) {
        # my_new_hm %>% select(-c("cluster_number", "u1", "u2"))
        pheatmap(my_new_hm %>% select(-c("cluster_number", "u1", "u2")), cluster_cols = F)
      } else {
        ggplot() + theme_void() + ggtitle("Select area on Umap to plot Heatmap")
      }
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

    # Server - UMAP interactive Tab -------------------------------------------
    output$umap2 <- renderPlot(
      ggplot(dr_umap, aes(x = u1, y = u2)) +
        geom_point(size = 1.0) +
        theme_bw() +
        theme(legend.text=element_text(size=8)) +
        guides(color = guide_legend(override.aes = list(size = 4)))
    )

    output$hm2 <- renderPlot({
      my_new_hm <- round(brushedPoints(dr_umap, input$umap2_brush), 2)

      if (dim(my_new_hm)[1] > 0) {

        pheatmap(my_new_hm %>% select(-c("cluster_number", "u1", "u2")), cluster_cols = F)
      } else {
        ggplot() + theme_void() + ggtitle("Select area on Umap to plot Heatmap")
      }
    })

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

    # Server - UMAP Marker Expression Tab ----
    output$umap3 <-
      renderPlot(
        ggplot(dr_umap, aes_string(
          x = "u1",
          y = "u2",
          color = input$markerSelect
        )) +
          geom_point(size = 1.0) +
          theme_bw() +
          scale_color_gradientn(input$markerSelect,
                                colours = colorRampPalette(rev(brewer.pal(
                                  n = 11, name = "Spectral"
                                )))(50))
      )

    # Server - Differential Abundance Tab ----
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

    # Do the differential abundance ----
    observeEvent(input$doDA, {
      # browser()
      req(reactVals$counts_table, reactVals$md)

      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df01Tree$cell
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

      reactVals$DA_result_table <- DA_df

    })

    # output$mynetworkid <- renderVisNetwork({
    #   visNetwork(graph$nodes, graph$edges, width = "100%") %>%
    #     visEdges(arrows = "from") %>%
    #     visHierarchicalLayout()
    # })

    ## -------------------------------------------------------------------
    # Interactive DA Tree ----
    # output$interactiveTree <- renderVisNetwork({
    #   visNetwork(graph$nodes, graph$edges, width = "100%") %>%
    #     visEdges(arrows = "from") %>%
    #     visHierarchicalLayout()
    # })

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
      countsTable['cell'] <- df01Tree$cell
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
      countsTable['cell'] <- df01Tree$cell
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
        geom_boxplot() +
        xlab("Condition") +
        ylab("Values") +
        ggtitle(reactVals$graph$nodes$label[myNode$selected])


    })

    # children <- graph$edges$from[graph$edges$to == myNode$selected]

    observeEvent(input$btnLoadDemoData, {


      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "loading Data...", value = 0)

      # browser()

      ## Load median expression and cell frequencies
      df <- read.csv("data/McCarthy_expr_median_400.csv")
      df_global <<- df

      # create initial master node of all Unassigned clusters
      nodes <- tibble(id = 1,
                      label = "Unassigned",
                      pm = list(""),
                      nm = list(""),
                      color = "blue"
      )

      edges <- data.frame(from = numeric(), to = numeric())
      reactVals$graph <- list(nodes = nodes, edges = edges)

      annotationlist <<- list("Unassigned")

      # cell_freq <<- read.csv("data/McCarthy_cluster_freq_400.csv") %>%
      #   mutate(total = sum(column_name)) %>%
      #   mutate(frequency = round(column_name / total * 100, 2))

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

      annotaionDF <<- data.frame("cell" = "unassigned",
                                 clusterSize = cell_freq$clustering_prop)
      at <<- reactiveValues(data = annotaionDF, dr_umap = dr_umap)

      reactVals$th <- kmeansTH(df01)


    })
    ## -------------------------------------------------------------------
    # Upload All Annotation Demo Data ----
    observeEvent(input$btnLoadAnnoData, {
      # browser()
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())

      progress$set(message = "loading Data...", value = 0)

      # browser()

      ## Load median expression and cell frequencies
      df <- read.csv("data/McCarthy_expr_median_400.csv")
      df_global <<- df

      # create initial master node of all Unassigned clusters
      nodes <- tibble(id = 1,
                      label = "Unassigned",
                      pm = list(""),
                      nm = list(""),
                      color = "blue"
      )

      edges <- data.frame(from = numeric(), to = numeric())
      reactVals$graph <- list(nodes = nodes, edges = edges)

      annotationlist <<- list("Unassigned")

      cell_freq <<- read.csv("data/McCarthy_cluster_freq_400.csv")

      # cell_freq <<- read.csv("data/McCarthy_cluster_freq_400.csv") %>%
      #   mutate(total = sum(column_name)) %>%
      #   mutate(frequency = round(column_name / total * 100, 2))

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

      annotaionDF <<- data.frame("cell" = "unassigned",
                                 clusterSize = cell_freq$clustering_prop)
      at <<- reactiveValues(data = annotaionDF, dr_umap = dr_umap)

      reactVals$th <- kmeansTH(df01)


      ## Load metadata
      md <<- read.csv("data/metadata_Fig1_Panel2.csv")
      md$X <- NULL
      reactVals$md<- md

      ## Load counts table
      ct <<- read.csv("data/props_table_400.csv")
      ct$X <- NULL
      reactVals$counts_table <- ct

      ## Load annotaiton Tree
      df_nodes <- read.csv("data/nodesTable_data.csv")
      df_edges <- read.csv("data/edgesTable_data.csv")
      df_anno <- read.csv("data/annTable_data.csv")

      df_nodes$pm[is.na(df_nodes$pm)] <- ""
      df_nodes$pm <- strsplit(df_nodes$pm, "\\|")

      df_nodes$nm[is.na(df_nodes$nm)] <- ""
      df_nodes$nm <- strsplit(df_nodes$nm, "\\|")

      reactVals$graph$nodes <- df_nodes
      reactVals$graph$edges <- df_edges

      df01Tree <<- df_anno

      annotationlist <<- as.list(df_nodes$label)

    })

  }

  shinyApp(
    ui = ui,
    server = server
  )
}

