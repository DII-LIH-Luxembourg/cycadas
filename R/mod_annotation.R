# Tree-Annotation tab ---------------------------------------------------------
# Build the annotation tree: pick a parent node and positive / negative markers,
# preview the matching clusters, and create or delete nodes. Also holds the
# CATALYST meta-cluster level and the CATALYST merge / export.

mod_annotation_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      width = 4,
      box(
        width = NULL,
        title = "Create Node",
        textOutput(ns("selection_freq")),
        textOutput(ns("selection_n")),
        textInput(ns("newNode"), "Set Name..."),
        pickerInput(ns("parentPicker"), label = "Select Parent:", choices = NULL,
                    options = picker_options(), multiple = F)
      ),
      box(
        width = NULL,
        title = "Select MetaCluster Level (CATALYST)",
        div(id = ns("metadiv"),
            pickerInput(ns("metaLevel"), label = "Select Cluster Level:", choices = NULL,
                        options = picker_options(), multiple = F))
      ),
      box(
        width = NULL,
        column(width = 4,
               checkboxGroupButtons(ns("treePickerPos"), label = "Positive:",
                                    choices = c("A"), direction = "vertical")),
        column(width = 4,
               checkboxGroupButtons(ns("treePickerNeg"), label = "Negative:",
                                    choices = c("A"), direction = "vertical"))
      ),
      box(width = NULL, title = "Create New Node",
          actionButton(ns("createNodeBtn"), "Create Node")),
      box(width = NULL, title = "Delete Node",
          actionButton(ns("deleteNodeBtn"), "Delete Node")),
      box(width = NULL, title = "Export Annotation",
          downloadButton(ns("exportAnnotationBtn"), "Export Annotation")),
      box(width = NULL, title = "Merge CATALYST Annotation",
          actionButton(ns("mergeCatalystBtn"), "Merge Annotation")),
      box(width = NULL, title = "Export CATALYST Annotation",
          downloadButton(ns("exportCatalystBtn"), "Export Annotation"))
    ),
    column(
      width = 8,
      box(width = NULL, title = "Annotation Tree", visNetworkOutput(ns("tree"))),
      box(width = NULL, title = "Heatmap", plotOutput(ns("heatmap"))),
      box(width = NULL, plotOutput(ns("umap")))
    )
  )
}

mod_annotation_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    show_error <- function(title, msg) {
      showModal(modalDialog(title = title, msg, easyClose = TRUE, footer = NULL))
    }

    # Selected node -----------------------------------------------------------
    # `selected_label` is the source of truth; the parent picker and clicks on
    # the tree both update it.
    selected_label <- reactiveVal(NULL)

    select_node <- function(label) {
      selected_label(label)
      updatePickerInput(session, "parentPicker", selected = label)
    }

    labels <- reactive({
      req(state$graph)
      state$graph$nodes$label
    })

    observeEvent(labels(), {
      sel <- isolate(selected_label())
      if (is.null(sel) || !sel %in% labels()) {
        sel <- labels()[1]
        selected_label(sel)
      }
      updatePickerInput(session, "parentPicker", choices = labels(), selected = sel)
    })

    observeEvent(state$graph, {
      if (is.null(state$graph)) {
        selected_label(NULL)
        updatePickerInput(session, "parentPicker", choices = character(0))
      }
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    observeEvent(input$parentPicker, selected_label(input$parentPicker))

    observeEvent(input$tree_click, {
      req(length(input$tree_click) > 0)
      nodes <- state$graph$nodes
      select_node(nodes$label[nodes$id == input$tree_click[1]])
    })

    current_node <- reactive({
      req(state$graph, selected_label())
      node <- state$graph$nodes[state$graph$nodes$label == selected_label(), ]
      req(nrow(node) == 1)
      node
    })

    # Marker pickers ----------------------------------------------------------
    # Markers already used by the selected node are disabled, and a marker
    # picked on one side is disabled on the other.
    update_marker_picker <- function(inputId, selected, disabled) {
      updateCheckboxGroupButtons(session, inputId = inputId, choices = state$markers,
                                 selected = selected, disabledChoices = disabled)
    }

    observeEvent(list(state$markers, current_node()), {
      node <- current_node()
      update_marker_picker("treePickerPos", NULL, unlist(node$pm))
      update_marker_picker("treePickerNeg", NULL, unlist(node$nm))
    })

    observeEvent(input$treePickerPos, {
      req(length(state$markers) > 0)
      update_marker_picker("treePickerNeg", input$treePickerNeg,
                           union(unlist(current_node()$nm), input$treePickerPos))
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    observeEvent(input$treePickerNeg, {
      req(length(state$markers) > 0)
      update_marker_picker("treePickerPos", input$treePickerPos,
                           union(unlist(current_node()$pm), input$treePickerNeg))
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    # Clusters of the selected node filtered by the picked markers ------------
    preview <- reactive({
      req(state$expr, state$th, selected_label())
      expr <- state$expr
      filterHM(expr[expr$cell == selected_label(), state$markers, drop = FALSE],
               input$treePickerPos, input$treePickerNeg, state$th)
    })

    # the plots follow the marker toggles with a short delay
    preview_plot <- debounce(preview, 400)

    output$selection_freq <- renderText({
      paste0(round(sum(state$expr[rownames(preview()), "freq"]), 3), "% in selection")
    })
    output$selection_n <- renderText({
      paste0(nrow(preview()), " Cluster in selection")
    })

    output$heatmap <- renderPlot({
      if (is.null(state$expr)) plot_message("No Data Available") else plot_cluster_heatmap(preview_plot())
    })

    output$umap <- renderPlot({
      req(state$umap)
      plot_umap_selection(state$umap, filterColor(state$expr, preview_plot()))
    })

    output$tree <- renderVisNetwork({
      req(state$graph, state$expr)
      graph <- state$graph
      graph$nodes <- color_nodes_by_usage(graph$nodes, state$expr$cell)
      tree_network(graph, session$ns("tree_click"), hierarchical = TRUE)
    })

    # Create node -------------------------------------------------------------
    observeEvent(input$createNodeBtn, {
      req(state$expr)
      name <- input$newNode
      parent <- selected_label()

      if (is.null(input$treePickerPos) && is.null(input$treePickerNeg)) {
        return(show_error("No Marker Selection", "Select positive and / or negative markers!"))
      }
      if (name == "") {
        return(show_error("Phenotype Name", "Set a Name for this Phenotype!"))
      }
      if (name %in% state$graph$nodes$label) {
        return(show_error("Naming Error!", "The Name for this Phenotype is already taken!"))
      }
      selection <- preview()
      if (nrow(selection) == 0) {
        return(show_error("No Result", "The Selection for this Phenotype is empty!"))
      }

      state$expr[rownames(selection), "cell"] <- name
      selected_label(name)
      state$graph <- add_node(state$graph, parent, name,
                              list(input$treePickerPos), list(input$treePickerNeg), "blue")
    })

    # Delete node -------------------------------------------------------------
    observeEvent(input$deleteNodeBtn, {
      node <- current_node()
      graph <- state$graph

      if (node_has_children(graph, node$id)) {
        return(show_error("Cannot delete Node", "The selected Node is not a Leaf Node!"))
      }

      # give the clusters of the node back to its parent
      parent_id <- graph$edges$to[graph$edges$from == node$id]
      parent_label <- graph$nodes$label[graph$nodes$id == parent_id]
      state$expr$cell[state$expr$cell == node$label] <- parent_label

      updateTextInput(session, "newNode", value = "")
      selected_label(parent_label)
      state$graph <- delete_leaf_node(graph, node$id)
    })

    # Export annotation -------------------------------------------------------
    output$exportAnnotationBtn <- downloadHandler(
      filename = function() {
        paste("cycadas_annotation_data_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file) {
        export_annotation_zip(file, state$expr, state$graph)
      },
      contentType = "application/zip"
    )

    # CATALYST ----------------------------------------------------------------
    observeEvent(state$sce, {
      if (is.null(state$sce)) {
        shinyjs::disable("metadiv")
        updatePickerInput(session, "metaLevel", choices = character(0))
      } else {
        updatePickerInput(session, "metaLevel",
                          choices = catalyst_meta_levels(state$sce),
                          selected = state$meta_level)
        shinyjs::enable("metadiv")
      }
    }, ignoreNULL = FALSE)

    observeEvent(input$metaLevel, {
      req(state$sce)
      # the picker also fires when its choices are set after an import
      if (identical(input$metaLevel, state$meta_level)) return()

      state$meta_level <- input$metaLevel
      load_dataset(state, catalyst_meta_level_data(state$sce, input$metaLevel, state$markers))
    })

    observeEvent(input$mergeCatalystBtn, {
      req(state$sce, state$expr)
      state$sce <- merge_catalyst_annotation(state$sce, state$expr, state$meta_level)
      showNotification("Annotation merged into the CATALYST object.", type = "message")
    })

    output$exportCatalystBtn <- downloadHandler(
      filename = function() {
        paste("Annotated_sce_", Sys.Date(), ".rds", sep="")
      },
      content = function(file) {
        saveRDS(state$sce, file)
      },
      contentType = "application/octet-stream"
    )
  })
}
