# Annotation tab --------------------------------------------------------------
# Build the annotation tree: pick a parent node and positive / negative markers,
# preview the matching clusters, and create or delete nodes. Also holds the
# CATALYST meta-cluster level and the CATALYST merge / export.

mod_annotation_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 400, open = "always", title = "Build a phenotype",
      tags$div(class = "section-label", "Parent node"),
      pickerInput(ns("parentPicker"), label = NULL, choices = NULL,
                  options = picker_options(), multiple = FALSE),
      uiOutput(ns("lineage")),
      tags$div(class = "section-label", "Markers"),
      marker_selector_input(ns("markers")),
      uiOutput(ns("summary")),
      textInput(ns("newNode"), "Phenotype name", placeholder = "e.g. CD4 T cells"),
      actionButton(ns("createNodeBtn"), "Create node", icon = icon("plus"),
                   class = "btn-primary w-100"),
      shinyjs::hidden(tags$div(
        id = ns("catalyst_panel"),
        tags$hr(),
        tags$div(class = "section-label", "CATALYST"),
        pickerInput(ns("metaLevel"), label = "Meta-cluster level", choices = NULL,
                    options = picker_options(), multiple = FALSE),
        tags$div(class = "d-flex gap-2",
                 actionButton(ns("mergeCatalystBtn"), "Merge annotation",
                              class = "btn-outline-primary btn-sm"),
                 downloadButton(ns("exportCatalystBtn"), "Export object",
                                class = "btn-outline-secondary btn-sm"))
      ))
    ),
    uiOutput(ns("notice")),
    card(
      full_screen = TRUE,
      card_title("Annotation tree",
                 actionButton(ns("deleteNodeBtn"), "Delete node", icon = icon("trash"),
                              class = "btn-outline-danger btn-sm"),
                 downloadButton(ns("exportAnnotationBtn"), "Export annotation",
                                class = "btn-outline-secondary btn-sm")),
      card_body(class = "p-1", visNetworkOutput(ns("tree"), height = "500px"))
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(full_screen = TRUE, card_header("Heatmap of the selection"),
           plotOutput(ns("heatmap"), height = "480px")),
      card(full_screen = TRUE, card_header("Selection on the UMAP"),
           plotOutput(ns("umap"), height = "480px"))
    )
  )
}

mod_annotation_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {

    show_error <- function(title, msg) {
      showModal(modalDialog(title = title, msg, easyClose = TRUE, footer = modalButton("OK")))
    }

    output$notice <- renderUI(if (is.null(state$expr)) no_data_notice())

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

    lineage <- reactive(lineage_markers(state$graph, current_node()$id))

    output$lineage <- renderUI({
      node <- current_node()
      path <- state$graph$nodes$label[match(node_lineage(state$graph, node$id), state$graph$nodes$id)]
      def <- lineage()
      tags$div(
        class = "lineage",
        tags$div(lapply(seq_along(path), function(i) {
          tagList(if (i > 1) tags$span(class = "sep", "\u203a"),
                  if (i == length(path)) tags$strong(path[i]) else path[i])
        })),
        if (length(def$pos) + length(def$neg) > 0) tags$div(marker_chips(def$pos, def$neg))
      )
    })

    # Marker selector ---------------------------------------------------------
    # Markers already defined by the lineage are locked.
    observeEvent(state$markers, {
      flagged <- if (!is.null(state$th)) state$th$cell[state$th$color == "red"] else character(0)
      update_marker_selector(session, "markers", markers = state$markers, flagged = flagged)
    }, ignoreNULL = FALSE)

    observeEvent(state$th, {
      update_marker_selector(session, "markers", flagged = state$th$cell[state$th$color == "red"],
                             clear = FALSE)
    })

    observeEvent(list(state$markers, current_node()), {
      def <- lineage()
      locked <- c(setNames(rep("pos", length(def$pos)), def$pos),
                  setNames(rep("neg", length(def$neg)), def$neg))
      update_marker_selector(session, "markers", locked = locked)
    })

    picked <- reactive(marker_selection(input$markers))

    # Clusters of the selected node filtered by the picked markers ------------
    preview <- reactive({
      req(state$expr, state$th, selected_label())
      expr <- state$expr
      filterHM(expr[expr$cell == selected_label(), state$markers, drop = FALSE],
               picked()$pos, picked()$neg, state$th)
    })

    output$summary <- renderUI({
      req(state$expr)
      n_parent <- sum(state$expr$cell == selected_label())
      tags$div(
        class = "selection-summary my-2",
        tags$div(tags$div(class = "value", nrow(preview())),
                 tags$div(class = "label", sprintf("of %d clusters in parent", n_parent))),
        tags$div(tags$div(class = "value",
                          sprintf("%.2f%%", sum(state$expr[rownames(preview()), "freq"]))),
                 tags$div(class = "label", "of all cells"))
      )
    })

    # the plots follow the marker toggles with a short delay
    preview_plot <- debounce(preview, 400)

    output$heatmap <- renderPlot({
      req(state$expr)
      plot_cluster_heatmap(preview_plot())
    })

    output$umap <- renderPlot({
      req(state$umap)
      plot_umap_selection(state$umap, filterColor(state$expr, preview_plot()))
    })

    # Tree --------------------------------------------------------------------
    output$tree <- renderVisNetwork({
      req(state$graph, state$expr)
      tree_network(state$graph, session$ns("tree_click"), expr = state$expr,
                   selected = isolate(current_node()$id))
    })

    observeEvent(current_node(), {
      visNetworkProxy(session$ns("tree")) %>%
        visUpdateNodes(tree_selection_style(state$graph$nodes, state$expr, current_node()$id))
    })

    # Create node -------------------------------------------------------------
    observeEvent(input$createNodeBtn, {
      req(state$expr)
      name <- trimws(input$newNode)
      parent <- selected_label()

      if (length(picked()$pos) + length(picked()$neg) == 0) {
        return(show_error("No marker selection", "Select positive and / or negative markers."))
      }
      if (name == "") {
        return(show_error("Phenotype name", "Set a name for this phenotype."))
      }
      if (name %in% state$graph$nodes$label) {
        return(show_error("Name taken", sprintf("A phenotype named '%s' already exists.", name)))
      }
      selection <- preview()
      if (nrow(selection) == 0) {
        return(show_error("Empty selection", "No cluster of the parent matches this marker selection."))
      }

      state$expr[rownames(selection), "cell"] <- name
      selected_label(name)
      state$graph <- add_node(state$graph, parent, name,
                              list(picked()$pos), list(picked()$neg), "blue")
      updateTextInput(session, "newNode", value = "")
      showNotification(sprintf("Created '%s' with %d clusters.", name, nrow(selection)),
                       type = "message")
    })

    # Delete node -------------------------------------------------------------
    observeEvent(input$deleteNodeBtn, {
      node <- current_node()

      if (node_has_children(state$graph, node$id)) {
        return(show_error("Cannot delete node",
                          sprintf("'%s' has child nodes. Only leaf nodes can be deleted.", node$label)))
      }
      parent_id <- state$graph$edges$to[state$graph$edges$from == node$id]
      showModal(modalDialog(
        title = "Delete node",
        sprintf("Delete '%s'? Its %d clusters are returned to '%s'.", node$label,
                sum(state$expr$cell == node$label),
                state$graph$nodes$label[state$graph$nodes$id == parent_id]),
        easyClose = TRUE,
        footer = tagList(modalButton("Cancel"),
                         actionButton(session$ns("confirmDelete"), "Delete", class = "btn-danger"))
      ))
    })

    observeEvent(input$confirmDelete, {
      removeModal()
      node <- current_node()
      graph <- state$graph
      req(!node_has_children(graph, node$id))

      # give the clusters of the node back to its parent
      parent_id <- graph$edges$to[graph$edges$from == node$id]
      parent_label <- graph$nodes$label[graph$nodes$id == parent_id]
      state$expr$cell[state$expr$cell == node$label] <- parent_label

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
      shinyjs::toggle("catalyst_panel", condition = !is.null(state$sce))
      if (!is.null(state$sce)) {
        updatePickerInput(session, "metaLevel",
                          choices = catalyst_meta_levels(state$sce),
                          selected = state$meta_level)
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
