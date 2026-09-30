# UI helpers ------------------------------------------------------------------

# truthy fileInput?
has_file <- function(x) !is.null(x) && !is.null(x$datapath) && nzchar(x$datapath)

# JavaScript and CSS in inst/app/www
cycadas_dependency <- function() {
  htmltools::htmlDependency(
    "cycadas", as.character(utils::packageVersion("cycadas")),
    src = c(file = app_file("www")),
    script = "cycadas.js", stylesheet = "cycadas.css"
  )
}

# Status light with a label ----
status_item <- function(ok, label, detail = NULL) {
  tags$li(
    tags$span(tags$span(class = paste("status-dot", if (ok) "ok" else "missing")), label),
    if (!is.null(detail)) tags$span(class = "detail", detail)
  )
}

# Status lights for a set of required files plus a summary line
file_status_ui <- function(present, labels, ready_msg, missing_msg) {
  tagList(
    tags$ul(class = "status-list", Map(status_item, present, labels)),
    tags$p(class = paste("small mt-2 mb-0", if (all(present)) "text-success" else "text-body-secondary"),
           if (all(present)) ready_msg else missing_msg)
  )
}

# Track which file inputs hold a file that has not been imported yet ----------
# A fileInput cannot be cleared on the server, so `input[[id]]` keeps its value
# after shinyjs::reset(). This tracker is set when a file is chosen and cleared
# by `$clear()`. Must be called inside a module server.
file_tracker <- function(input, ids) {
  present <- do.call(reactiveValues, setNames(as.list(rep(FALSE, length(ids))), ids))
  for (id in ids) {
    local({
      id <- id
      observeEvent(input[[id]], present[[id]] <- has_file(input[[id]]))
    })
  }
  list(
    present = function() vapply(ids, function(id) isTRUE(present[[id]]), logical(1)),
    ready = function() all(vapply(ids, function(id) isTRUE(present[[id]]), logical(1))),
    clear = function() for (id in ids) present[[id]] <- FALSE
  )
}

# Dropdown options shared by the pickers
picker_options <- function(...) {
  list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3",
       `live-search` = TRUE, ...)
}

# Card header with a title and optional actions on the right
card_title <- function(title, ...) {
  card_header(class = "d-flex align-items-center gap-2", tags$span(title),
              if (length(list(...))) tags$span(class = "card-actions ms-auto", ...))
}

# Shown on analysis tabs before any data is loaded ----
no_data_notice <- function() {
  tags$div(
    class = "no-data",
    icon("database"),
    tags$h5("No data loaded"),
    tags$p(class = "text-body-secondary",
           "Import cluster data or load the demo data on the Workspace tab.")
  )
}

# Marker selector input -------------------------------------------------------
# Value: list(pos = character, neg = character). Filled by
# update_marker_selector(); the behaviour is in inst/app/www/cycadas.js.
marker_selector_input <- function(id) {
  tags$div(
    id = id, class = "marker-selector",
    tags$div(class = "marker-tools",
             tags$input(type = "search", class = "form-control form-control-sm marker-search",
                        placeholder = "Filter markers"),
             tags$span(class = "marker-count"),
             tags$button(type = "button", class = "btn btn-link marker-clear", hidden = NA, "Clear")),
    tags$p(class = "marker-empty small text-body-secondary", "No markers loaded."),
    tags$div(class = "marker-grid"),
    tags$div(class = "marker-legend",
             HTML("&minus; below threshold &middot; + above threshold &middot; "),
             tags$span(class = "dot", HTML("&#9679;")), " not bimodal")
  )
}

# `locked` is a named character vector: marker = "pos" / "neg"
update_marker_selector <- function(session, id, markers = NULL, locked = NULL,
                                   flagged = NULL, clear = TRUE) {
  msg <- list(clear = clear)
  if (!is.null(markers)) msg$markers <- I(markers)
  if (!is.null(locked)) msg$locked <- if (length(locked)) as.list(locked) else setNames(list(), character(0))
  if (!is.null(flagged)) msg$flagged <- I(flagged)
  session$sendInputMessage(id, msg)
}

# Selected markers as list(pos, neg) with character(0) for none
marker_selection <- function(value) {
  list(pos = as.character(unlist(value$pos)), neg = as.character(unlist(value$neg)))
}

# Marker definition chips, e.g. CD3+ CD8-
marker_chips <- function(pos, neg) {
  tagList(
    lapply(pos, function(m) tags$span(class = "marker-chip pos", paste0(m, "+"))),
    lapply(neg, function(m) tags$span(class = "marker-chip neg", paste0(m, "\u2212")))
  )
}

# Zoom / fit / pan / export buttons for a tree widget --------------------------
# Place inside a container together with visNetworkOutput(widget_id).
tree_toolbar <- function(widget_id) {
  btn <- function(action, icon_name, title) {
    tags$button(type = "button", class = "btn btn-light btn-sm", title = title,
                `aria-label` = title,
                onclick = sprintf("cycadasTree.run('%s', '%s')", widget_id, action),
                icon(icon_name))
  }
  tags$div(
    class = "tree-toolbar",
    tags$div(class = "btn-group",
             btn("zoomOut", "magnifying-glass-minus", "Zoom out"),
             btn("zoomIn", "magnifying-glass-plus", "Zoom in")),
    tags$div(class = "btn-group",
             btn("fit", "expand", "Fit whole tree"),
             btn("reset", "rotate-left", "Reset view"),
             btn("focus", "crosshairs", "Center on selected node")),
    tags$div(class = "btn-group",
             btn("left", "arrow-left", "Move left"),
             btn("up", "arrow-up", "Move up"),
             btn("down", "arrow-down", "Move down"),
             btn("right", "arrow-right", "Move right")),
    tags$div(class = "btn-group",
             tags$button(type = "button", class = "btn btn-light btn-sm",
                         title = "Download the whole tree as SVG (vector graphic)",
                         onclick = sprintf("cycadasTree.exportTree('%s', 'svg')", widget_id),
                         icon("download"), "SVG"),
             tags$button(type = "button", class = "btn btn-light btn-sm",
                         title = "Download the whole tree as PNG (3x resolution)",
                         onclick = sprintf("cycadasTree.exportTree('%s', 'png')", widget_id),
                         "PNG"))
  )
}

# visNetworkOutput with the toolbar on top
tree_output <- function(widget_id, height) {
  tags$div(class = "tree-container",
           tree_toolbar(widget_id),
           visNetworkOutput(widget_id, width = "100%", height = height))
}

# Phenotype composition list ---------------------------------------------------
# Indented rows in tree order with a bar for the total share (light) and the
# remaining share (dark). Clicking a row sets `input_id` to the node id.
composition_list <- function(comp, selected = NULL, input_id) {
  rows <- lapply(seq_len(nrow(comp)), function(i) {
    r <- comp[i, ]
    tags$div(
      class = paste("comp-row", if (isTRUE(r$id == selected)) "selected"),
      `data-id` = r$id,
      title = sprintf("%s: %.2f%% of cells, remaining %.2f%%, %d clusters",
                      r$phenotype, r$total, r$remaining, r$clusters_total),
      onclick = sprintf("Shiny.setInputValue('%s', [%s], {priority: 'event'})", input_id, r$id),
      tags$span(class = "comp-name", style = sprintf("padding-left: %.1frem", r$depth * 0.8),
                r$phenotype),
      tags$span(class = "comp-bar",
                tags$span(class = "comp-total", style = sprintf("width: %.1f%%", r$total)),
                tags$span(class = "comp-own", style = sprintf("width: %.1f%%", r$remaining))),
      tags$span(class = "comp-value", sprintf("%.1f%%", r$total)),
      tags$span(class = "comp-parent",
                if (is.na(r$of_parent)) "" else sprintf("%.0f%%", r$of_parent))
    )
  })
  tags$div(
    class = "composition",
    tags$div(class = "comp-row comp-head",
             tags$span(class = "comp-name", "Phenotype"),
             tags$span(class = "comp-bar", ""),
             tags$span(class = "comp-value", title = "Share of all cells, incl. subpopulations", "Cells"),
             tags$span(class = "comp-parent", title = "Share of the parent phenotype", "Parent")),
    rows,
    tags$div(class = "comp-legend",
             tags$span(class = "swatch total"), "incl. subpopulations",
             tags$span(class = "swatch own ms-2"), "remaining")
  )
}
