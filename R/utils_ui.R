# UI helpers ------------------------------------------------------------------

# truthy fileInput?
has_file <- function(x) !is.null(x) && !is.null(x$datapath) && nzchar(x$datapath)

# pretty status items ----
status_item <- function(ok, label) {
  col <- if (ok) "#28a745" else "#dc3545"  # green/red
  icon <- if (ok) "✔" else "✖"
  tags$div(
    style = "margin:4px 0;",
    tags$span(style = sprintf(
      "display:inline-block;width:10px;height:10px;border-radius:50%%;background:%s;margin-right:6px;", col)),
    tags$span(icon, style = sprintf("color:%s;margin-right:6px;", col)),
    tags$span(label)
  )
}

# Status lights for a set of required files plus a summary line
file_status_ui <- function(present, labels, ready_msg, missing_msg) {
  tagList(
    Map(status_item, present, labels),
    tags$div(style = "margin-top:8px;",
             if (all(present))
               tags$span(ready_msg, style = "color:#28a745;")
             else
               tags$span(missing_msg, style = "color:#dc3545;")
    )
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
picker_options <- function() {
  list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
}

# Collapsible dashboard box with the settings used on the Workspace tab
settings_box <- function(title, ..., status = "info", collapsed = FALSE) {
  box(title = title, collapsible = TRUE, solidHeader = TRUE, status = status,
      width = NULL, collapsed = collapsed, ...)
}
