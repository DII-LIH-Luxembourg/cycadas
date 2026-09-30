# Home tab --------------------------------------------------------------------
# Renders inst/app/HOME.md.

mod_home_ui <- function(id) {
  fluidRow(align = "left",
           box(width = NULL,
               includeMarkdown(app_file("HOME.md")),
               tags$hr()
           )
  )
}

# Path of a file shipped in inst/app
app_file <- function(...) {
  system.file("app", ..., package = "cycadas", mustWork = TRUE)
}
