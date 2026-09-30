# Home tab --------------------------------------------------------------------
# Renders inst/app/HOME.md.

mod_home_ui <- function(id) {
  tags$div(class = "container-lg py-3",
           card(card_body(class = "home-content", includeMarkdown(app_file("HOME.md")))))
}

# Path of a file shipped in inst/app
app_file <- function(...) {
  system.file("app", ..., package = "cycadas", mustWork = TRUE)
}
