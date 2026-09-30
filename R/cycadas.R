#' Start the CyCadas app
#'
#' @export
cycadas <- function() {
  # images referenced by HOME.md and the navbar
  addResourcePath("www", app_file("www"))
  shinyApp(ui = app_ui(), server = app_server)
}
