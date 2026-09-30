#' Start the CyCadas app
#'
#' @export
cycadas <- function() {
  shinyApp(ui = app_ui(), server = app_server)
}
