#' Start the CyCadas app
#'
#' @param dir Optional folder with `median_expr.csv` and `cluster_freq.csv`
#'   (see fct_folder.R). The app opens it at start and gets a "Save" button
#'   that writes the workspace and the cluster labels into it; a saved
#'   workspace there is reopened.
#' @export
cycadas <- function(dir = NULL) {
  # images referenced by HOME.md and the navbar
  addResourcePath("www", app_file("www"))
  if (!is.null(dir)) dir <- normalizePath(dir, mustWork = TRUE)
  shinyApp(ui = app_ui(folder = !is.null(dir)),
           server = function(input, output, session) app_server(input, output, session, dir = dir))
}
