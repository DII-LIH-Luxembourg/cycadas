# Workspace tab ---------------------------------------------------------------
# Lays out the import and workspace sub-modules.

mod_data_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(width = 4,
             mod_import_gigasom_ui(ns("gigasom")),
             mod_import_catalyst_ui(ns("catalyst")),
             mod_import_remotesom_ui(ns("remotesom"))
      ),
      column(width = 4,
             mod_import_optional_ui(ns("optional"))
      ),
      column(width = 4,
             mod_workspace_ui(ns("workspace"))
      )
    )
  )
}

mod_data_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    mod_import_gigasom_server("gigasom", state)
    mod_import_catalyst_server("catalyst", state)
    mod_import_remotesom_server("remotesom", state)
    mod_import_optional_server("optional", state)
    mod_workspace_server("workspace", state)
  })
}
