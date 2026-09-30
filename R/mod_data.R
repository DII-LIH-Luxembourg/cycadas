# Workspace tab ---------------------------------------------------------------
# Lays out the import and workspace sub-modules.

mod_data_ui <- function(id) {
  ns <- NS(id)
  tags$div(
    class = "container-xl py-3",
    layout_columns(
      col_widths = c(7, 5),
      tagList(
        navset_card_underline(
          title = "Import cluster data",
          nav_panel("GigaSOM / FlowSOM", mod_import_gigasom_ui(ns("gigasom"))),
          nav_panel("CATALYST", mod_import_catalyst_ui(ns("catalyst"))),
          nav_panel("RemoteSOM", mod_import_remotesom_ui(ns("remotesom"))),
          nav_panel("Demo data", mod_workspace_demo_ui(ns("workspace")))
        ),
        card(card_header("Optional files"),
             mod_import_optional_ui(ns("optional")))
      ),
      mod_workspace_ui(ns("workspace"))
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
