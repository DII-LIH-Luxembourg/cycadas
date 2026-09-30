# Application server ----------------------------------------------------------
# Creates the per-session state and starts one module per tab. Modules only
# communicate through `state`.

app_server <- function(input, output, session) {

  state <- new_app_state()

  mod_data_server("settings", state)
  mod_explore_server("umap_reactive", state)
  mod_thresholds_server("thresholds", state)
  mod_annotation_server("treeannotation", state)
  mod_da_server("DA_tab", state)
  mod_da_tree_server("DA_tree", state)
}
