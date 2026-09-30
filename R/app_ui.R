# Application UI --------------------------------------------------------------
# Every tab is a module; the tab value doubles as the module id.

cycadas_theme <- function() {
  bs_theme(version = 5, preset = "shiny", primary = "#2c7be5")
}

app_ui <- function() {
  page_navbar(
    id = "tabs",
    title = tags$img(src = "www/logo_navbar.png", alt = "CyCadas"),
    window_title = "CyCadas",
    theme = cycadas_theme(),
    fillable = FALSE,
    header = tagList(useShinyjs(), cycadas_dependency()),
    nav_panel("Home", value = "home", mod_home_ui("home")),
    nav_panel("Workspace", value = "settings", mod_data_ui("settings")),
    nav_panel("Explore", value = "umap_reactive", mod_explore_ui("umap_reactive")),
    nav_panel("Thresholds", value = "thresholds", mod_thresholds_ui("thresholds")),
    nav_panel("Annotation", value = "treeannotation", mod_annotation_ui("treeannotation")),
    nav_menu(
      "Differential abundance",
      nav_panel("Results table", value = "DA_tab", mod_da_ui("DA_tab")),
      nav_panel("Interactive tree", value = "DA_tree", mod_da_tree_ui("DA_tree"))
    ),
    nav_spacer(),
    nav_item(uiOutput("dataset_badge"))
  )
}
