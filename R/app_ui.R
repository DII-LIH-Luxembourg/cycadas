# Application UI --------------------------------------------------------------
# Every tab is a module; the tab name doubles as the module id.

app_ui <- function() {
  dashboardPage(
    dashboardHeader(title = "CyCadas"),
    dashboardSidebar(
      sidebarMenu(id = "tabs",
                  menuItem("Home", tabName = "home"),
                  menuItem("Workspace", tabName = "settings"),
                  menuItem("Explore Data", tabName = "umap_reactive"),
                  menuItem("Thresholds", tabName = "thresholds"),
                  menuItem("Tree-Annotation", tabName = "treeannotation"),
                  menuItem("Differential Abundance", tabName = "DA_tab"),
                  menuItem("DA interactive Tree", tabName = "DA_tree")
      )
    ),
    dashboardBody(
      useShinyjs(),
      tabItems(
        tabItem(tabName = "home", mod_home_ui("home")),
        tabItem(tabName = "settings", mod_data_ui("settings")),
        tabItem(tabName = "umap_reactive", mod_explore_ui("umap_reactive")),
        tabItem(tabName = "thresholds", mod_thresholds_ui("thresholds")),
        tabItem(tabName = "treeannotation", mod_annotation_ui("treeannotation")),
        tabItem(tabName = "DA_tab", mod_da_ui("DA_tab")),
        tabItem(tabName = "DA_tree", mod_da_tree_ui("DA_tree"))
      )
    )
  )
}
