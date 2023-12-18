ui <- dashboardPage(
  dashboardHeader(title = "Cluster Annotation"),

  # Tab menu layout ---------------------------------------------------------
  sidebar <- dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("Load", tabName = "settings"),
      menuItem("Thresholds", tabName = "thresholds"),
      menuItem("Tree-Annotation", tabName = "treeannotation"),
      menuItem("UMAP interactive", tabName = "umap_reactive"),
      menuItem("UMAP Marker expression", tabName = "UMAP_Marker_expression"),
      menuItem("Differential Abundance", tabName = "DA_tab"),
      menuItem("DA interactive Tree", tabName = "DA_tree"),
      menuItem("Help", tabName = "help")
    )
  ),
  dashboardBody(
    #Custom CSS
    tags$head(tags$style("#thTableBox{height:400px !important;}")),
    #Custom CSS
    tags$head(tags$style("#mytreebox{height:700px !important;}")),
    # Load/Settings Tab ----------------------------------------------------------------
    tabItems(
      tabItem(tabName = "settings",Settings_UI(id="Settings")),
      tabItem(tabName = "thresholds",Threshold_UI(id="threshold")),
      tabItem(tabName = "treeannotation",
              fluidRow(
                 column(width = 4,
                        selectizeInput(inputId = ("parentPicker"),label = "Select Parent:",choices=NULL),
                        TreeAnnotation_UI(id="TreeAnnotation"),
                        box(width = NULL,title = "Delete Node",actionButton(("deleteNodeBtn"), "Delete Node"))),
                 column(width = 8,
                        box(width = NULL,title = "Annotation Tree",visNetwork_UI(id="visNetwork")),
                        box(width = NULL,title = "Heatmap",Heatmap_UI(id="Heatmap")),
                        box(width = NULL,title = "Umap",Umap_UI(id="Umap")))
                 )
              ),
      # umap Reactive  Tab ----------------------------------------------------
      tabItem(tabName = "umap_reactive",umap_interactive_UI(id="umap_interactive")),
      # Umap Marker expression Tab --------------------------------------------
      tabItem(tabName = "UMAP_Marker_expression",UMAP_ME_UI(id="UMAP_ME")),
      # Differential Abundance Tab --------------------------------------------
      tabItem(tabName = "DA_tab",Differential_Abundance_UI(id="Differential_Abundance")),
      tabItem(tabName = "DA_tree",
              fluidRow(column(width = 8,
                              box(id="mytreebox",
                                  width = NULL,
                                  title = "Interactive DA Tree",
                                  visNetworkOutput("interactiveTree", width = "100%", height = "700px")
                              )
              ),DA_interactive_UI(id="DA_interactive")))
              
              
    ) # tabItems
  ) # dashboardBody
)

