
ui <- dashboardPage(
  dashboardHeader(title = "CyCadas"),
  # Tab menu layout ---------------------------------------------------------
  sidebar <- dashboardSidebar(
    sidebarMenu(id = "tabs",
                menuItem("Home", tabName = "home"),
                menuItem("Load", tabName = "settings"),
                menuItem("Thresholds", tabName = "thresholds"),
                menuItem("Tree-Annotation", tabName = "treeannotation"),
                menuItem("UMAP interactive", tabName = "umap_reactive"),
                menuItem("UMAP Marker expression", tabName = "UMAP_Marker_expression"),
                menuItem("Differential Abundance", tabName = "DA_tab"),
                menuItem("DA interactive Tree", tabName = "DA_tree")
    )
  ),
  dashboardBody(
    #Custom CSS
    tags$head(tags$style("#thTableBox{height:400px !important;}")),
    #Custom CSS
    tags$head(tags$style("#mytreebox{height:700px !important;}")),
    # Load/Settings Tab ----------------------------------------------------------------
    tabItems(
      tabItem(tabName = "home",
              fluidRow(align = "center",
                tags$img(src="./images/logo4_small2.png"),
                  ),
              fluidRow(align = "left",
                       box(width = NULL,
                              HTML(markdown::markdownToHTML(knit("HOME.md", quiet=T),fragment.only = T)),
                              tags$hr()
                       )
                       )

              ),
      tabItem(tabName = "settings",
              fluidPage(
                fluidRow(
                  column(width = 4,
                         box(title = "Requried",collapsible = TRUE,solidHeader = TRUE,status = "info",width = NULL,collapsed = F,
                             fileInput("fMarkerExpr","Upload Marker Expressions",
                                       placeholder = "Choose CSV File",
                                       multiple = FALSE,
                                       accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                             fileInput("cluster_freq","Upload Cluster Frequencies",
                                       placeholder = "Choose CSV File",
                                       multiple = FALSE,
                                       accept = c("text/csv","text/comma-separated-values,text/plain",".csv"))
                         )),
                  column(width = 4,
                         box(title = "Optional - Marker-thresholds",collapsible = TRUE,solidHeader = TRUE,status = "warning",width = NULL,collapsed = T,
                             fileInput("fTH","Choose CSV File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                             tags$hr(),
                             checkboxInput("header", "Header", TRUE),
                             radioButtons("sep","Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ",",inline = T)),
                         box(title = "Optional - Annotation Tree",collapsible = TRUE,solidHeader = TRUE,status = "warning",width = NULL,collapsed = T,
                             fileInput("fNodes","Choose Nodes File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                             fileInput("fEdges","Choose Edges File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                             # fileInput("fAnno","Choose Annotation File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                             actionButton("btnImportTree", "Import")),
                         box(title = "Optional - Metadata",collapsible = TRUE,solidHeader = TRUE,status = "warning",width = NULL,collapsed = T,
                             fileInput("metadata","Choose CSV File",multiple = F,accept = c("text/csv","text/comma-separated-values,text/plain",".csv"))),
                         box(title = "Optional - Count Table",collapsible = TRUE,solidHeader = TRUE,status = "warning",width = NULL,collapsed = T,
                             fileInput("counts_table","Choose CSV File",multiple = F,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")))
                  ),
                  column(width = 4,
                         box(title = "Load Demo Data",collapsible = TRUE,solidHeader = TRUE,status = "success",width = NULL,collapsed = F,
                             splitLayout(cellWidths = c("50%", "50%"),
                                         actionButton("btnLoadDemoData", "Unannotated"),
                                         actionButton("btnLoadAnnoData", "Annotated")))
                  )
                )
              )
              ),

      # Tree-Annotation Tab ----------------------------------------------------
      tabItem(tabName = "treeannotation",
              fluidRow(column(
                width = 4,
                box(
                  width = NULL,
                  title = "Create Node",
                  textOutput("progressBox"),
                  textOutput("progressBox2"),
                  textInput("newNode", "Set Name..."),
                  pickerInput(
                    inputId = "parentPicker",
                    label = "Select Parent:",
                    choices = NULL,
                    options = list(
                      `actions-box` = TRUE,
                      size = 10,
                      `selected-text-format` = "count > 3"
                    ),
                    multiple = F
                  )
                ),
                box(
                  width = NULL,
                  column(
                    width = 4,
                    checkboxGroupButtons(
                      inputId = "treePickerPos",
                      label = "Positive:",
                      choices = c("A"),
                      direction = "vertical"
                    )
                  ),
                  column(
                    width = 4,
                    checkboxGroupButtons(
                      inputId = "treePickerNeg",
                      label = "Negative:",
                      choices = c("A"),
                      direction = "vertical"
                    )
                  )

                ),

                box(
                  width = NULL,
                  # switchInput(
                  #   inputId = "treeSwitch",
                  #   label = "Change Layout",
                  #   value = FALSE,
                  #   onLabel = "Hierarchical",
                  #   offLabel = "Loose"
                  # ),
                  title = "Create New Node",
                  actionButton("createNodeBtn", "Create Node")
              ),
              # box(
              #   width = NULL,
              #   title = "Update Node",
              #   actionButton("updateNodeBtn", "Update Node"),
              #
              # ),
              box(
                width = NULL,
                title = "Delete Node",
                actionButton("deleteNodeBtn", "Delete Node")
              ),
              # box(
              #   width = NULL,
              #   title = "Update Node",
              #   pickerInput(
              #     inputId = "updateNodePicker",
              #     label = "Select Node:",
              #     choices = NULL,
              #     options = list(
              #       `actions-box` = TRUE,
              #       size = 10,
              #       `selected-text-format` = "count > 3"
              #     ),
              #     multiple = F
              #   ),
              #   textInput("renameNode", "Rename..."),
              #   # select markers
              #   pickerInput(
              #     inputId = "updatePickerPos",
              #     label = "Update Positive Markers",
              #     choices = NULL,
              #     options = list(
              #       `actions-box` = TRUE,
              #       size = 10,
              #       `selected-text-format` = "count > 3"
              #     ),
              #     multiple = TRUE
              #   ),
              #   pickerInput(
              #     inputId = "updatePickerNeg",
              #     label = "Update Negative Markers",
              #     choices = NULL,
              #     options = list(
              #       `actions-box` = TRUE,
              #       size = 10,
              #       `selected-text-format` = "count > 3"
              #     ),
              #     multiple = TRUE
              #   ),
              #   actionButton("updateNodeBtn", "Update Node"),
              #   actionButton("deleteNodeBtn", "Delete Node")
              # ),
              box(
                width = NULL,
                title = "Export Annotation",
                # actionButton("exportAnnotationBtn", "Export Annotation")
                downloadButton("exportAnnotationBtn", "Export Annotation")
              )
              # box(
              #   width = NULL,
              #   title = "Export Tree as Image",
              #   actionButton("exportTreeGraphics", "Export Tree Image")
              # )
              ), # end column
              column(
                width = 8,
                box(
                  width = NULL,
                  title = "Annotation Tree",
                  visNetworkOutput("mynetworkid")
                ),
                box(
                  width = NULL,
                  title = "Heatmap",
                  plotOutput("hm_tree")
                ),
                box(width = NULL,
                    plotOutput("umap_tree"))
              ) # end column
              ) # end row

              ),
      # Thresholds Tab --------------------------------------------------------
      tabItem(tabName = "thresholds",
              fluidRow(column(
                width = 6,
                box(
                  width = NULL,
                  title = "Marker Expression:",
                  plotOutput("plot", click = "plot_click")
                ),
                box(width = NULL,
                    title = "Histogram:",
                    plotOutput("plot2"))
              ),
              column(
                width = 6,
                # box(
                #   width = NULL,
                #   id = "transformationBox",
                #   title = "Data transformation",
                #   radioButtons(
                #     "radio",
                #     label = NULL,
                #     choices = list("0 to 1" = "1"),
                #     selected = "1"
                #   )
                # ),
                box(width = NULL,
                    selectInput(
                      "th_method",
                      "Threshold Method:",
                      c("K-means" = "km",
                        "GMM mid" = "gmm_mid",
                        "GMM-high" = "gmm_high",
                        "GMM-low" = "gmm_low"),
                      selected = "km",
                      multiple = FALSE,
                      selectize = TRUE,
                      width = NULL,
                      size = NULL
                    )),
                box(width = NULL,
                    DTOutput('table'))
              ))),
      # umap Reactive  Tab ----------------------------------------------------
      tabItem(tabName = "umap_reactive",
              fluidRow(column(width = 6,
                              box(
                                width = NULL,
                                plotOutput("umap2", brush = "umap2_brush")
                              )),
                       column(width = 6,
                              box(
                                width = NULL,
                                plotOutput("hm2")
                              ))),
              fluidRow(column(width = 12,
                              box(
                                width = NULL,
                                DTOutput("umap_data")
                              )))),
      # Umap Marker expression Tab --------------------------------------------
      tabItem(tabName = "UMAP_Marker_expression",
              fluidRow(column(width = 10,
                              box(
                                width = NULL,
                                plotOutput("umap3")
                              )),
                       column(width = 2,
                              box(
                                width = NULL,
                                selectInput("markerSelect", "Select:", choices = NULL)
                              )))),
      # Differential Abundance Tab --------------------------------------------
      tabItem(tabName = "DA_tab",
              fluidRow(column(width = 6, # left column
                              
                              box(width = NULL,
                                  fluidRow(
                                    column(width = 5,
                                           box(
                                             width = NULL,
                                             title = "Metadata preview",
                                             tableOutput("md_table"))
                                           ),
                                    column(width = 5,
                                           box(
                                             width = NULL,
                                             title = "Counts Table preview",
                                             tableOutput("counts_table"))
                                    )  
                                  ),
                              ),

                              box(
                                width = NULL,
                                fluidRow(column(width = 6,
                                                 box(
                                                   width = NULL,
                                                   selectInput("correction_method", "Select:", choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                                                                           "fdr", "none"))
                                                 )
                                                ),
                                         column(width = 6,
                                                 box(
                                                   width = NULL,
                                                   title = "Do Analysis",
                                                   actionButton("doDA", "Calculate")
                                                 )
                                         )
                                
                                )
                              ),
                              box(
                                width = NULL,
                                fluidRow(column(width = 6,
                                                box(
                                                  width = NULL,
                                                  title = "Export DA Result",
                                                  downloadButton("exportDA", "Download")
                                                )
                                        ),
                                        
                                        column(width = 6,
                                               box(
                                                 width = NULL,
                                                 title = "Export Proportion Table",
                                                 downloadButton("exportProp", "Download")
                                               )
                                        ),
                                
                                )  
                              ),
                              
                              ),
                       column(width = 6, # right column
                              box(
                                width = NULL,
                                title = "DA Result",
                                tableOutput("DA_result_table")
                              )
                              )
                       )),
      tabItem(tabName = "DA_tree",
              fluidRow(column(width = 8,
                              box(id="mytreebox",
                                width = NULL,
                                title = "Interactive DA Tree",
                                visNetworkOutput("interactiveTree", width = "100%", height = "700px")
                              )
                              ),
                       column(width = 4,
                              box(
                                width = NULL,
                                tableOutput("DA_interactive_table")
                              ),
                              box(
                                width = NULL,
                                plotOutput("boxplot")
                              )
                              )

                       ))# tabItem
    ) # tabItems
  ) # dashboardBody
)
