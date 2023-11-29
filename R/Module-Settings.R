# library(ggplot2)
# library(shiny)


Settings_UI <- function(id) {
  ns <- NS(id)
  fluidRow(column(width = 12,
    tabBox(width = NULL,title = "Data upload",id = "tab-required",height = "450px",
           tabPanel("Required",
                    "Required",
                    
        box(title = "Upload Marker-expression",collapsible = TRUE,solidHeader = TRUE,status = "info",width = NULL,collapsed = F,
          fileInput("fMarkerExpr","Choose CSV File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
          tags$hr(),
          checkboxInput("header", "Header", TRUE),
          radioButtons("sep","Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ",",inline = T)),
        
        box(title = "Upload cluster_freq",collapsible = TRUE,solidHeader = TRUE,status = "info",width = NULL,collapsed = F,
            fileInput("cluster_freq","Choose CSV File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
          tags$hr(),
          checkboxInput("header", "Header", TRUE),
          radioButtons("sep","Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ",",inline = T))
        ),
        
      tabPanel("Optional",
               "Optional",
               
        box(title = "Upload Marker-thresholds",collapsible = TRUE,solidHeader = TRUE,status = "info",width = NULL,collapsed = F,
          fileInput("fTH","Choose CSV File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
          tags$hr(),
          checkboxInput("header", "Header", TRUE),
          radioButtons("sep","Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ",",inline = T)),
        
        box(title = "Upload Annotation Tree",collapsible = TRUE,solidHeader = TRUE,status = "info",width = NULL,collapsed = F,
          fileInput("fNodes","Choose Nodes File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
          fileInput("fEdges","Choose Edges File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
          fileInput("fAnno","Choose Annotation File",multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
          actionButton("btnImportTree", "Import")),
        
        box(title = "Upload Metadata",collapsible = TRUE,solidHeader = TRUE,status = "info",width = NULL,collapsed = F,
            fileInput("metadata","Choose CSV File",multiple = F,accept = c("text/csv","text/comma-separated-values,text/plain",".csv"))),
        box(title = "Upload Count Table",collapsible = TRUE,solidHeader = TRUE,status = "info",width = NULL,collapsed = F,
            fileInput("counts_table","Choose CSV File",multiple = F,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")))),
      
      tabPanel("Demo Data",
               box(title = "Load Cluster Expression Demo Data",width = NULL,actionButton("btnLoadDemoData", "Load")),
               box(title = "Load Annotated Demo Data",width = NULL,actionButton("btnLoadAnnoData", "Load")))
    )
  ))
}

Settings_Server1 <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {

# browser()
# Create a Progress object
progress <- shiny::Progress$new()
# Make sure it closes when we exit this reactive, even if there's an error
on.exit(progress$close())

progress$set(message = "loading Data Cluster Expression Demo Data...", value = 0)

# browser()

## Load median expression and cell frequencies
df <- read.csv("data/demo_data/median_expr_1600.csv")
df_global <<- df

reactVals$graph <- initTree()

annotationlist <<- list("Unassigned")

cell_freq <- read.csv("data/demo_data/cluster_freq_1600.csv")

progress$set(message = "loading Data Cluster Expression Demo Data...", value = 0.2)

labels_row <-
  paste0(rownames(df), " (", cell_freq$clustering_prop , "%)")

marker_names <- rownames(df)

set.seed(1234)

my_umap <- umap(df)
dr_umap <<- data.frame(
  u1 = my_umap$layout[, 1],
  u2 = my_umap$layout[, 2],
  my_umap$data,
  cluster_number = 1:length(my_umap$layout[, 1]),
  check.names = FALSE
)

progress$set(message = "loading Data Cluster Expression Demo Data...", value = 0.4)

df01 <<- df %>% normalize01()
myDF <<- df %>% normalize01()

df01Tree <<- df %>% normalize01()
allMarkers <<- colnames(df)
df01Tree$cell <<- "Unassigned"

selectedMarkers <<- colnames(df)
posPickerList <<- colnames(df)


progress$set(message = "loading Data Cluster Expression Demo Data...", value = 0.6)

## load only at start to fill the picker list
updatePickerInput(
  session,
  inputId = "myPickerPos",
  label = "Select Positive Markers",
  choices = colnames(df01),
  # choices = NULL,
  options = list(
    `actions-box` = TRUE,
    size = 10,
    `selected-text-format` = "count > 3"
  )
)
updatePickerInput(
  session,
  inputId = "myPickerNeg",
  label = "Select Negative Markers",
  choices = colnames(df01),
  options = list(
    `actions-box` = TRUE,
    size = 10,
    `selected-text-format` = "count > 3"
  )
)
updateSelectInput(session, "markerSelect", "Select:", colnames(df01))

updateCheckboxGroupButtons(
  session,
  inputId = "treePickerPos",
  choices = colnames(df01),
  selected = NULL
)
updateCheckboxGroupButtons(
  session,
  inputId = "treePickerNeg",
  choices = colnames(df01),
  selected = NULL
)

annotaionDF <<- data.frame("cell" = "unassigned",
                           clusterSize = cell_freq$clustering_prop)
at <<- reactiveValues(data = annotaionDF, dr_umap = dr_umap)

reactVals$th <- kmeansTH(df01)


})}
Settings_Server2 <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      ## -------------------------------------------------------------------
      # Upload All Annotation Demo Data ----

        # browser()
        
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "loading Data Annotation Demo Data...", value = 0)
        
        # browser()
        
        ## Load median expression and cell frequencies
        df <- read.csv("data/demo_data/median_expr_1600.csv")
        # !----------- TEST bimodal check !-----------------
        #
        df$testCol <- rnorm(nrow(df), 5.0, 1.0)
        
        df_global <<- df
        
        reactVals$graph <- initTree()
        
        annotationlist <<- list("Unassigned")
        
        progress$set(message = "loading Data Annotation Demo Data...", value = 0.1)
        
        cell_freq <<- read.csv("data/demo_data/cluster_freq_1600.csv")
        
        progress$set(message = "loading Data Annotation Demo Data...", value = 0.2)
        
        labels_row <-
          paste0(rownames(df), " (", cell_freq$clustering_prop , "%)")
        
        marker_names <- rownames(df)
        
        set.seed(1234)
        
        my_umap <- umap(df)
        dr_umap <<- data.frame(
          u1 = my_umap$layout[, 1],
          u2 = my_umap$layout[, 2],
          my_umap$data,
          cluster_number = 1:length(my_umap$layout[, 1]),
          check.names = FALSE
        )
        progress$set(message = "loading Data Annotation Demo Data...", value = 0.4)
        
        df01 <<- df %>% normalize01()
        myDF <<- df %>% normalize01()
        
        df01Tree <<- df %>% normalize01()
        allMarkers <<- colnames(df)
        df01Tree$cell <<- "Unassigned"
        
        selectedMarkers <<- colnames(df)
        posPickerList <<- colnames(df)
        
        progress$set(message = "loading Data Annotation Demo Data...", value = 0.6)
        ## load only at start to fill the picker list
        updatePickerInput(session,inputId = "myPickerPos",label = "Select Positive Markers",choices = colnames(df01),
                          options = list(`actions-box` = TRUE,size = 10,`selected-text-format` = "count > 3"))
        updatePickerInput(session,inputId = "myPickerNeg",label = "Select Negative Markers",choices = colnames(df01),
                          options = list(`actions-box` = TRUE,size = 10,`selected-text-format` = "count > 3"))
        updateSelectInput(session, "markerSelect", "Select:", colnames(df01))
        updateCheckboxGroupButtons(session,inputId = "treePickerPos",choices = colnames(df01),selected = NULL)
        updateCheckboxGroupButtons(session,inputId = "treePickerNeg",choices = colnames(df01),selected = NULL)
        
        annotaionDF <<- data.frame("cell" = "unassigned",clusterSize = cell_freq$clustering_prop)
        at <<- reactiveValues(data = annotaionDF, dr_umap = dr_umap)
        
        
        # reactVals$th <- kmeansTH(df01)
        
        ## Load demo thresholds
        th <- read.csv("data/demo_data/MarkerThresholds.csv")
        
        th$X <- NULL
        th$color <- "blue"
        
        th[th$bi_mod < 0.555, "color"] <- "red"
        
        rownames(th) <- th$cell
        reactVals$th<- th
        
        ## Load metadata
        md <<- read.csv("data/demo_data/metadata.csv")
        md$X <- NULL
        reactVals$md<- md
        
        ## Load counts table
        ct <<- read.csv("data/demo_data/cluster_counts_1600.csv")
        ct$X <- NULL
        reactVals$counts_table <- ct
        
        ## Load annotaiton Tree
        df_nodes <- read.csv("data/demo_data/nodesTable_data.csv")
        df_edges <- read.csv("data/demo_data/edgesTable_data.csv")
        df_anno <- read.csv("data/demo_data/annTable_data.csv")
        
        df_nodes$pm[is.na(df_nodes$pm)] <- ""
        df_nodes$pm <- strsplit(df_nodes$pm, "\\|")
        
        df_nodes$nm[is.na(df_nodes$nm)] <- ""
        df_nodes$nm <- strsplit(df_nodes$nm, "\\|")
        
        # browser()
        reactVals$graph$nodes <- df_nodes
        reactVals$graph$edges <- df_edges
        
        df01Tree <<- df_anno
        
        reactVals$hm <- df01Tree
        
        annotationlist <<- as.list(df_nodes$label)

    }
  )
}
Settings_Server3 <- function(id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {

      # Load Annotation Tree ----
        # browser()
        
        req(input$fNodes)
        req(input$fEdges)
        req(input$fAnno)
      
        req(input$fMarkerExpr)
        req(input$cluster_freq)
        
        df_nodes <- read.csv(input$fNodes$datapath)
        df_edges <- read.csv(input$fEdges$datapath)
        df_anno <- read.csv(input$fAnno$datapath)
        
        df_nodes$pm[is.na(df_nodes$pm)] <- ""
        df_nodes$pm <- strsplit(df_nodes$pm, "\\|")
        
        df_nodes$nm[is.na(df_nodes$nm)] <- ""
        df_nodes$nm <- strsplit(df_nodes$nm, "\\|")
        
        reactVals$graph$nodes <- df_nodes
        reactVals$graph$edges <- df_edges
        
        df01Tree <<- df_anno
        
        annotationlist <<- as.list(df_nodes$label)
        
    })}