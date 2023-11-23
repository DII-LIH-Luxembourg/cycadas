library(ggplot2)
library(shiny)


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


chartServer <- function(id, x, y, title) {
  moduleServer(
    id = id,
    module = function(input, output, session) {

    }
  )
}