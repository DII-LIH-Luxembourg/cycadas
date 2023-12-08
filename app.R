
# List of packages you want to check and install if needed
packages_to_install <- c("shiny", "DT", "ggplot2", "matrixStats", "tidyverse", "stats",
                         "pheatmap", "Ckmeans.1d.dp", "umap", "RColorBrewer", "shinydashboard",
                         "shinyWidgets", "visNetwork", "glue", "purrr", "reshape2", "mousetrap")

# Check if each package is already installed, and install if not
for (package in packages_to_install) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }

}





source("modules/utils.R")
source("modules/Module-Settings.R")
source("modules/Module-Threshold.R")
# source("modules/Module-parentPicker.R")
source("modules/Module-TreeAnnotation.R")
source("modules/Module-Heatmap.R")
source("modules/ui.R")
# cycadas <- function() {

  # browser()



  # constructs a string of positive or negative markers
  ph_name <<- ""

  # Server function ----
  server = function(input, output, session) {
    
    reactVals <- reactiveValues(th = NULL,
                                myTH = NULL,
                                md = NULL,
                                counts_table = NULL,
                                DA_result_table = NULL,
                                DA_interactive_table = NULL,
                                graph = NULL,
                                hm = NULL,
                                merged_prop_table = NULL,
                                ParentValue=NULL)

    
    # Save the DA Table ----
    output$exportDA <- downloadHandler(
      filename = function() {
        paste("DA_Table_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(reactVals$DA_result_table, file)
      }
    )
    
    get_merged_prop_table <- function() {
      
      req(reactVals$counts_table)
      
      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df01Tree$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)
      
      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL
      
      props_table <- t(t(countsTable) / colSums(countsTable)) * 100
      
      return(props_table)
      
    }
    
    # Save the Merged Proportion Table ----
    output$exportProp <- downloadHandler(
      
      filename = function() {
        paste("Merged_Proportions_Table_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(get_merged_prop_table(), file)
      }
    )
    

    # Load the marker threshold file ---------------------------------------
    observeEvent(input$fTH,{

      # browser()

      th <<- read.csv(input$fTH$datapath)
      th$X <- NULL
      th$color <- "blue"

      th[th$bi_mod < 0.555, "color"] <- "red"

      rownames(th) <- th$cell

      reactVals$th<- th

    })

    # Load the metadata ----------------------------------------------------
    observeEvent(input$metadata,{
      md <<- read.csv(input$metadata$datapath)
      md$X <- NULL
      reactVals$md<- md
    })

    # Load the counts talbe -------------------------------------------------
    observeEvent(input$counts_table,{
      ct <<- read.csv(input$counts_table$datapath)
      ct$X <- NULL
      reactVals$counts_table <- ct
    })

    ## create Phenotype name from picker e.g. CD4+CD8+CD19 --------------------
    observeEvent(input$myPickerPos, {
      output$r1 <- renderText(setPhenotypeName(input$myPickerPos, "pos", ph_name))
    })
    observeEvent(input$myPickerNeg, {
      output$r2 <- renderText(setPhenotypeName(input$myPickerNeg, "neg", ph_name))
    })

    ## Annotation Table -------------------------------------------------------
    observeEvent(c(input$myPickerPos, input$myPickerNeg, reactVals$th),{
      req(input$fMarkerExpr)
      req(input$cluster_freq)

      tmp=filterHM(df01,input$myPickerPos, input$myPickerNeg, reactVals$th)

      # update myDF to the filtered HM
      myDF <<- tmp

      mydf=annotaionDF[rownames(tmp),]
      rownames(mydf)=rownames(tmp)

      # update heatmap plot ----
      output$hm <- renderPlot({
        if (dim(tmp)[1] < 2) {
          pheatmap(tmp, cluster_cols = F, cluster_rows = F)
        } else {
          message(dim(tmp)[1])
          pheatmap(tmp, cluster_cols = F)
        }
      })

      # update umap plot ----
      myColor <- filterColor(df_global,tmp)
      Temp1<<-myColor
      Temp2<<-dr_umap

      output$umap <-
        renderPlot(
          ggplot(dr_umap, aes(
            x = u1, y = u2, color = myColor
          )) +
            geom_point(size = 1.0) +
            theme_bw() +
            theme(legend.text = element_text(size =
                                               8)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
        )
      # update table ----
      output$tableAnnotation <- DT::renderDT(
        at$data,
        extensions = c('Buttons', 'Scroller'),
        options = list(
          dom = 'Bfrtip',
          buttons = c('csv'),
          paging = FALSE
        )
      )

      # update valueboxes ----
      mysum <- sum(mydf$clusterSize)
      output$progressBox <- renderValueBox({
        valueBox(paste0(mysum, "%"), "Selected", icon = icon("list"),color = "purple")
      })
      output$progressBox2 <- renderValueBox({
        valueBox(dim(mydf)[1], "Cluster", icon = icon("list"),color = "purple")
      })
    })

    # Set Phenotype names -----------------------------------------------------
    observeEvent(input$btnSetType, {

      # if cluster has already been assigned, add the new name after
      for (n in rownames(myDF)) {
        if (at$data[n, "cell"] == "unassigned") {
          at$data[n, "cell"] <- input$phenotype
        } else(
          at$data[n, "cell"] <- paste0(at$data[n, "cell"], "_", input$phenotype)
        )

      }
    })

    # Server - UMAP interactive Tab -------------------------------------------
    output$umap2 <- renderPlot(
      ggplot(dr_umap, aes(x = u1, y = u2)) +
        geom_point(size = 1.0) +
        theme_bw() +
        theme(legend.text=element_text(size=8)) +
        guides(color = guide_legend(override.aes = list(size = 4)))
    )

    output$hm2 <- renderPlot({
      my_new_hm <- round(brushedPoints(dr_umap, input$umap2_brush), 2)

      if (dim(my_new_hm)[1] > 0) {

        pheatmap(my_new_hm %>% select(-c("cluster_number", "u1", "u2")), cluster_cols = F)
      } else {
        ggplot() + theme_void() + ggtitle("Select area on Umap to plot Heatmap")
      }
    })

    output$umap_data <-
      DT::renderDT(server = FALSE, {
        DT::datatable(
          round(brushedPoints(dr_umap, input$umap2_brush), 2),
          filter = 'top',
          extensions = 'Buttons',
          options = list(
            scrollY = 600,
            scrollX = TRUE,
            dom = '<"float-left"l><"float-right"f>rt<"row"<"col-sm-4"B><"col-sm-4"i><"col-sm-4"p>>',
            lengthMenu =  list(c(10, 25, 50,-1),
                               c('10', '25', '50', 'All')),
            scrollCollapse = TRUE,
            lengthChange = TRUE,
            widthChange = TRUE,
            rownames = TRUE
          )
        )
      })

    # Server - UMAP Marker Expression Tab ----
    output$umap3 <-
      renderPlot(
        ggplot(dr_umap, aes_string(
          x = "u1",
          y = "u2",
          color = input$markerSelect
        )) +
          geom_point(size = 1.0) +
          theme_bw() +
          scale_color_gradientn(input$markerSelect,
                                colours = colorRampPalette(rev(brewer.pal(
                                  n = 11, name = "Spectral"
                                )))(50))
      )

    # Server - Differential Abundance Tab ----
    output$md_table <-
      renderTable(
        reactVals$md[1:5, ]
      )

    output$counts_table <-
      renderTable(
        reactVals$counts_table[1:5, 1:5]
      )

    output$DA_result_table <-
      renderTable(
        reactVals$DA_result_table
      )

    # Do the differential abundance ----
    observeEvent(input$doDA, {
      # browser()
      req(reactVals$counts_table, reactVals$md)

      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df01Tree$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL

      props_table <- t(t(countsTable) / colSums(countsTable)) * 100

      # mm <- match(md$sample_id, colnames(props_table))
      mm <- match(colnames(props_table), md$sample_id)

      tmp_cond <- md$condition[mm]

      DA_df <- data.frame()
      for (i in 1:nrow(props_table)) {
        foo <- pairwise.wilcox.test(as.numeric(props_table[i,]), tmp_cond, p.adjust.method=input$correction_method)

        df <- subset(melt(foo$p.value), value!=0)
        df$cell <- rownames(countsTable)[i]
        DA_df <- rbind(DA_df, as.data.frame(df))

      }

      # browser()
      # change the name of nodes:
      # if node has children add "_remaining"
      my_list <- lapply(DA_df$cell, function(x) {
            # get the id
            nid <- reactVals$graph$nodes$id[reactVals$graph$nodes$label == x]
            # check if the id has children
            if(nid %in% reactVals$graph$edges$to) {
              return (x <- paste0(x, "_remaining"))
            }
            else {
              return(x)
            }
          })

      DA_df$list <- unlist(my_list)
      
      colnames(DA_df) <- c("Cond1", "Cond2", "p-value", "Cell", "Naming")

      reactVals$DA_result_table <- DA_df

    })

    output$interactiveTree <- renderVisNetwork({
      visNetwork(reactVals$graph$nodes, reactVals$graph$edges, width = "100%") %>%
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    })

    myNode <- reactiveValues(selected = '')

    observeEvent(input$current_node_id, {
      myNode$selected <<- input$current_node_id
    })

    output$selectedNode <- renderText({myNode$selected})

    output$DA_interactive_table <- renderTable({

      # browser()
      children <- all_my_children(reactVals$graph, myNode$selected)

      if (!is.null(children)) {
        selected_labels <- c(reactVals$graph$nodes$label[children], reactVals$graph$nodes$label[myNode$selected])
      } else {
        selected_labels <- reactVals$graph$nodes$label[myNode$selected]
      }

      # browser()
      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df01Tree$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL

      # countsTable <- countsTable[selected_labels,]

      props_table <- t(t(countsTable) / colSums(countsTable)) * 100

      # check if any of the nodes is empty and therefore not in the props table
      # this cause the
      # make sure that sleecte labels are larger than on, otherwise it is
      # not a vector and treat it as single value
      if (length(selected_labels) > 1) {
        filter_vec <- selected_labels %in% rownames(props_table)
        props_table <- props_table[selected_labels[filter_vec],]
      }
      else {
        props_table <- t(as.data.frame(props_table[selected_labels,]))
      }

      ## subest the prps tabel based on the selected node
      # props_table <- props_table[selected_labels, ]
      props_table <- t(as.data.frame(colSums(props_table)))

      # mm <- match(md$sample_id, colnames(props_table))
      mm <- match(colnames(props_table), md$sample_id)

      tmp_cond <- md$condition[mm]

      DA_df <- data.frame()

      foo <- pairwise.wilcox.test(props_table, tmp_cond, p.adjust.method="none")

      df <- subset(melt(foo$p.value), value!=0)
      df$cell <- reactVals$graph$nodes$label[myNode$selected]
      DA_df <- rbind(DA_df, as.data.frame(df))

      colnames(DA_df) <- c("Var1", "Var2", "p-value", "Cell")
      reactVals$DA_interactive_table <- DA_df

    })

    output$boxplot <- renderPlot({

      # browser()
      children <- all_my_children(reactVals$graph, myNode$selected)

      if (!is.null(children)) {
        selected_labels <- c(reactVals$graph$nodes$label[children], reactVals$graph$nodes$label[myNode$selected])
      } else {
        selected_labels <- reactVals$graph$nodes$label[myNode$selected]
      }

      # browser()
      countsTable <- reactVals$counts_table
      ## aggregate the clusters by name:
      countsTable['cell'] <- df01Tree$cell
      # merge and aggregate by cell
      countsTable <- aggregate(. ~ cell, countsTable, sum)

      rownames(countsTable) <- countsTable$cell
      countsTable$cell <- NULL

      # countsTable <- countsTable[selected_labels,]

      props_table <- t(t(countsTable) / colSums(countsTable)) * 100

      # check if any of the nodes is empty and therefore not in the props table
      # this cause the
      # make sure that sleecte labels are larger than on, otherwise it is
      # not a vector and treat it as single value
      if (length(selected_labels) > 1) {
        filter_vec <- selected_labels %in% rownames(props_table)
        props_table <- props_table[selected_labels[filter_vec],]
      }
      else {
        props_table <- t(as.data.frame(props_table[selected_labels,]))
      }

      ## subest the prps tabel based on the selected node
      # props_table <- props_table[selected_labels, ]

      props_table <- as.data.frame(colSums(props_table))


      # mm <- match(md$sample_id, rownames(props_table))
      mm <- match(rownames(props_table), md$sample_id)

      props_table$cond <- md$condition[mm]

      names(props_table) <- c("value", "cond")

      ggplot(props_table, aes(x = cond, y = value, fill=cond))+
        geom_boxplot() +
        xlab("Condition") +
        ylab("Values") +
        ggtitle(reactVals$graph$nodes$label[myNode$selected])


    })

    # children <- graph$edges$from[graph$edges$to == myNode$selected]

    # observeEvent(input$btnLoadDemoData, {
    # 
    # })
    
   
# Load data based on pressing the unannotated data button
    observe({Settings_Server1(id="Settings",reactVals=reactVals)}) %>% 
      bindEvent(input$btnLoadDemoData)

# Load data based on pressing the annotated data button ----
    
    observe({Settings_Server2(id="Settings",reactVals=reactVals)}) %>% 
      bindEvent(input$btnLoadAnnoData)
    
# Load data based on pressing the Import Tree button  ----
    
    observeEvent(input$btnImportTree, {Settings_Server3(id="Settings")})
    
# Lead threshold Tab    ----
    
    observe(threshold_Server(id="threshold",reactVals))
    
# Run the Tree Annotation Tab based on pressing the Annotated demo button  ----
    observe({
      req(exists("cellfreq"))
      req(exists("df01"))
      TreeAnnotation_Server(id="TreeAnnotation",
                            reactVals=reactVals,
                            df01 = df01)}) %>%
      bindEvent(c(input$btnLoadAnnoData,input$btnLoadDemoData))

# Update Parent picker based on annotation list once Annotated demo is pressed ----
    
    observe({
      req(exists("annotationlist"))
      print("update Parent Picker")
      updateSelectizeInput(session, ("parentPicker"), choices = annotationlist, server = TRUE)
      })  %>%
    bindEvent(c(input$btnLoadAnnoData,input$btnLoadDemoData))
    
# Update Heatmap in Tree annotation tab based on parent picker ----
    
    observe({
      req(!is.null(reactVals$hm))
      Heatmap_Server(id="Heatmap",reactVals=reactVals,filter=input$parentPicker)}) %>% 
      bindEvent(input$parentPicker)
      
  }

  shinyApp(
    ui = ui,
    server = server
  )
# }

