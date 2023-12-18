
DA_interactive_UI <- function(id) {
  ns <- NS(id)
  column(width = 4,
         box(
           width = NULL,
           tableOutput(ns("DA_interactive_table"))
         ),
         box(
           width = NULL,
           plotOutput(ns("boxplot"))
         )
         )

  }

DA_interactive_Server <- function(id,current_node_id) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      print("Run DA interactive Tab")

      output$selectedNode <- renderText({current_node_id})
      
      output$DA_interactive_table <- renderTable({
        
        children <- all_my_children(session$userData$vars$graph, current_node_id)
        
        if (!is.null(children)) {
          selected_labels <- c(session$userData$vars$graph$nodes$label[children], session$userData$vars$graph$nodes$label[current_node_id])
        } else {
          selected_labels <- session$userData$vars$graph$nodes$label[current_node_id]
        }

        countsTable <- session$userData$vars$counts_table
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
        df$cell <- session$userData$vars$graph$nodes$label[current_node_id]
        DA_df <- rbind(DA_df, as.data.frame(df))
        
        colnames(DA_df) <- c("Var1", "Var2", "p-value", "Cell")
        session$userData$vars$DA_interactive_table <- DA_df
        
      })

      output$boxplot <- renderPlot({

        # browser()
        children <- all_my_children(session$userData$vars$graph, current_node_id)

        if (!is.null(children)) {
          selected_labels <- c(session$userData$vars$graph$nodes$label[children], session$userData$vars$graph$nodes$label[current_node_id])
        } else {
          selected_labels <- session$userData$vars$graph$nodes$label[current_node_id]
        }

        # browser()
        countsTable <- session$userData$vars$counts_table
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
          ggtitle(session$userData$vars$graph$nodes$label[current_node_id])


      })
      
    })}