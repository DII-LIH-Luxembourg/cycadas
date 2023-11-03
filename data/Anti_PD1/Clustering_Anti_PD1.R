

library(SummarizedExperiment)
library(Rtsne)
library(umap)
library(ggplot2)
library(diffcyt)

###########
# Load data
###########

library(HDCytoData)
d_SE <- Krieg_Anti_PD_1_SE()

###############
# Preprocessing
###############

# extract data
d_sub <- assay(d_SE[, colData(d_SE)$marker_class == "type"])
sample_ids <- rowData(d_SE)$sample_id
dim(d_sub)
stopifnot(nrow(d_sub) == length(sample_ids))

md <- metadata(d_SE)$experiment_info

md$batch_id <- NULL
colnames(md) <- c("condition", "sample_id")

write.csv(md, "./data/Anti_PD1/metadata.csv")

panel <- data.frame(colData(d_SE))
write.csv(panel, "./data/Anti_PD1/panel.csv")

# transform data
cofactor <- 5
d_sub <- asinh(d_sub / cofactor)
dim(d_sub)
summary(d_sub)



library(FlowSOM)
fsom <- ReadInput(d_sub, transform = FALSE, scale = FALSE)
# set.seed(566961715)
set.seed(42)
som <- BuildSOM(fsom, xdim=20, ydim=20, rlen=40)
## Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]
codes <- som$map$codes
# cell_clustering1 <- code_clustering1[cell_clustering_som]
cell_clustering1 <- cell_clustering_som

expr_median <- som$map$medianValues

# Calculate cluster frequencies
clustering_table <- as.numeric(table(som$map$mapping[,1]))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
sum(clustering_prop)
df_prop <- as.data.frame(clustering_prop)
df_prop$cluster <- rownames(df_prop)
head(as.data.frame(expr_median))
head(df_prop)

write.csv(expr_median, "./data/Anti_PD1/expr_median_400.csv", row.names = F)
write.csv(df_prop, "./data/Anti_PD1/cluster_freq_400.csv")

## ----------------------------------------------------------------------------
## Generate the Proportion Table
## ----------------------------------------------------------------------------
counts_table <- table(cell_clustering1, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
# counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(props, "./data/Anti_PD1/props_table_400.csv")






















