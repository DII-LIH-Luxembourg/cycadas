![Alt text](./logo4_small2.png)

## Cytometry Cluster Annotation and Differential Abundance Suite

#### Efficient and reproducible cytometry data

Aims:

• facilitating the process of cluster annotation while reducing user bias • improving reproducibility

Key features:

• defining the threshold of positive/negative marker expression

• interactive inspection of cluster phenotypes

• automatic merging of populations

• differential abundance analysis

## Installation Instructions



``` r
library(devtools)

# Install all required packages
devtools::install_github("DII-LIH-Luxembourg/cycadas", dependencies = TRUE)

# start the cycadas shiny app
cycadas()
```

### Demo dataset

To enable tool exploration, we provide the demo dataset that can be loaded (**Load** tab → **Demo Data**) either as  cluster expression data only (**Load Cluster Expression Demo Data**, allowing the user to create the annotation) or as annotated data (**Load Annotated Demo Data** which include the annotation tree).

*This demo dataset is generated from the publicly available mass cytometry data of patients with idiopathic Parkinson’s disease and healthy controls (Capelle, C.M. et al., Nat Commun, 2023) that were clustered with GigaSOM to generate 1600 clusters.*

### Input data

All the data tables loaded into CyCadas are uploaded in the **Load** tab → **Required**. The tool requires output of the clustering algorithm in the form of 2 data tables:

• marker expression (mean or median expression of each marker in each cluster),

• cluster frequency (proportion of each cluster within the dataset).

## Prepare the Data input

### Median Expression and Cluster Frequencies from FlowSOM (R Code):

``` r
# within your clustering workflow create sample_ids according to the metadata files:

sample_ids <- rep(metadata$sample_id, fsApply(fcs, nrow))

library(FlowSOM)
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)

set.seed(42)
som <- BuildSOM(fsom, colsToUse = lineage_markers, xdim=20, ydim=20, rlen=40)

expr_median <- som$map$medianValues

# Calculate cluster frequencies
clustering_table <- as.numeric(table(som$map$mapping[,1]))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
df_prop <- as.data.frame(clustering_prop)
df_prop$cluster <- rownames(df_prop)

write.csv(expr_median, "expr_median.csv", row.names = F)
write.csv(df_prop, "cluster_freq.csv")

## ----------------------------------------------------------------------------
## Generate the Proportion Table
## ----------------------------------------------------------------------------
counts_table <- table(som$map$mapping[,1], sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100

props <- as.data.frame.matrix(props_table)

write.csv(props, "proportion_table.csv")
```

### Median Expression and Cluster Frequencies from GigaSOM (Julia Code):

``` julia
Random.seed!(42)

som = initGigaSOM(di, gridSize, gridSize, seed = 13) # set a seed value
som = trainGigaSOM(som, di, epochs = nEpochs)

mapping_di = mapToGigaSOM(som, di)
mapping = gather_array(mapping_di)

# get the cluster frequencies
clusterFreq = dcount(mc, mapping_di)
df = DataFrame(cluster = 1:length(df), clustering_prop = clusterFreq)
df.clustering_prop = df.clustering_prop ./ sum(df.clustering_prop)
CSV.write("cluster_freq_$mc.csv", df)

files = distributeFCSFileVector(:fileIDs, md[:, :file_name])

# Get the count table per fileID
count_tbl = dcount_buckets(mc, mapping_di, size(md, 1), files)
ct = DataFrame(count_tbl, :auto)
# ct = convert(DataFrame, count_tbl)
rename!(ct, md.sample_id)
# export the count talbe
CSV.write("cluster_counts_$mc.csv", ct)

# Get the median expression per cluster
expr_tbl = dmedian_buckets(di, mc, mapping_di, cols)

et  = DataFrame(expr_tbl, :auto)
rename!(et, lineage_markers)

# export median marker expression
CSV.write("median_expr_$mc.csv", et)
```

Detailed workflow for each method can be found in the data section.

### Data exploration

The **UMAP interactive** tab allows the preview of marker expression in the clusters selected by the user on the UMAP:

![](./www/umap_interactive.png)

In the **UMAP Marker expression** tab, user can investigate the expression level of the selected marker across all the clusters.

![](./www/umap_marker_expression.png)

### Thresholds

In the **Thresholds** tab, the estimation of threshold value defining negative and positive marker expression of each marker is based on 1-dimensional k-means clustering. The bimodality for every marker is assessed and the bimodal coefficient values are reported. The blue threshold line indicates that data meets the bimodal distribution criteria, otherwise it is colored red. The threshold value can be manually adjusted by clicking on the scatterplot.

*Expression of CD8a with blue threshold line indicating the bimodal distribution:*

![](./www/thresholds_bimodal_cd8.png)

*Expression of TCRgd with red threshold line indicating that this marker expression does not follow the bimodal distribution:*

![](./www/thresholds_notbimodal_tcrgd.png)

### Annotation

The **Annotation** tab allows performing the annotation in a tree-based hierarchical process - initially, the main cell types are defined, followed by the identification of their subtypes (with the level of detail defined by the user).

All the clusters are initially defined as "unassigned". Then, upon the selection of positive and negative markers defining the population, clusters characterized by given expression pattern are re-assigned from the parent node to the child node.

*Scheme depicting the process of building the annotation tree:*

![](./www/annotation_building_tree.png)

*Cropped fragment of the completed annotation tree:*

![](./www/annotation_tree.png)

Upon selection of the node, heatmap displaying the expression of all the markers in all the clusters belonging to this node is shown.

*Heatmap depicting phenotype of clusters annotated as CD8+ TEM cells:*

![](./www/annotation_heatmap.png)

### Differential abundance analysis

In the **Differential Abundance** tab, a pairwise Wilcoxon test on all the nodes is performed upon selecting the desired multiple testing correction method:

![](./www/differential_abundance.png)

**DA Interactive Tree** allows exploration of abundance of all the defined subpopulations across the conditions by selecting the node on the annotation tree.

*Upon clicking on the desired node...*

![](./www/DA_interactive_tree_tree.png)

*... proportion of the selected celltype across the condition is plotted.*

![](./www/DA_interactive_tree_plot.png)

### Data export

Differential abundance analysis results, as well as proportion table (% of defined cell populations across all the samples) can be exported in the **Differential Abundance** tab.

![](./www/differential_abundance_export.png)

Files enabling the continuation of the analysis - modified threshold values, as well as annotation tree structure, can be exported from the **Thresholds** and **Annotation** tabs, respectively, and re-loaded (**Load** tab) to continue the analysis.

*Exporting annotation tree:*

![](./www/annotation_export.png)

*Exporting threshold values:*

![](./www/thresholds_export.png)
