![Alt text](./cycadas_logo.png)

## Cytometry Cluster Annotation and Differential Abundance Suite

#### Efficient and reproducible cytometry data analysis

## Installation Instructions

To run: - Checkout branch develop - load all required settings by: pkgload::load_all(".") - call cycadas() to run

## Working with Cycadas

### Starting a new project and exploring data:

**Cycadas** requires two CSV files being an output of clustering algorithm -- 1) median expression of markers in the clusters and 2) frequency of each cluster (Load -- Required). To have a first look on the data, explore the UMAP interactive and UMAP Marker expression tabs. In UMAP interactive tab, select the desired part of the UMAP to visualise marker expression in selected clusters on the heatmap. In the UMAP Marker expression, select the marker of interest to evaluate its expression in all the clusters.

### Selecting threshold for each marker:

As the cluster annotation (phenotype definition) relies on a defined set of negative/positive markers, each marker should follow a bi-modal distribution. An initial separation (threshold) between positve (high) and negative (low) expression is estimated after uploading the data-set and indicated by a blue vertical line in the scatter plot and histogram. If a marker doesn't show a bimodal distribution, it's marked with a red line. In such cases the adjustments have to be made manually, or this marker may not be considered as phenotype marker.

The threshold value can be adjusted by clicking at the desired position within the scatter plot. In such cases any already defined annotation is automatically re-calculated.

Threshold settings can be exported as csv and re-used by loading.

### Annotation

Select the desired combination of positive and negative markers in the rolling menu. The classification of the clusters as negative or positive for a specific marker is occurring based on the threshold values set in the [Thresholds](#thresholds) tab. Clusters characterised by the defined (negative/positive) expression for the chosen markers are then selected, and:

• heatmap displays normalized expression of all the markers in the selected clusters

• UMAP represents all the clusters, with the selected clusters being highlighted in blue
