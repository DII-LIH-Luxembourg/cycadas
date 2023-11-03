# cycadas

### ToRun:
pkgload::load_all(".")

cycadas()
=======
![Alt text](./cycadas_logo.png)

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

To run:

-   Checkout branch develop

<!-- -->

-   load all required settings by: pkgload::load_all(".")

-   call cycadas() to run

## Working with CyCadas

### Starting a new project and exploring data:

**CyCadas** requires two CSV files being an output of clustering algorithm -- 1) median expression of markers in the clusters and 2) frequency of each cluster (Load -- Required). To have a first look on the data, explore the UMAP interactive and UMAP Marker expression tabs. In UMAP interactive tab, select the desired part of the UMAP to visualise marker expression in selected clusters on the heatmap. In the UMAP Marker expression, select the marker of interest to evaluate its expression in all the clusters.

### Selecting threshold for each marker:

As the cluster annotation (phenotype definition) relies on a defined set of negative/positive markers, each marker should follow a bi-modal distribution. An initial separation (threshold) between positve (high) and negative (low) expression is estimated after uploading the data-set and indicated by a blue vertical line in the scatter plot and histogram. If a marker doesn't show a bimodal distribution, it is marked with a red line. In such cases the adjustments have to be made manually, or this marker should be re-considered as phenotypic marker.

The threshold value can be adjusted by clicking at the desired position within the scatter plot. In such cases any already defined annotation is automatically re-calculated.

Threshold settings can be exported as csv and re-used by loading.

### Annotation

Select the desired combination of positive and negative markers in the rolling menu. The classification of the clusters as negative or positive for a specific marker is occurring based on the threshold values set in the Thresholds tab. Clusters characterised by the defined (negative/positive) expression for the chosen markers are then selected, and:

• heatmap displays normalized expression of all the markers in the selected clusters, • UMAP represents all the clusters, with the selected clusters being highlighted in blue.

A new node containing clusters characterised by a selected phenotype will be generated.

Annotation is a hierarchical process, where desired sub-populations can be differentiated from the main cell phenotypes.

### Differential abundance analysis

Upon uploading: • metadata table (defining the condition of each sample, for example patient vs healthy), • proportion table (defining the cluster composition of each sample), differential abundance analysis (comparing cellular composition in samples of different conditions) is performed (pairwise Wilcoxon test). Multiple testing correction methods (Bonferroni, Hochberg, FDR, among others) can be selected.

### Demo dataset

CyCadas contains not annotated and annotated demo datasets of clustered mass cytometry data from acute infected infividuals and uninfected controls (<https://doi.org/10.1016/j.celrep.2022.110815>) that can be loaded to enable tool exploration.
