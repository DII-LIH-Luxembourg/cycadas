# Development

## New features
- `cycadas(dir)` opens the cluster data in a folder (`median_expr.csv`,
  `cluster_freq.csv`, optionally `cluster_counts.csv` and `metadata.csv`) and
  adds a **Save** button that writes `workspace.rds` and `labels.csv`
  (cluster, label, path in the tree) into it; a saved workspace there is
  reopened. Used by CytoSuite to run CyCadas on its clusterings.

# Version 1.4.0

## Bug fixes
- Fixed an issue where uploading an invalid or incorrectly formatted data file could cause the app to crash.

## New features
- Added **import function** for *remotesom* data.
- Added a **“Clear Workspace”** button to reset the app state.
- Added **JavaScript enhancements** to improve the file upload experience.


# Version 1.3.0

## Changes
- Merged the UMAP exploration views into a single page: **“Explore Data”**.


# Version 1.2.0

## New features
- Added a **Workspace Status Box** to display current session information.


# Version 1.1.0

## New features
- Added buttons for **saving and loading** the workspace state.
