# Plots -----------------------------------------------------------------------
# Plot builders used by the modules. They take plain data and return a
# ggplot object, or draw directly for base / pheatmap graphics.

theme_cycadas <- function(legend_text = 12, axis_text = 12, axis_title = 20) {
  theme_classic() +
    theme(legend.text = element_text(size = legend_text),
          legend.title = element_text(size = 20),
          axis.text = element_text(size = axis_text),
          axis.title = element_text(size = axis_title))
}

# Empty plot with a message
plot_message <- function(msg = "No Data Available") {
  plot(1, type = "n", main = msg)
}

# Heatmap of cluster expression ------------------------------------------------
plot_cluster_heatmap <- function(hm) {
  if (is.null(hm)) {
    plot_message("No Data Available")
  } else if (nrow(hm) == 0) {
    plot_message("Node is empty!")
  } else if (nrow(hm) < 2) {
    pheatmap(hm, cluster_cols = F, cluster_rows = F)
  } else {
    pheatmap(hm, cluster_cols = F)
  }
}

# UMAP plots ------------------------------------------------------------------
plot_umap <- function(umap) {
  ggplot(umap, aes(x = u1, y = u2)) +
    geom_point(size = 1.0) +
    theme_cycadas(legend_text = 8) +
    guides(color = guide_legend(override.aes = list(size = 4)))
}

plot_umap_marker <- function(umap, marker) {
  ggplot(umap, aes(x = u1, y = u2, color = .data[[raw_markers(marker)]])) +
    geom_point(size = 1.0) +
    theme_bw() +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 20),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20)) +
    scale_color_gradientn(marker,
                          colours = colorRampPalette(rev(brewer.pal(
                            n = 11, name = "Spectral"
                          )))(50))
}

# `selection` is the vector returned by filterColor()
plot_umap_selection <- function(umap, selection) {
  cbind(umap, ClusterSelection = selection) %>%
    dplyr::mutate(ClusterSelection = fct_relevel(as.factor(ClusterSelection),
                                                 "other clusters", "selected phenotype")) %>%
    arrange(desc(ClusterSelection)) %>%
    ggplot(aes(x = u1, y = u2, color = ClusterSelection)) +
    geom_point(size = 1.0) +
    theme_cycadas() +
    guides(color = guide_legend(override.aes = list(size = 4)))
}

# Threshold plots -------------------------------------------------------------
# `points` is data.frame(x, y): the 0-1 marker expression and a random jitter
threshold_vline <- function(threshold, color) {
  geom_vline(xintercept = threshold, linetype = "dotted", color = color, linewidth = 1.5)
}

plot_threshold_scatter <- function(points, threshold, color) {
  ggplot(points, aes(x = x, y = y)) +
    geom_point(size = 1) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 18),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    labs(x = "Scale 0 to 1") +
    threshold_vline(threshold, color)
}

plot_threshold_histogram <- function(points, threshold, color) {
  ggplot(points, aes(x = x)) +
    geom_histogram(bins = 80) +
    labs(x = "Scale 0 to 1") +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 18)) +
    threshold_vline(threshold, color)
}

# DA boxplot ------------------------------------------------------------------
# `props` is the data.frame(value, cond) returned by node_proportions()
plot_da_boxplot <- function(props, title) {
  ggplot(props, aes(x = cond, y = value, fill = cond)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    xlab("Condition") +
    ylab("Proportion") +
    geom_signif(
      comparisons = get_pairs(levels(props$cond)),
      map_signif_level = F,
      textsize = 6,
      step_increase = 0.1
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 22),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.position = "none") +
    ggtitle(title)
}
