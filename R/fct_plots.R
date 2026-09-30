# Plots -----------------------------------------------------------------------
# Plot builders used by the modules. They take plain data and return a
# ggplot object, or draw directly for base / pheatmap graphics.

COLOR_SELECTED <- "#2c7be5"
COLOR_OTHER <- "#d0d5dd"

theme_cycadas <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#eef0f3"),
          axis.title = element_text(color = "#495057"),
          legend.position = "bottom")
}

# Empty plot with a message
plot_message <- function(msg = "No data available") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, size = 5, color = "#6c757d") +
    theme_void()
}

# Heatmap of cluster expression ------------------------------------------------
plot_cluster_heatmap <- function(hm) {
  if (is.null(hm)) {
    return(plot_message("No data available"))
  }
  if (nrow(hm) == 0) {
    return(plot_message("No clusters in this selection"))
  }
  pheatmap(hm, cluster_cols = F, cluster_rows = nrow(hm) >= 2,
           show_rownames = nrow(hm) <= 50, border_color = NA, fontsize = 11,
           treeheight_row = 30, silent = FALSE)
}

# UMAP plots ------------------------------------------------------------------
plot_umap <- function(umap) {
  ggplot(umap, aes(x = u1, y = u2)) +
    geom_point(size = 1.2, color = COLOR_SELECTED, alpha = .7) +
    theme_cycadas()
}

plot_umap_marker <- function(umap, marker) {
  ggplot(umap, aes(x = u1, y = u2, color = .data[[raw_markers(marker)]])) +
    geom_point(size = 1.2) +
    theme_cycadas() +
    theme(legend.position = "right") +
    scale_color_gradientn(marker,
                          colours = colorRampPalette(rev(brewer.pal(
                            n = 11, name = "Spectral"
                          )))(50))
}

# `selection` is the vector returned by filterColor()
plot_umap_selection <- function(umap, selection) {
  cbind(umap, ClusterSelection = selection) %>%
    dplyr::mutate(ClusterSelection = factor(ClusterSelection,
                                            levels = c("other clusters", "selected phenotype"))) %>%
    arrange(desc(ClusterSelection)) %>%
    ggplot(aes(x = u1, y = u2, color = ClusterSelection)) +
    geom_point(size = 1.2) +
    scale_color_manual(NULL, values = c("other clusters" = COLOR_OTHER,
                                        "selected phenotype" = COLOR_SELECTED),
                       labels = c("other clusters" = "Other clusters",
                                  "selected phenotype" = "Selection"),
                       drop = FALSE) +
    theme_cycadas() +
    guides(color = guide_legend(override.aes = list(size = 4)))
}

# Threshold plots -------------------------------------------------------------
# `points` is data.frame(x, y): the 0-1 marker expression and a random jitter
# `color` is the threshold table color: "blue" bimodal, "red" not bimodal.
# When `estimated` is given and differs, it is drawn as a thin grey line; both
# lines are labelled, each label on the side facing away from the other line.
COLOR_ESTIMATE <- "#868e96"

threshold_vline <- function(threshold, color, estimated = NA) {
  col <- if (identical(color, "red")) "#e8590c" else COLOR_SELECTED
  layers <- list(geom_vline(xintercept = threshold, linetype = "dashed", linewidth = 1, color = col))
  if (!is.na(estimated) && abs(threshold - estimated) > 5e-4) {
    right <- threshold > estimated
    layers <- c(
      geom_vline(xintercept = estimated, linewidth = .6, color = COLOR_ESTIMATE),
      annotate("text", x = estimated, y = Inf, label = "estimate", vjust = 1.5,
               hjust = if (right) 1.1 else -0.1, size = 3.5, color = COLOR_ESTIMATE),
      layers,
      annotate("text", x = threshold, y = Inf, label = "threshold", vjust = 1.5,
               hjust = if (right) -0.1 else 1.1, size = 3.5, color = col, fontface = "bold")
    )
  }
  layers
}

plot_threshold_scatter <- function(points, threshold, color, estimated = NA) {
  ggplot(points, aes(x = x, y = y, color = x >= threshold)) +
    geom_point(size = 1.2, alpha = .7) +
    scale_color_manual(values = c(`FALSE` = "#adb5bd", `TRUE` = COLOR_SELECTED), guide = "none") +
    theme_cycadas() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    labs(x = "Expression (scaled 0 to 1)") +
    threshold_vline(threshold, color, estimated)
}

plot_threshold_histogram <- function(points, threshold, color, estimated = NA) {
  ggplot(points, aes(x = x)) +
    geom_histogram(bins = 80, fill = "#adb5bd") +
    labs(x = "Expression (scaled 0 to 1)", y = "Clusters") +
    theme_cycadas() +
    threshold_vline(threshold, color, estimated)
}

# DA boxplot ------------------------------------------------------------------
# `props` is the data.frame(value, cond) returned by node_proportions()
plot_da_boxplot <- function(props, title) {
  ggplot(props, aes(x = cond, y = value, fill = cond)) +
    geom_boxplot(outlier.shape = NA, alpha = .6, width = .6) +
    geom_jitter(width = 0.15, size = 1.8, alpha = .8) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Proportion (%)") +
    geom_signif(
      comparisons = get_pairs(levels(props$cond)),
      map_signif_level = F,
      textsize = 4.5,
      step_increase = 0.1
    ) +
    theme_cycadas() +
    theme(legend.position = "none")
}
