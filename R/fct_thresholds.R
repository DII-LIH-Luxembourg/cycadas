# Marker thresholds -----------------------------------------------------------
# The threshold table ("th") has one row per marker and the columns
#   cell       marker name (also used as rowname)
#   threshold  cut-off on the 0-1 scale
#   color      "blue" for bimodal markers, "red" otherwise
#   bi_mod     bimodality coefficient
#   estimated  threshold estimated by kmeansTH(), kept when the threshold is
#              changed by hand; NA if unknown

# Bimodality cut-off as stated in Pfister et al., 2013
BIMODALITY_CUTOFF <- 0.555

# Threshold of one marker: the midpoint of two groups from 1-D k-means or a
# two-component Gaussian mixture, whichever separates the values better
# (silhouette). A marker with fewer than 2 distinct values gets 0.5; when the
# mixture can't be fitted (e.g. all clusters but one at the same value), the
# k-means midpoint is used.
estimate_threshold <- function(x) {
  if (length(unique(x)) < 2) return(0.5)

  set.seed(42)
  z <- Ckmeans.1d.dp(x, 2)
  midpoint_kmeans <- mean(z$centers)
  kmeans_avg_silhouette <- mean(silhouette(z$cluster, dist(x))[, 3])

  gmm <- tryCatch({
    fit <- Mclust(x, G = 2, verbose = FALSE)
    if (is.null(fit) || length(unique(fit$classification)) < 2) NULL
    else list(midpoint = mean(fit$parameters$mean),
              silhouette = mean(silhouette(fit$classification, dist(x))[, 3]))
  }, error = function(e) NULL)

  if (is.null(gmm) || !is.finite(gmm$midpoint) || kmeans_avg_silhouette >= gmm$silhouette) {
    midpoint_kmeans
  } else {
    gmm$midpoint
  }
}

# Estimate the threshold values -----------------------------------------------
kmeansTH <- function(df, th_mode="km") {
  th <- data.frame(cell = colnames(df), threshold = 0.0, color = "blue", bi_mod = 0)

  for (m in th$cell) {
    # check for bi-modal distribution, if not color = red to indicate
    bi_mod_value <- bimodality_coefficient(df[, m])
    if (!is.finite(bi_mod_value)) bi_mod_value <- 0   # constant marker
    th[th$cell == m, "bi_mod"] <- round(bi_mod_value, 3)
    if (bi_mod_value < BIMODALITY_CUTOFF) {
      th[th$cell == m, "color"] <- "red"
    }

    th[th$cell == m, "threshold"] <- estimate_threshold(df[, m])
  }
  th$estimated <- th$threshold
  rownames(th) <- th$cell
  return(th)
}

# Re-estimate the thresholds with a specific method ---------------------------
# Currently not reachable from the UI.
updateTH <- function(df, th, th_mode) {

  set.seed(42)

  for (m in th$cell) {
    set.seed(42)

    th[th$cell == m, "threshold"] <- switch(
      th_mode,
      km = {
        z <- Ckmeans.1d.dp(df[, m], 2)
        round(ahist(z, style="midpoints", data=df[, m], plot=FALSE)$breaks[2:2], 3)
      },
      gmm_mid = {
        fit <- normalmixEM(df[, m], k = 2)
        round(mean(fit$mu), 3)
      },
      gmm_high = {
        fit <- normalmixEM(df[, m], k = 2)
        round(fit$mu[2] - fit$sigma[2], 3)
      },
      gmm_low = {
        fit <- normalmixEM(df[, m], k = 2)
        round(fit$mu[1] + fit$sigma[1], 3)
      },
      mclust = {
        model <- Mclust(df[, m], G = 2)
        mean(model$parameters$mean)
      },
      kde = {
        dens <- density(asinh(df[, m]))
        # Fit Gaussian Mixture Model (GMM) to the peaks and take the
        # midpoint between the two cluster means
        peaks <- data.frame(x = dens$x, y = dens$y)
        fit <- Mclust(peaks, G = 2)
        mean(sort(fit$parameters$mean))
      },
      stop("Unknown threshold method: ", th_mode)
    )
  }
  return(th)
}

# Prepare an uploaded threshold table -----------------------------------------
prepare_thresholds <- function(th) {

  th$X <- NULL
  th$color <- ifelse(th$bi_mod < BIMODALITY_CUTOFF, "red", "blue")
  if (is.null(th$estimated)) th$estimated <- NA_real_
  rownames(th) <- th$cell

  return(th)
}

# Fill in missing estimates -------------------------------------------------
# Thresholds from a CSV or an older workspace have no estimate; run
# kmeansTH() for those markers.
add_estimates <- function(th, expr, markers) {

  if (is.null(th$estimated)) th$estimated <- NA_real_
  missing <- intersect(th$cell[is.na(th$estimated)], markers)
  if (length(missing) > 0) {
    est <- kmeansTH(expr[, missing, drop = FALSE])
    th$estimated[match(missing, th$cell)] <- est$threshold[match(missing, est$cell)]
  }
  th
}

# TRUE for markers whose threshold differs from the estimate
threshold_modified <- function(th) {
  !is.na(th$estimated) & abs(th$threshold - th$estimated) > 5e-4
}

# Set the thresholds of `cells` back to their estimates
reset_thresholds <- function(th, cells = th$cell) {
  idx <- th$cell %in% cells & !is.na(th$estimated)
  th$threshold[idx] <- th$estimated[idx]
  th
}

# Threshold table as written by the download and read by prepare_thresholds()
write_thresholds <- function(th, file) {
  write.csv(th[, c("cell", "threshold", "estimated", "bi_mod")], file)
}
