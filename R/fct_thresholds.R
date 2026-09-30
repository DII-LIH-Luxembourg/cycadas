# Marker thresholds -----------------------------------------------------------
# The threshold table ("th") has one row per marker and the columns
#   cell       marker name (also used as rowname)
#   threshold  cut-off on the 0-1 scale
#   color      "blue" for bimodal markers, "red" otherwise
#   bi_mod     bimodality coefficient

# Bimodality cut-off as stated in Pfister et al., 2013
BIMODALITY_CUTOFF <- 0.555

# Estimate the threshold values -----------------------------------------------
kmeansTH <- function(df, th_mode="km") {
  th <- data.frame(cell = colnames(df), threshold = 0.0, color = "blue", bi_mod = 0)

  for (m in th$cell) {
    # check for bi-modal distribution, if not color = red to indicate
    bi_mod_value <- bimodality_coefficient(df[, m])
    th[th$cell == m, "bi_mod"] <- round(bi_mod_value, 3)
    if (bi_mod_value < BIMODALITY_CUTOFF) {
      th[th$cell == m, "color"] <- "red"
    }

    ######## K-Means
    set.seed(42)

    z <- Ckmeans.1d.dp(df[, m], 2)
    midpoint_kmeans <- mean(z$centers)

    kmeans_silhouette <- silhouette(z$cluster, dist(df[, m]))
    kmeans_avg_silhouette <- mean(kmeans_silhouette[, 3])

    ###########################################################################
    # Perform GMM clustering ##################################################
    gmm_result <- Mclust(df[, m], G = 2)

    # Midpoint between the means of the Gaussian components
    midpoint_gmm <- mean(gmm_result$parameters$mean)
    gmm_silhouette <- silhouette(gmm_result$classification, dist(df[, m]))
    gmm_avg_silhouette <- mean(gmm_silhouette[, 3])
    ###########################################################################

    if (kmeans_avg_silhouette >= gmm_avg_silhouette) {
      th[th$cell == m, "threshold"] <- midpoint_kmeans
    } else {
      th[th$cell == m, "threshold"] <- midpoint_gmm
    }
  }
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
  rownames(th) <- th$cell

  return(th)
}
