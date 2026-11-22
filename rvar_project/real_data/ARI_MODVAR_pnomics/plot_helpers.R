# plot_helpers.R
# Helper functions for MOD-VAR network decomposition figures
# Requires: igraph, purrr

# =========================================================================
# Single signed, directed network (circle layout)
# =========================================================================
plot_signed_graph <- function(A,
                              labels = colnames(A),
                              vertex.size = 30,
                              vertex.label.cex = 2,
                              edge.width = 2,
                              range = max(abs(A), na.rm = TRUE),
                              pal_floor = 0,
                              vertex.color = "white",
                              min_abs = 0,
                              margin = 0) {
  
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Please install the 'igraph' package.")
  if (!is.matrix(A)) stop("A must be a matrix.")
  if (nrow(A) != ncol(A)) stop("A must be square.")
  if (length(range) != 1 || !is.numeric(range) || range <= 0)
    stop("'range' must be a single positive number.")
  if (pal_floor < 0 || pal_floor > 1)
    stop("'pal_floor' must be between 0 and 1.")
  
  n <- nrow(A)
  if (is.null(labels)) labels <- paste0("X", seq_len(n))
  
  edges <- data.frame(from = integer(0), to = integer(0),
                      weight = numeric(0), sign = integer(0))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (abs(A[i, j]) > min_abs) {
        edges <- rbind(edges, data.frame(from = j, to = i,
                                         weight = abs(A[i, j]),
                                         sign = sign(A[i, j])))
      }
    }
  }
  
  g <- igraph::graph_from_data_frame(edges, directed = TRUE,
                                     vertices = data.frame(name = seq_len(n)))
  
  if (nrow(edges) > 0) {
    floor_col <- \(hi) colorRampPalette(c("#FFFFFF", hi))(101)[1 + round(pal_floor * 100)]
    blue_pal  <- colorRampPalette(c(floor_col("#093e78"), "#093e78"))(101)
    red_pal   <- colorRampPalette(c(floor_col("#A30707"), "#A30707"))(101)
    alpha <- pmin(edges$weight / range, 1)
    idx   <- pmax(1, round(alpha * 100) + 1)
    igraph::E(g)$color      <- ifelse(edges$sign > 0, blue_pal[idx], red_pal[idx])
    igraph::E(g)$width      <- 1 + 6 * alpha
    igraph::E(g)$arrow.size <- 0.2 + 2 * alpha
    igraph::E(g)$curved     <- 0.15
  }
  
  igraph::V(g)$label       <- labels
  igraph::V(g)$size        <- vertex.size
  igraph::V(g)$color       <- vertex.color
  igraph::V(g)$frame.color <- "black"
  igraph::V(g)$label.color <- "black"
  
  plot(g, layout = igraph::layout_in_circle(g),
       vertex.label.cex = vertex.label.cex, margin = margin)
}

# =========================================================================
# Internal: matrices for one row (GCN, Shared, moderators, [Individual])
# =========================================================================
decomp_mats <- function(model, mod_mat, mod_names, x,
                        include_individual = FALSE, drop_diag = TRUE) {
  zero_diag <- \(A) { if (drop_diag) diag(A) <- 0; A }
  mods  <- lapply(seq_along(mod_names),
                  \(m) zero_diag(mod_mat[x, m] * model$moderator_coeffs[[m]]))
  comps <- c(list(zero_diag(model$joint_coeffs)), mods)
  if (include_individual) comps <- c(comps, list(zero_diag(model$idiographic_coeffs[[x]])))
  c(list(zero_diag(model$bysubject_coeffs[[x]])), comps)
}

# Internal: color range across a set of rows
#   "gcn"   = largest |edge| among the GCNs (original behavior)
#   "all"   = largest |edge| in any panel
#   "panel" = NA; each panel scaled to its own max
decomp_range <- function(rows, range_by = c("gcn", "all", "panel")) {
  range_by <- match.arg(range_by)
  switch(range_by,
         gcn   = max(sapply(rows, \(r) max(abs(r[[1]]), na.rm = TRUE)), 1e-8),
         all   = max(sapply(rows, \(r) max(sapply(r, \(A) max(abs(A), na.rm = TRUE)))), 1e-8),
         panel = NA_real_
  )
}

# =========================================================================
# One row: GCN = comp1 + comp2 + ...
# =========================================================================
decomp_row <- function(mats, titles, first_row = FALSE, range, labels,
                       pal_floor = 0, vertex_cex = 2,
                       title_cex = 2, title_col = "#093e78", op_cex = 4) {
  for (i in seq_along(mats)) {
    if (i == 2) {
      plot.new(); text(0.5, 0.5, "=", cex = op_cex)
    } else if (i > 2) {
      plot.new(); text(0.5, 0.5, "+", cex = op_cex)
    }
    r <- if (is.na(range)) max(abs(mats[[i]]), 1e-8) else range
    plot_signed_graph(mats[[i]], labels = labels, range = r,
                      vertex.label.cex = vertex_cex, pal_floor = pal_floor)
    if (first_row) title(main = titles[i], col.main = title_col, cex.main = title_cex)
  }
}

# =========================================================================
# Full figure: one row per subject
# =========================================================================
plot_decomp <- function(model, mod_mat, mod_names, subjects,
                        labels = seq_len(nrow(model$joint_coeffs)),
                        include_individual = FALSE,
                        drop_diag = TRUE,
                        range_by = c("gcn", "all", "panel"),
                        pal_floor = 0) {
  
  rows   <- lapply(subjects, \(x) decomp_mats(model, mod_mat, mod_names, x,
                                              include_individual, drop_diag))
  titles <- c("GCN", "Shared", mod_names, if (include_individual) "Individual")
  range  <- decomp_range(rows, range_by)
  
  n_comp   <- length(titles) - 1
  n_panels <- 2 * n_comp + 1
  widths   <- c(4, rep(c(0.5, 4), n_comp))
  
  layout(matrix(seq_len(length(subjects) * n_panels),
                nrow = length(subjects), byrow = TRUE),
         widths = widths, heights = rep(1, length(subjects)))
  
  for (k in seq_along(subjects)) {
    decomp_row(rows[[k]], titles, first_row = (k == 1),
               range = range, labels = labels, pal_floor = pal_floor)
  }
}

# =========================================================================
# Check that shared + moderators + individual = bysubject for each subject
# =========================================================================
check_decomp <- function(model, mod_mat) {
  n <- length(model$bysubject_coeffs)
  stopifnot(nrow(mod_mat) == n)
  purrr::map_lgl(seq_len(n), \(x) {
    recon <- model$joint_coeffs +
      Reduce(`+`, lapply(seq_len(ncol(mod_mat)),
                         \(m) mod_mat[x, m] * model$moderator_coeffs[[m]])) +
      model$idiographic_coeffs[[x]]
    isTRUE(all.equal(recon, model$bysubject_coeffs[[x]], check.attributes = FALSE))
  })
}

# =========================================================================
# Draw and save full figure to PNG
# =========================================================================
save_decomp_png <- function(file, ..., n_subjects, width = 24,
                            row_height = 3.5, res = 300) {
  png(file, width = width, height = row_height * n_subjects, units = "in", res = res)
  on.exit(dev.off())
  par(mar = c(1, 0, 2, 0))
  plot_decomp(...)
}

# =========================================================================
# 16:9 slide: decomposition rows + key column
# =========================================================================
# range: NULL = computed from these subjects; pass a number (from
#        slide_range()) to share one color scale across several slides.
save_decomp_slide <- function(file, model, mod_mat, mod_names, subjects,
                              node_labels,
                              row_labels = NULL,
                              include_individual = FALSE,
                              title = NULL,
                              title_main_cex = 2.6,
                              title_height = 0.25,
                              drop_diag = TRUE,
                              range = NULL,
                              range_by = c("all", "gcn"),
                              pal_floor = 0.25,
                              width = 13.333, height = 7.5, res = 300,
                              title_cex = 2.2, vertex_cex = 1.8,
                              key_cex = 1.5, op_cex = 4) {
  
  range_by <- match.arg(range_by)
  rows   <- lapply(subjects, \(x) decomp_mats(model, mod_mat, mod_names, x,
                                              include_individual, drop_diag))
  titles <- c("GCN", "Shared", mod_names, if (include_individual) "Individual")
  if (is.null(range)) range <- decomp_range(rows, range_by)
  
  n_rows     <- length(subjects)
  n_comp     <- length(titles) - 1
  n_panels   <- 2 * n_comp + 1
  has_rowlab <- !is.null(row_labels)
  if (has_rowlab) stopifnot(length(row_labels) == n_rows)
  
  ## Layout: [row label] + panels + key column spanning all rows
  panel_ids <- matrix(seq_len(n_rows * n_panels), nrow = n_rows, byrow = TRUE)
  next_id   <- max(panel_ids)
  if (has_rowlab) {
    lab_ids   <- next_id + seq_len(n_rows)
    next_id   <- max(lab_ids)
    panel_ids <- cbind(lab_ids, panel_ids)
  }
  key_id <- next_id + 1
  lay    <- cbind(panel_ids, rep(key_id, n_rows))
  widths <- c(if (has_rowlab) 1.2, 4, rep(c(0.5, 4), n_comp), 5.5)
  heights <- rep(1, n_rows)
  
  has_title <- !is.null(title)
  if (has_title) {
    title_id <- key_id + 1
    lay      <- rbind(rep(title_id, ncol(lay)), lay)
    heights  <- c(title_height, heights)
  }
  
  png(file, width = width, height = height, units = "in", res = res, pointsize = 12)
  on.exit(dev.off())
  layout(lay, widths = widths, heights = heights)
  
  ## Panels
  par(mar = c(0.5, 0, 3, 0))
  for (k in seq_len(n_rows)) {
    decomp_row(rows[[k]], titles, first_row = (k == 1),
               range = range, labels = seq_along(node_labels),
               pal_floor = pal_floor, vertex_cex = vertex_cex,
               title_cex = title_cex, title_col = "black", op_cex = op_cex)
  }
  
  ## Row labels
  if (has_rowlab) {
    par(mar = c(0, 0, 0, 0))
    for (k in seq_len(n_rows)) {
      plot.new()
      text(0.5, 0.5, row_labels[k], srt = 90, cex = key_cex, font = 2)
    }
  }
  
  ## Key
  par(mar = c(0, 1, 3, 0))
  plot.new()
  n  <- length(node_labels)
  ys <- seq(0.92, 0.30, length.out = n)
  text(0, 1.0, "Nodes", adj = c(0, 0.5), cex = key_cex * 1.1, font = 2, xpd = NA)
  text(0.02, ys, seq_len(n), adj = c(0, 0.5), cex = key_cex, font = 2)
  for (i in seq_len(n)) {
    text(0.14, ys[i], node_labels[[i]], adj = c(0, 0.5), cex = key_cex)
  }
  
  segments(0.02, 0.18, 0.14, 0.18, col = "#093e78", lwd = 6)
  text(0.18, 0.18, "Positive lagged effect", adj = c(0, 0.5), cex = key_cex)
  segments(0.02, 0.08, 0.14, 0.08, col = "#A30707", lwd = 6)
  text(0.18, 0.08, "Negative lagged effect", adj = c(0, 0.5), cex = key_cex)
  
  ## Title
  if (has_title) {
    par(mar = c(0, 0, 0, 0))
    plot.new()
    text(0.5, 0.5, title, cex = title_main_cex, font = 2)
  }
  
  invisible(range)
}

# =========================================================================
# One color range across several slides
# =========================================================================
slide_range <- function(model, mod_mat, mod_names, subjects,
                        include_individual = FALSE, drop_diag = TRUE,
                        range_by = c("all", "gcn")) {
  range_by <- match.arg(range_by)
  rows <- lapply(subjects, \(x) decomp_mats(model, mod_mat, mod_names, x,
                                            include_individual, drop_diag))
  decomp_range(rows, range_by)
}

# =========================================================================
# Preview in a separate window sized for the figure (avoids RStudio pane limits)
# =========================================================================
preview_decomp <- function(..., n_subjects, width = 24, row_height = 3.5, res = 150) {
  graphics.off()
  f <- tempfile(fileext = ".png")
  png(f, width = width, height = row_height * n_subjects, units = "in", res = res)
  par(mar = c(1, 0, 2, 0))
  ok <- tryCatch({ plot_decomp(...); TRUE },
                 finally = dev.off())
  if (isTRUE(ok)) utils::browseURL(f)
  invisible(f)
}

# =========================================================================
# Grid of networks on a 16:9 slide, with a compact domain legend
# =========================================================================
save_network_grid <- function(file, mats, titles, nrow, ncol,
                              labels,
                              vertex_colors = "white",
                              legend_domains = NULL,     # named vector: label -> color
                              main = NULL,
                              range = NULL,
                              min_abs = 0,
                              pal_floor = 0.25,
                              width = 13.333, height = 7.5, res = 300,
                              main_cex = 2.6, panel_title_cex = 2.2,
                              vertex_cex = 1.4, vertex_size = 38,
                              legend_cex = 1.8, legend_width = 3,
                              margin = -0.05) {
  
  stopifnot(length(mats) == length(titles), length(mats) <= nrow * ncol)
  if (is.null(range)) range <- max(sapply(mats, \(A) max(abs(A), na.rm = TRUE)), 1e-8)
  
  # layout: panel grid + legend column (+ title row)
  ids <- matrix(seq_len(nrow * ncol), nrow = nrow, byrow = TRUE)
  leg_id <- nrow * ncol + 1
  lay    <- cbind(ids, rep(leg_id, nrow))
  widths <- c(rep(4, ncol), legend_width)
  heights <- rep(1, nrow)
  if (!is.null(main)) {
    lay     <- rbind(rep(leg_id + 1, ncol(lay)), lay)
    heights <- c(0.18 * nrow, heights)
  }
  
  png(file, width = width, height = height, units = "in", res = res, pointsize = 12)
  on.exit(dev.off())
  layout(lay, widths = widths, heights = heights)
  
  # panels
  par(mar = c(0.2, 0.2, 3, 0.2))
  for (p in seq_len(nrow * ncol)) {
    if (p <= length(mats)) {
      plot_signed_graph(mats[[p]], labels = labels, range = range,
                        vertex.size = vertex_size, vertex.label.cex = vertex_cex,
                        vertex.color = vertex_colors, pal_floor = pal_floor,
                        min_abs = min_abs, margin = margin)
      title(main = titles[p], cex.main = panel_title_cex)
    } else {
      plot.new()
    }
  }
  
  # legend
  par(mar = c(0, 0.5, 3, 0))
  plot.new()
  y <- 0.85
  if (!is.null(legend_domains)) {
    for (k in seq_along(legend_domains)) {
      points(0.06, y, pch = 21, bg = legend_domains[k], col = "black", cex = legend_cex * 2)
      text(0.18, y, names(legend_domains)[k], adj = c(0, 0.5), cex = legend_cex)
      y <- y - 0.09
    }
    y <- y - 0.06
  }
  segments(0.01, y, 0.12, y, col = "#093e78", lwd = 7)
  text(0.18, y, "Positive", adj = c(0, 0.5), cex = legend_cex)
  y <- y - 0.09
  segments(0.01, y, 0.12, y, col = "#A30707", lwd = 7)
  text(0.18, y, "Negative", adj = c(0, 0.5), cex = legend_cex)
  
  # title
  if (!is.null(main)) {
    par(mar = c(0, 0, 0, 0))
    plot.new()
    text(0.5, 0.5, main, cex = main_cex, font = 2)
  }
  
  invisible(range)
}
