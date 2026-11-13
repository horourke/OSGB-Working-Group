

##############################################################
##############################################################
##############################################################

plot_signed_graph <- function(A,
                              labels = colnames(A),
                              vertex.size = 30,
                              edge.width = 2) {
  
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Please install the 'igraph' package.")
  
  if (!is.matrix(A))
    stop("A must be a matrix.")
  
  if (nrow(A) != ncol(A))
    stop("A must be square.")
  
  n <- nrow(A)
  
  if (is.null(labels))
    labels <- paste0("X", seq_len(n))
  
  ## Build edge list
  edges <- data.frame(
    from   = integer(0),
    to     = integer(0),
    weight = numeric(0),
    sign   = integer(0)
  )
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (A[i, j] != 0) {
        edges <- rbind(
          edges,
          data.frame(
            from   = j,
            to     = i,
            weight = abs(A[i, j]),
            sign   = sign(A[i, j])
          )
        )
      }
    }
  }
  
  g <- igraph::graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = data.frame(name = labels)
  )
  
  ## ---------- Edge colors ----------
  w <- edges$weight
  
  # Normalize magnitudes to [0,1]
  if (max(w) == min(w)) {
    alpha <- rep(1, length(w))
  } else {
    alpha <- (w - min(w)) / (max(w) - min(w))
  }
  
  # Create palettes
  blue_pal <- colorRampPalette(c("#DCEEFF", "#0000CC"))(101)
  red_pal  <- colorRampPalette(c("#FFDCDC", "#CC0000"))(101)
  
  idx <- pmax(1, round(alpha * 100) + 1)
  
  edge_cols <- ifelse(edges$sign > 0,
                      blue_pal[idx],
                      red_pal[idx])
  
  ## ---------- Graph attributes ----------
  igraph::E(g)$color <- edge_cols
  
  # Optional: scale width with magnitude as well
  igraph::E(g)$width <- 1 + 4 * alpha
  
  igraph::E(g)$arrow.size <- 1.2
  igraph::E(g)$curved <- 0.15
  
  igraph::V(g)$size <- vertex.size
  igraph::V(g)$color <- "white"
  igraph::V(g)$frame.color <- "black"
  igraph::V(g)$label.color <- "black"
  
  plot(
    g,
    layout = igraph::layout_in_circle(g),
    vertex.label.cex = 2
  )
}


##############################################################
##############################################################
##############################################################
load("711_Data.RData")

library(tseries)
library(MTS)
library(ppcor)
library(TSA)
library(forecast)
library(ggplot2)

library(magrittr)
library(readxl)
library(tidyverse)
library(plot.matrix)

library(glmnet)
library(mvtnorm)
library(multivar)
library(BigVAR)
library(expm)
library(gridExtra)


mod_coeffs <- lapply(
  cv.model$moderator_coeffs, 
  function(x) {
    colnames(x) <- colnames(dfn_list[[1]])
    rownames(x) <- colnames(dfn_list[[1]])
    return(x)
  }
)
par()$mar
par(mfrow = c(1,4),
    mar = c(2.1, 2.1, 4.1, 2.1))

labels <- colnames(dfn_list[[1]])
mod_names <- 
  c("Age", "Gender", "LC Volume", "LC-CNR")

lapply(1:4,
  function(x) {
    plot_signed_graph(
      mod_coeffs[[x]],
      labels = 1:9,
      vertex.size = 30,
      edge.width = 2) 
    title(main = mod_names[x],
          col.main = "blue", 
          cex.main = 2)
  }
)





par(mfrow = c(2,2),
    mar = c(2.1, 2.1, 4.1, 2.1))

labels <- colnames(dfn_list[[1]])
mod_names <- 
  c("Age", "Biological Sex", "LC-Volume", "LC-CNR")

lapply(1:4,
       function(x) {
         plot_signed_graph(
           mod_coeffs[[x]],
           labels = 1:9,
           vertex.size = 30,
           edge.width = 3) 
         title(main = mod_names[x],
               col.main = "blue", 
               cex.main = 2)
       }
)












##############################################################
##############################################################
##############################################################

plot_signed_graph <- function(A,
                              labels = colnames(A),
                              vertex.size = 30,
                              vertex.label.cex = 2,
                              edge.width = 2,
                              range = max(abs(A), na.rm = TRUE)) {
  
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Please install the 'igraph' package.")
  
  if (!is.matrix(A))
    stop("A must be a matrix.")
  
  if (nrow(A) != ncol(A))
    stop("A must be square.")
  
  if (length(range) != 1 || !is.numeric(range) || range <= 0)
    stop("'range' must be a single positive number.")
  
  n <- nrow(A)
  
  if (is.null(labels))
    labels <- paste0("X", seq_len(n))
  
  ## Build edge list
  edges <- data.frame(
    from   = integer(0),
    to     = integer(0),
    weight = numeric(0),
    sign   = integer(0)
  )
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (A[i, j] != 0) {
        edges <- rbind(
          edges,
          data.frame(
            from   = j,
            to     = i,
            weight = abs(A[i, j]),
            sign   = sign(A[i, j])
          )
        )
      }
    }
  }
  
  g <- igraph::graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = data.frame(name = labels)
  )
  
  ## ---------- Edge colors ----------
  ##
  ## Map -range -> red
  ##      0      -> white
  ##     +range -> blue
  ##
  
  if (nrow(edges) > 0) {
    
    blue_pal <- colorRampPalette(
      c("#FFFFFF", "#0000CC")
    )(101)
    
    red_pal <- colorRampPalette(
      c("#FFFFFF", "#CC0000")
    )(101)
    
    ## Normalize absolute magnitude relative to common range
    alpha <- pmin(edges$weight / range, 1)
    
    idx <- pmax(1, round(alpha * 100) + 1)
    
    edge_cols <- ifelse(
      edges$sign > 0,
      blue_pal[idx],
      red_pal[idx]
    )
    
    igraph::E(g)$color <- edge_cols
    
    ## Scale width using the same common range
    igraph::E(g)$width <- 1 + 6 * alpha
    
    igraph::E(g)$arrow.size <- 0.2 + 2 * alpha
    
    igraph::E(g)$curved <- 0.15
  }
  
  ## ---------- Vertex attributes ----------
  
  igraph::V(g)$size <- vertex.size
  igraph::V(g)$color <- "white"
  igraph::V(g)$frame.color <- "black"
  igraph::V(g)$label.color <- "black"
  
  ## ---------- Plot ----------
  
  plot(
    g,
    layout = igraph::layout_in_circle(g),
    vertex.label.cex = vertex.label.cex
  )
}


decomp_row <- function(gcn, shared, mod1, mod2, mod3,
                       first_row = FALSE,
                       range) {
  
  ## GCN
  plot_signed_graph(
    gcn,
    labels = 1:9,
    range = range
  )
  
  if (first_row)
    title(main = "GCN",
          col.main = "blue", 
          cex.main = 2)
  
  
  ## =
  plot.new()
  text(0.5, 0.5, "=", cex = 4)
  
  
  ## Shared
  plot_signed_graph(
    shared,
    labels = 1:9,
    range = range
  )
  
  if (first_row)
    title(main = "Shared",
          col.main = "blue", 
          cex.main = 2)
  
  
  ## +
  plot.new()
  text(0.5, 0.5, "+", cex = 4)
  
  
  ## Age
  plot_signed_graph(
    mod1,
    labels = 1:9,
    range = range
  )
  
  if (first_row)
    title(main = "Age",
          col.main = "blue", 
          cex.main = 2)
  
  
  ## +
  plot.new()
  text(0.5, 0.5, "+", cex = 4)
  
  
  ## LC-Volume
  plot_signed_graph(
    mod2,
    labels = 1:9,
    range = range
  )
  
  if (first_row)
    title("LC-Volume",
          col.main = "blue", 
          cex.main = 2)
  
  
  ## +
  plot.new()
  text(0.5, 0.5, "+", cex = 4)
  
  
  ## CNR
  plot_signed_graph(
    mod3,
    labels = 1:9,
    range = range
  )
  
  if (first_row)
    title("LC-CNR",
          col.main = "blue", 
          cex.main = 2)
}


plot_decomp <- function(lcdf_norm, cv.model, subjects) {
  
  n_subjects <- length(subjects)
  
  ## ------------------------------------------------------------
  ## Calculate subject-specific GCNs and remove diagonal
  ## ------------------------------------------------------------
  
  gcns <- lapply(subjects, function(x) {
    
    A <- cv.model$bysubject_coeffs[[x]]
    
    diag(A) <- 0
    
    A
  })

  ## ------------------------------------------------------------
  ## Determine common range from GCNs
  ## ------------------------------------------------------------
  
  range <- max(
    sapply(gcns, function(A) max(abs(A), na.rm = TRUE))
  )
  
  ## ------------------------------------------------------------
  ## Layout:
  ##
  ## GCN = Shared + Age + LC-Volume + LC-CNR
  ##
  ## 9 columns:
  ## plot, =, plot, +, plot, +, plot, +, plot
  ## ------------------------------------------------------------
  
  lay <- matrix(
    seq_len(n_subjects * 9),
    nrow = n_subjects,
    byrow = TRUE
  )
  
  layout(
    lay,
    widths = c(4, 0.7, 4, 0.7, 4, 0.7, 4, 0.7, 4),
    heights = rep(1, n_subjects)
  )
  
  
  ## ------------------------------------------------------------
  ## Shared structure
  ## ------------------------------------------------------------
  
  shared <- cv.model$joint_coeffs
  diag(shared) <- 0
  
  
  ## ------------------------------------------------------------
  ## Generate rows
  ## ------------------------------------------------------------
  
  lcdf_norm <- lcdf_norm %>% as.matrix()
  
  for (k in seq_along(subjects)) {
    
    x <- subjects[k]
    
    ## Subject-specific moderator contributions
    
    age <- lcdf_norm[x, 1] *
      cv.model$moderator_coeffs[[1]]
    
    lc_volume <- lcdf_norm[x, 3] *
      cv.model$moderator_coeffs[[3]]
    
    cnr <- lcdf_norm[x, 4] *
      cv.model$moderator_coeffs[[4]]
    
    
    ## Draw row
    
    decomp_row(
      gcn = gcns[[k]],
      shared = shared,
      mod1 = age,
      mod2 = lc_volume,
      mod3 = cnr,
      first_row = (k == 1),
      range = range
    )
  }
}

lcdf_norm %>%
  select(-Gender) %>%
  as.matrix() %>%
  pairs()



which.min(lcdf_norm$Age)
which.max(lcdf_norm$CNR)
which.max(lcdf_norm$LC_vol)
which.max(lcdf_norm$LC_vol[-15]) + 1


par(mar = c(1, 0, 2, 0))
plot_decomp(
  lcdf_norm, 
  cv.model, 
  subjects = c(11,14,19))
  #subjects = c(11,14,15))

