

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
    layout = igraph::layout_in_circle(g)
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




