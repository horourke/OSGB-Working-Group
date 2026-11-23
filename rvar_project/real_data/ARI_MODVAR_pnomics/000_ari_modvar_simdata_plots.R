##########################################################################
# Simulated MOD-VAR example for plots: 5 variables, 2 people, 1 moderator
##########################################################################

library(tidyverse)
library(Rcpp)

set.seed(2026)

#  MOD-VAR helpers + plot helper
helper_files <- c("matrix_fista.cpp", "auxfunct.r", "adaweights.r",
                  "bic.modvar.r", "cv.modvar.r", "ada.modvar.r")
missing <- c(helper_files, "plot_helpers.R")[!file.exists(c(helper_files, "plot_helpers.R"))]
if (length(missing) > 0) stop("Not in ", getwd(), ": ", toString(missing))

Rcpp::sourceCpp("matrix_fista.cpp")
walk(helper_files[-1], source)
source("plot_helpers.R")
dir.create("output", showWarnings = FALSE)

# settings 
yvars    <- c("Y1", "Y2", "Y3", "Y4", "Y5")
p        <- length(yvars)
n_subj   <- 2
T_len    <- 100          # time points per person
n_test   <- 6
lambdas1 <- 10^seq(0, -5, length.out = 10)
ratios   <- 10^seq(2, -2, length.out = 10)

strip <- \(m) matrix(as.numeric(m), nrow = nrow(m))

# true coefficient matrices (row = outcome at t, col = predictor at t-1)
B_shared <- matrix(c(
  0.35,  0.00,  0.20,  0.00,  0.00,
  0.20,  0.30,  0.00,  0.00, -0.15,
  0.00, -0.25,  0.30,  0.00,  0.00,
  0.00,  0.00,  0.20,  0.35,  0.00,
  0.15,  0.00,  0.00, -0.20,  0.30
), p, p, byrow = TRUE)

B_mod <- matrix(c(
  0.00,  0.00,  0.00,  0.00,  0.00,
  0.00,  0.00,  0.20,  0.00,  0.00,
 -0.25,  0.00,  0.00,  0.00,  0.00,
  0.00,  0.00,  0.00,  0.00,  0.00,
  0.00,  0.00,  0.00,  0.20,  0.00
), p, p, byrow = TRUE)

B_idio <- list(
  { A <- matrix(0, p, p); A[1, 2] <- 0.15; A },   # person 1: Y2 -> Y1
  { A <- matrix(0, p, p); A[4, 5] <- 0.15; A }    # person 2: Y5 -> Y4
)

# moderator: person 1 = 0, person 2 = 1
X_raw <- c(0, 1)

B_true <- map(seq_len(n_subj), \(i) B_shared + X_raw[i] * B_mod + B_idio[[i]])

# stationarity check: all eigenvalues inside the unit circle
stopifnot(all(map_dbl(B_true, \(B) max(Mod(eigen(B)$values))) < 1))

# simulate t-1 lag data ----------------------------------------------------
############################################################################
# burn-in period t = 100 for time series
sim_var1 <- function(B, T_len, burn = 100, sd = 1) {
  p <- nrow(B)
  y <- matrix(0, T_len + burn, p)
  for (t in 2:(T_len + burn)) {
    y[t, ] <- B %*% y[t - 1, ] + rnorm(p, 0, sd)
  }
  y[(burn + 1):(T_len + burn), ]
}

sim_df <- map(seq_len(n_subj), \(i) {
  y <- sim_var1(B_true[[i]], T_len)
  colnames(y) <- yvars
  as_tibble(y) |> mutate(id = 900 + i, t = row_number(), .before = 1)
}) |>
  list_rbind()

mod_raw <- tibble(id = 900 + seq_len(n_subj), X = X_raw)

# give unique ID
id_vec <- unique(sim_df$id)

# create list w/ dataset for each ID
sim_y_list <- lapply(
  id_vec,
  function(data, subj) {
    data %>%
      dplyr::filter(id == subj) %>%
      dplyr::arrange(t) %>%
      dplyr::select(-t, -id)
  },
  data = sim_df
)

split_list <- map(sim_y_list, \(d) {
  d       <- as.matrix(d)
  n_train <- nrow(d) - n_test
  train   <- d[seq_len(n_train), , drop = FALSE]
  test    <- d[(n_train + 1):nrow(d), , drop = FALSE]
  mu <- colMeans(train)
  s  <- apply(train, 2, sd)
  list(train = scale(train, mu, s), test = scale(test, mu, s))
})

dfn_list   <- map(split_list, \(x) strip(x$train))
Y_forecast <- map(split_list, \(x) strip(x$test))

mod_df <- mod_raw |> slice(match(id_vec, id))
stopifnot(all(mod_df$id == id_vec))
lcdf_norm <- mod_df |> select(-id) |> as.matrix() |> scale() |> strip()

# ---- Fit MI model --------------------------------------------------------
MI.sim <- cv.modvar(unname(dfn_list), X = lcdf_norm,
                    lambdas1 = lambdas1, ratios = ratios,
                    multi = TRUE, cv.type = "rolling")
saveRDS(MI.sim, "output/sim_MI_model.rds")

table(check_decomp(MI.sim, lcdf_norm))

### matrices to plot
zero_diag <- \(A) { diag(A) <- 0; A }

# moderator effect as group difference (X = 1 vs X = 0)
mod_diff_est <- MI.sim$moderator_coeffs[[1]] * diff(range(lcdf_norm[, 1]))

est_comp  <- lapply(list(MI.sim$joint_coeffs, mod_diff_est), zero_diag)
true_comp <- lapply(list(B_shared, B_mod), zero_diag)

est_gcn  <- lapply(MI.sim$bysubject_coeffs, zero_diag)
true_gcn <- lapply(B_true, zero_diag)

vertex_colors <- c("#00000000", "#00000000", "#00000000", "#00000000", "#00000000")
rng <- max(sapply(c(est_comp, true_comp, est_gcn, true_gcn), \(A) max(abs(A))), 1e-8)

# components, true (top) vs estimated (bottom) 
save_network_grid("output/sim_components.png",
                  mats   = c(true_comp, est_comp),
                  titles = c("True shared", "True X effect",
                             "Estimated shared", "Estimated X effect"),
                  nrow = 2, ncol = 2,
                  labels = yvars, vertex_colors = vertex_colors,
                  range = rng, main = "Simulated MI model: components",
                  panel_title_cex = 2.6, legend_cex = 2.2, legend_width = 3)

# person networks, true (top) vs estimated (bottom)
save_network_grid("output/sim_gcn.png",
                  mats   = c(true_gcn, est_gcn),
                  titles = c("True: X = 0", "True: X = 1",
                             "Estimated: X = 0", "Estimated: X = 1"),
                  nrow = 2, ncol = 2,
                  labels = yvars, vertex_colors = vertex_colors,
                  range = rng, main = "Simulated MI model: person networks",
                  panel_title_cex = 2.6, legend_cex = 2.2, legend_width = 3)