# may need to install these on UCR office comp
install.packages("tseries")
install.packages("MTS")
install.packages("ppcor")
install.packages("TSA")
install.packages("multivar")
install.packages("forecast")
install.packages("BigVAR")

##############
### README ###
##############
# Code currently uses the "scriptdir" folder with script helper functions to run MODVAR.
# Jose may update this.
# If he updates functions and this code breaks, the helper function scripts archived at the time 
# of writing this code (11/14/26) are in the "datadir" folder as a backup.

# These helper functions enable the full MOD-VAR analysis with:
# Shared moderated components, fully idiographic moderated components, and random components.

#change default package install directory
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))

#set wd
datadir <- "C:/Users/holly/GitHub/OSGB-Working-Group/rvar_project/real_data/ARI_MODVAR_pnomics"
scriptdir <- "C:/Users/holly/GitHub/OSGB-Working-Group/rvar_project/experiments_cluster2/041_modvar"

# which needed just for MODVAR?
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

## Model fitting packages:
library(glmnet)
library(mvtnorm)
library(multivar)
library(BigVAR)
library(expm)
library(gridExtra)

# read in data
library(openxlsx)
ari_data <- read.xlsx(paste0(datadir,"/ari_full_n54.xlsx"))
View(ari_data)

# Set NAs and compute new vars
ari_data <- ari_data |>
  mutate(across(where(is.character), \(x) na_if(x, "NA"))) |>
  rename(workingmem       = va_acc,
         attn_intensity   = m.rt,
         attn_consistency = cv.rt,
         exec_fn          = sum_flanker) |>
  mutate(f0m1      = gender - 1,
         rem_prop  = rem / sleep_time,
         deep_prop = deep / sleep_time)

ls(ari_data)

# Check gender 
ari_data %>%
  distinct(id, gender) %>%
  count(gender)

  # drop ppl with other or non-biary
  gender_drop_ids <- ari_data %>%
  filter(gender %in% c(3, 4)) %>%
  distinct(id) %>%
  pull(id)

ari_data <- ari_data %>%
  filter(!id %in% gender_drop_ids)

message("Dropped ", length(gender_drop_ids), " participant(s) with gender = 3 or 4; ",
        n_distinct(ari_data$id), " remaining")

# # standardize variables on all data
# ari_stdall <- ari_data |>
#   mutate(across(where(is.numeric) & !id & !gender & !gpa & !employed, \(x) as.numeric(scale(x))))
# head(ari_stdall)

# # standardize variables on individual data
# ari_stdid <- ari_time |>
#    group_by(id) |> 
#    mutate(across(where(is.numeric) & !gender & !gpa & !employed, \(x) as.numeric(scale(x)))) 

### Set up data for Y

# select which standardized dataset 
ari_working <- ari_data

# Code time as variable t 
ari_time <- ari_working %>%
  arrange(id, true_date) %>% # It's good practice to arrange by ID and Time first
group_by(id) %>%
  mutate(t = row_number()) %>%
  ungroup()

# subset data to contain only id, t, and Y vars
# all possible y vars here
# yvars <- c("depression", "anxiety", "distress", "somatization",
#            "agreeableness", "conscientiousness", "emotional_stability", "openness",
#            "workingmem", "attn_intensity", "attn_consistency", "exec_fn",
#            "hrv", "hrv_deep", "o2_mean", 
#            "sleep_time", "rem_prop", "deep_prop", "awake_time", "subjective_sleep")

# subset yvars with a few removed
yvars <- c("distress",
           "conscientiousness", "openness",
           "rem_prop", "deep_prop", "sleep_time",
           "o2_mean",
           "workingmem", "attn_intensity", "attn_consistency")

# check yvars for low variability
rel_sd_cut  <- 0.10   # SD < 10% of the typical within-person SD
modal_cut   <- 0.80   # >= 80% of time points at a single value
distinct_cut <- 3     # fewer than 3 distinct values

var_check <- ari_time %>%
  select(id, all_of(yvars)) %>%
  pivot_longer(all_of(yvars), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(id, variable) %>%
  summarise(
    sd         = sd(value),
    n_distinct = n_distinct(value),
    prop_modal = max(table(value)) / n(),
    .groups = "drop"
  ) %>%
  group_by(variable) %>%
  mutate(rel_sd = sd / median(sd, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    zero_sd    = is.na(sd) | sd == 0,
    low_rel_sd = !zero_sd & rel_sd < rel_sd_cut,
    sticky     = prop_modal >= modal_cut,
    few_values = n_distinct < distinct_cut,
    flag       = zero_sd | low_rel_sd | sticky | few_values
  )

# one row per participant, listing flagged variables
low_var_by_id <- var_check %>%
  filter(flag) %>%
  group_by(id) %>%
  summarise(
    n_low_var    = n(),
    n_zero_sd    = sum(zero_sd),
    zero_sd_vars = paste(variable[zero_sd], collapse = ", "),
    low_var_vars = paste(variable, collapse = ", "),
    .groups = "drop"
  )

message(nrow(low_var_by_id), " of ", n_distinct(ari_time$id),
        " participants have at least one low-variability variable; ",
        sum(low_var_by_id$n_zero_sd > 0), " have at least one zero-SD variable")
print(low_var_by_id, n = Inf)

# which variables are driving it, and by which criterion
var_check %>%
  filter(flag) %>%
  group_by(variable) %>%
  summarise(
    n_ids      = n(),
    zero_sd    = sum(zero_sd),
    low_rel_sd = sum(low_rel_sd),
    sticky     = sum(sticky),
    few_values = sum(few_values),
    .groups = "drop"
  ) %>%
  arrange(desc(n_ids))

# check how many participants have 25 or more complete rows on all variables
min_complete <- 25
low_n_by_id <- ari_time %>%
  group_by(id) %>%
  summarise(
    n_joint_complete = sum(if_all(all_of(yvars), ~ !is.na(.x))),
    .groups = "drop"
  ) %>%
  filter(n_joint_complete < min_complete)
message("Dropping ", nrow(low_n_by_id), " of ", n_distinct(ari_time$id),
        " participants with < ", min_complete, " joint-complete rows, ", (n_distinct(ari_time$id)-nrow(low_n_by_id)), " participants remaining")
print(low_n_by_id, n = Inf)

# create filtered dataset with only participants with min complete time points 
ari_time_filtered <- ari_time %>%
  filter(!id %in% low_n_by_id$id,
         !id %in% low_var_by_id$id)

# see how many ppl are left
n_distinct(ari_time_filtered$id)

# Select variables for analysis
ari_y <- ari_time_filtered %>%
  select(id, t, yvars)
ls(ari_y)

# Select complete cases only
ari_y_complete <- ari_y[complete.cases(ari_y), ]
View(ari_y_complete)
# Check there are no infinite values
stopifnot(sum(!is.finite(as.matrix(ari_y_complete[yvars]))) == 0)

# create unique ID for complete data
id_vec <- unique(ari_y_complete$id)

# check length of id_vec
length(unique(id_vec))

# create list with each subject's data as dataframes in list
ari_y_list <- lapply(
  id_vec,
  function(data, subj) {
    data %>%
      dplyr::filter(id == subj) %>%
      dplyr::arrange(t) %>%
      dplyr::select(-t,-id)
  },
  data = ari_y_complete
)

# check length of each complete dataset
map_int(ari_y_list, nrow)

stopifnot(all(map_int(ari_y_list, nrow) >= min_complete))

##################################################
##################################################
##################################################
# MOD-VAR (MO / IO / MI) on ARI data, scored against a persistence forecast
# Requires ari_y_list, id_vec, and ari_data from the prep script, saved with:
#   save(ari_y_list, id_vec, ari_data, file = "data/processed/ari_prepped.RData")

library(expm)          # %^%
library(Rcpp)
library(plot.matrix)   # plot() method for matrices

helper_files <- c("matrix_fista.cpp", "auxfunct.r", "adaweights.r",
                  "bic.modvar.r", "cv.modvar.r", "ada.modvar.r")
missing <- helper_files[!file.exists(helper_files)]
if (length(missing) > 0) stop("Not in ", getwd(), ": ", toString(missing))

Rcpp::sourceCpp("matrix_fista.cpp")
walk(helper_files[-1], source)

dir.create("output", showWarnings = FALSE)

# ---- Settings ------------------------------------------------------------
n_test   <- 6                                  # held-out points per person
horizon  <- 4                                   # forecast steps; must be < n_test
lambdas1 <- 10^seq(0, -5, length.out = 10)
ratios   <- 10^seq(2, -2, length.out = 10)

strip <- \(m) matrix(as.numeric(m), nrow = nrow(m))   # drop dimnames + scale attrs

stopifnot(length(ari_y_list) == length(id_vec), n_test > horizon)

# ---- Train/test split, scaled on training stats only ---------------------
split_list <- map(ari_y_list, \(d) {
  d       <- as.matrix(d)
  n_train <- nrow(d) - n_test
  train   <- d[seq_len(n_train), , drop = FALSE]
  test    <- d[(n_train + 1):nrow(d), , drop = FALSE]
  mu <- colMeans(train)
  s  <- apply(train, 2, sd)
  s_test  <- apply(test, 2, sd)
  list(train          = scale(train, mu, s),
       test           = scale(test, mu, s),
       zero_sd        = any(s == 0 | is.na(s)),
       zero_sd_test   = any(s_test == 0 | is.na(s_test)),
       zero_vars      = colnames(d)[s == 0 | is.na(s)],
       zero_vars_test = colnames(d)[s_test == 0 | is.na(s_test)])
})

# ---- Warn about zero training / test variance ----------------------------
zero_var_tbl <- tibble(
  id             = id_vec,
  zero_train     = map_lgl(split_list, "zero_sd"),
  zero_test      = map_lgl(split_list, "zero_sd_test"),
  zero_vars      = map_chr(split_list, \(x) paste(x$zero_vars, collapse = ", ")),
  zero_vars_test = map_chr(split_list, \(x) paste(x$zero_vars_test, collapse = ", "))
) |>
  filter(zero_train | zero_test)

if (nrow(zero_var_tbl) > 0) {
  warning(sum(zero_var_tbl$zero_train), " participant(s) with zero training variance; ",
          sum(zero_var_tbl$zero_test), " with zero test variance; ",
          sum(zero_var_tbl$zero_train & zero_var_tbl$zero_test), " with both (of ",
          length(id_vec), " total)")
  print(zero_var_tbl, n = Inf)
}

# use if the code directly above gives a warning about zero variance in test or training
if (any(map_lgl(split_list, "zero_sd"))) {
  keep       <- !map_lgl(split_list, "zero_sd")
  split_list <- split_list[keep]
  ari_y_list <- ari_y_list[keep]
  id_vec     <- id_vec[keep]
  message(length(id_vec), " participants retained after dropping zero training variance")
}

# convert lists to matrices
dfn_list   <- map(split_list, \(x) strip(x$train))
Y_forecast <- map(split_list, \(x) strip(x$test))

# ---- Between-person moderators, row-aligned with id_vec ------------------
# TODO remove any gender not 0/1
mod_df <- ari_data |>
  filter(id %in% id_vec) |>
  summarise(f0m1     = mean(f0m1, na.rm = TRUE),
            employed = mean(employed, na.rm = TRUE),
            .by = id) |>
  slice(match(id_vec, id))

# Check employment and gender frequencies
  table(mod_df$f0m1)
  table(mod_df$employed)

stopifnot(all(mod_df$id == id_vec), !anyNA(mod_df))

lcdf_norm <- mod_df |> select(-id) |> as.matrix() |> scale() |> strip()
stopifnot(!anyNA(lcdf_norm), nrow(lcdf_norm) == length(dfn_list))

############################################################################
# ---- Fit models ----------------------------------------------------------
############################################################################
dfn_mats <- unname(dfn_list)

MO.model <- cv.modvar(dfn_mats, X = lcdf_norm,
                      lambdas1 = lambdas1, ratios = ratios,
                      multi = FALSE, cv.type = "rolling")
saveRDS(MO.model, "output/ari_MO_model.rds")

IO.model <- cv.modvar(dfn_mats, X = NULL,
                      lambdas1 = lambdas1, ratios = ratios,
                      multi = TRUE, cv.type = "rolling")
saveRDS(IO.model, "output/ari_IO_model.rds")

MI.model <- cv.modvar(dfn_mats, X = lcdf_norm,
                      lambdas1 = lambdas1, ratios = ratios,
                      multi = TRUE, cv.type = "rolling")
saveRDS(MI.model, "output/ari_MI_model.rds")

models <- list(
  `MOD-VAR: MOD`    = MO.model,
  `MOD-VAR: FI`     = IO.model,
  `MOD-VAR: MOD+FI` = MI.model
)

# ---- Reload saved models if not re-running -------------------------------------------------
MO.model <- readRDS("output/ari_MO_model.rds")
IO.model <- readRDS("output/ari_IO_model.rds")
MI.model <- readRDS("output/ari_MI_model.rds")

models <- list(
  `MOD-VAR: MOD`    = MO.model,
  `MOD-VAR: FI`     = IO.model,
  `MOD-VAR: MOD+FI` = MI.model
)

# 1. Does the individual-only model (IO) estimate individual effects?
summary(map_int(IO.model$idiographic_coeffs, \(A) sum(abs(A) > 1e-8)))

# 2. Where did CV land for MI?
MI.model$lambda1
MI.model$ratios
which(MI.model$eval.mat == min(MI.model$eval.mat), arr.ind = TRUE)
dim(MI.model$eval.mat)

# =========================================================================
# DECOMPOSITION FIGURE: GCN = Shared + Gender + Employed + Individual
# =========================================================================
# Requires: plot_signed_graph() (second version, with `range`), MI.model,
# lcdf_norm, mod_df, id_vec, yvars
source("plot_helpers.R")

zero_diag <- \(A) { diag(A) <- 0; A }

# node labels (short) and domains, in yvars order; edit to match your yvars
node_info <- tibble::tribble(
  ~var,                ~short,        ~domain,
  "distress",          "Dist",        "Mood/personality",
  "conscientiousness", "Cons",        "Mood/personality",
  "openness",          "Open",        "Mood/personality",
  "rem_prop",          "REM",         "Sleep",
  "deep_prop",         "Deep",        "Sleep",
  "sleep_time",        "Slp",         "Sleep",
  "o2_mean",           "O\u2082",     "Physiology",
  "workingmem",        "WM",          "Cognition",
  "attn_intensity",    "AtI",         "Cognition",
  "attn_consistency",  "AtC",         "Cognition"
) |>
  slice(match(yvars, var))
stopifnot(nrow(node_info) == length(yvars), !anyNA(node_info$short))

# separate plots by gender
mod_names <- c("Gender", "Employed")
subj_F <- match(c(7007, 7010), id_vec)
subj_M <- match(c(7014, 7041), id_vec)
stopifnot(!anyNA(c(subj_F, subj_M)))

# one color scale shared by both slides
rng <- slide_range(MI.model, lcdf_norm, mod_names, c(subj_F, subj_M))

save_decomp_slide("output/ari_MI_slide_F.png", MI.model, lcdf_norm, mod_names,
                  subjects = subj_F, node_labels = node_labels,
                  row_labels = c("Unemployed", "Employed"), range = rng,
                  key_cex = 1.9
                  )

save_decomp_slide("output/ari_MI_slide_M.png", MI.model, lcdf_norm, mod_names,
                  subjects = subj_M, node_labels = node_labels,
                  row_labels = c("Unemployed", "Employed"), range = rng,
                  key_cex = 1.9
                  )

# reconstruction check
recon_ok <- check_decomp(MI.model, lcdf_norm)
table(recon_ok)

# node key
node_key <- tibble(node = seq_along(yvars), variable = yvars)
print(node_key, n = Inf)

# pick subjects (positions in id_vec)
mod_names <- c("Gender", "Employed")           # order = columns of lcdf_norm
print(mod_df)
# unemployed F, employed F, unemployed M, employed M
subjects  <- match(c(7007, 7010, 7014, 7041), id_vec)
stopifnot(!anyNA(subjects))

# draw
preview_decomp(MI.model, lcdf_norm, mod_names, subjects,
               n_subjects = length(subjects))

# save
save_decomp_png("output/ari_MI_decomposition.png",
                MI.model, lcdf_norm, mod_names, subjects,
                n_subjects = length(subjects))

### broken out by slide, new figs --------------------------------------------
##############################################################################
domain_cols   <- c("Mood/personality" = "#f8f1d7", "Sleep" = "#d5e7ee",
                   "Physiology" = "#efe6ef", "Cognition" = "#d6dfd5")
vertex_colors <- unname(domain_cols[node_info$domain])

# moderator effects as group differences (binary moderators, z-scored)
gender_diff <- MI.model$moderator_coeffs[[1]] * diff(range(lcdf_norm[, 1]))
employ_diff <- MI.model$moderator_coeffs[[2]] * diff(range(lcdf_norm[, 2]))

comp_mats <- lapply(list(MI.model$joint_coeffs, gender_diff, employ_diff), zero_diag)
gcn_ids   <- c(7007, 7010, 7014, 7041)   # F-unemp, F-emp, M-unemp, M-emp
gcn_mats  <- lapply(match(gcn_ids, id_vec), \(x) zero_diag(MI.model$bysubject_coeffs[[x]]))

# one color scale for both slides
rng <- max(sapply(c(comp_mats, gcn_mats), \(A) max(abs(A))))

# Slide 1: components
save_network_grid("output/slide_components.png",
                  mats = comp_mats,
                  titles = c("Shared", "Men vs. women", "Employed vs. not"),
                  nrow = 1, ncol = 3,
                  labels = node_info$short, vertex_colors = vertex_colors,
                  legend_domains = domain_cols, range = rng,
                  panel_title_cex = 3,     # was 2.2
                  legend_cex = 2.4,        # was 1.8
                  legend_width = 3.8)      # wider column for larger key text)

# Slide 2: resulting networks, 2 x 2
save_network_grid("output/slide_gcn_grid.png",
                  mats = gcn_mats,
                  titles = c("Woman, not employed", "Woman, employed",
                             "Man, not employed", "Man, employed"),
                  nrow = 2, ncol = 2,
                  labels = node_info$short, vertex_colors = vertex_colors,
                  legend_domains = domain_cols, range = rng,
                  panel_title_cex = 3,     # was 2.2
                  legend_cex = 2.4,        # was 1.8
                  legend_width = 3.8)      # wider column for larger key text)

############################################################################
# ---- Tuning Diagnostics ----------------------------------------------------------
############################################################################

# ---- Tuning diagnostics: is the CV optimum on the grid edge? -------------
iwalk(models, \(m, nm) {
  idx <- which(m$eval.mat == min(m$eval.mat), arr.ind = TRUE)[1, ]
  on_edge <- any(idx == 1) || idx[1] == nrow(m$eval.mat) || idx[2] == ncol(m$eval.mat)
  message(nm, ": optimum at [", idx[1], ", ", idx[2], "]",
          if (on_edge) "  <-- ON GRID EDGE, widen lambdas1/ratios" else "")
})

par(mfrow = c(1, 3))
iwalk(models, \(m, nm) plot(log10(m$eval.mat), main = nm, breaks = 20))

# ---- Forecast evaluation -------------------------------------------------
eval_forecast <- function(Y_forecast, B_est_list, horizon) {
  n    <- length(Y_forecast)
  msfe <- matrix(NA_real_, nrow = horizon, ncol = n)
  for (h in seq_len(horizon)) {
    for (k in seq_len(n)) {
      Y <- Y_forecast[[k]]
      r <- nrow(Y)
      if (r <= h) next
      B <- B_est_list[[k]] %^% h
      Z <- Y[(h + 1):r, , drop = FALSE]
      X <- Y[1:(r - h), , drop = FALSE]
      msfe[h, k] <- sqrt(mean((Z - X %*% t(B))^2))
    }
  }
  rownames(msfe) <- paste0("msfe_step", seq_len(horizon))
  msfe
}

B_persist    <- map(dfn_list, \(x) diag(ncol(x)))
msfe_persist <- eval_forecast(Y_forecast, B_persist, horizon)

rel_error <- \(B_list) (eval_forecast(Y_forecast, B_list, horizon) - msfe_persist) / msfe_persist

perf_data <- imap(models, \(m, nm) {
  rel_error(m$bysubject_coeffs) |>
    t() |>
    as_tibble() |>
    mutate(id = id_vec) |>
    pivot_longer(starts_with("msfe_step"), names_to = "step",
                 names_prefix = "msfe_step", names_transform = as.integer,
                 values_to = "rel_error") |>
    mutate(type = nm)
}) |>
  list_rbind()

perf_summary <- perf_data |>
  summarise(median_rel_error = median(rel_error, na.rm = TRUE),
            prop_beats_persistence = mean(rel_error < 0, na.rm = TRUE),
            n = sum(!is.na(rel_error)),
            .by = c(type, step))
write_csv(perf_summary, "output/ari_modvar_perf_summary.csv")
print(perf_summary, n = Inf)

p <- ggplot(perf_data, aes(x = factor(step), y = rel_error, fill = type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
  coord_cartesian(ylim = c(-0.5, 1)) +   # zooms without dropping data (ylim() drops)
  labs(
    title    = "MOD-VAR forecast error relative to persistence",
    subtitle = "Below 0 = more accurate than carrying the last value forward",
    x        = "Forecast horizon (time points ahead)",
    y        = "Relative RMSE vs. persistence (proportion)",
    fill     = "Model"
  ) +
  theme_minimal()

ggsave("output/ari_modvar_rel_error.png", p, width = 8, height = 5, dpi = 300)
p



























eval_forecast <- function(Y_forecast, B_est_list, range, horizon) {
  
  n <- length(Y_forecast)
  
  #msfe <- numeric(horizon)
  msfe <- matrix(0, ncol = n, nrow = horizon)
  for (h in 1:horizon) {
    
    bysubj_msfe <- numeric(n)
    for (k in 1:n) {
      Y <- Y_forecast[[k]]
      B <- B_est_list[[k]] %^% h
      
      Z <- Y[h + 1:(range - h),]
      X <- Y[1:(range - h),]
      
      err_mat <- Z - X %*% t(B)
      
      bysubj_msfe[k] <- sqrt(mean(err_mat^2))
      
    }
    msfe[h,] <- bysubj_msfe
    
  }
  
  rownames(msfe) <- paste0("msfe_step", 1:h)
  
  return(msfe)
  
}
##################################################
##################################################
##################################################

Tsamp <- 100
dfn_list <- df_list %>% 
  lapply(function(x) {
    x %>%
      dplyr::select(-t,-sub) %>%
      mutate_if(is.numeric, scale) %>%
      return()
  })
dfn_list[[1]] %>% apply(2,mean)  
dfn_list[[1]] %>% apply(2,var)  



models <- list()
results <- list()
Y_forecast <- list()
for(i in 1:length(df_list)) {
  models[[i]] <- df_list[[i]][1:Tsamp,] %>% 
    tibble() %>%
    dplyr::select(-sub,-t) %>% 
    as.matrix() %>%
    BigVAR::constructModel( 
      p = 1,
      gran = c(50,10),
      struct = "Basic",
      cv = "Rolling",
      verbose = TRUE,
      ownlambdas = FALSE,
      model.controls=list(intercept=FALSE),
      linear = FALSE)
  results[[i]] <- BigVAR::cv.BigVAR(models[[i]])
  print(coef(results[[i]]))
  
  Y_forecast[[i]] <- df_list[[i]][-(1:Tsamp),] %>% 
    tibble() %>%
    dplyr::select(-sub,-t) %>% 
    as.matrix() 

}

results[[1]]@Loss

str(results[[1]])
apply(results[[1]]@InSampMSFE, 2, mean)
apply(results[[1]]@OOSMSFE, 2, mean)
results[[1]]
results[[1]]@LambdaGrid
results[[1]]@InSampMSFE

B_est_list <- lapply(results,coef) %>% 
  lapply(as.matrix) 
eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) %>%
  apply(1, mean)
##RW Out of Sample Loss
## [1] 0.294


B_est_list2 <- rep(ncol(df) - 2, length(df_list)) %>%
  lapply(function(x) {return(diag(x))})
results %>% lapply(function(x) {})



par(mar = c(5.1, 4.1, 4.1, 4.1))
((eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  plot(breaks = 10)


((eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  apply(1, mean) %>%
  plot()


## Benefits compared to persistent forecasting.
((eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
    geom_boxplot() +
    geom_hline(yintercept = 0)
  
## We observe that for the 1st step prediction, the model performs
## somewhat poorly (~1% worse than persistence forecasting), but
## as we increase the time horizon, the estimation improves
## significantly (on average, 20% relative improvement)


##################################################
##################################################
##################################################
##################################################
##################################################
##################################################

dfn_list <- df_list %>% 
  lapply(function(x) {
    x %>%
      dplyr::filter(t %in% 1:Tsamp) %>%
      dplyr::select(-t,-sub) %>%
      mutate_if(is.numeric, scale) %>%
      return()
  })

dfn_list[[1]] %>% dim()
dfn_list[[1]] %>% apply(2,mean)  
dfn_list[[1]] %>% apply(2,var)  


modelM50 <- multivar::constructModel(data = dfn_list, lassotype = "standard")
fitM50 <- multivar::cv.multivar(modelM50)
fitM50$mats %>% str()
B_estM50_list <- (fitM50$mats)$total
str(fitM50)

dim(fitM50$MSFE)
str(fitM50)

((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  plot(breaks = 10)


((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  apply(1, mean) %>%
  plot()



## Benefits compared to persistent forecasting.
rfeM50 <- ((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
p1 <- rfeM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  ggtitle("Multi VAR")

## Benefits compared to persistent forecasting.
rfeVAR <- (eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)
p2 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  ggtitle("VAR")



grid.arrange(p1, p2, ncol = 2)



##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
## MOD-VAR with only individual structure:
Rcpp::sourceCpp("041_modvar/matrix_fista.cpp")
source("041_modvar/auxfunct.r")
source("041_modvar/adaweights.r")
source("041_modvar/bic.modvar.r")
source("041_modvar/cv.modvar.r")
source("041_modvar/ada.modvar.r")


#lambdas1  <- 10^(seq(3, 0, length.out = 10)) ## Surprisingly good!
# ratios    <- 10^(seq(2, -2, length.out = 10))

# lambdas1  <- 10^(seq(2, -1, length.out = 10)) # Improved early prediction,
# ratios    <- 10^(seq(2, -2, length.out = 10)) # worsen later prediction. 

# lambdas1  <- 10^(seq(2, -1, length.out = 10)) # Best so far!
# ratios    <- 10^(seq(2, 0, length.out = 10)) 

# lambdas1  <- 10^(seq(3, -3, length.out = 50)) # Bad! WTF
# ratios    <- 10^(seq(3, -3, length.out = 50)) 



lambdas1  <- 10^(seq(2, -2, length.out = 30))
ratios    <- 10^(seq(2, 0, length.out = 30))
cv.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = NULL,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = TRUE,
    cv.type = "rolling")


save.image("BrainNetworksGeneData/102_Data.RData")
rm(list = ls())
##################################################
##################################################
load("BrainNetworksGeneData/102_Data.RData")


par(mfrow = c(1,2),
    mar = c(5.1, 4.1, 4.1, 4.1))
plot(log(cv.model$eval.mat, 10), breaks = 30)
plot(cv.model$eval.mat == min(cv.model$eval.mat))


cv.model$eval.mat
cv.model$lambda1
cv.model$ratios

par(mfrow = c(2,2))

cv.model$bysubject_coeffs[1:4] %>%
  lapply(plot)
(fitM50$mats)$total[1:4] %>%
  lapply(plot)
B_est_list [1:4] %>%
  lapply(plot)

cv.model$moderator_coeffs
cv.model$bysubject_coeffs[[1]]
cv.model$joint_coeffs
cv.model$idiographic_coeffs[[1]]


cv.model$bysubject_coeffs[[1]] == cv.model$idiographic_coeffs[[1]]


B_estMOD50_list <- cv.model$bysubject_coeffs

str(B_estMOD50_list)

((eval_forecast(Y_forecast, B_estMOD50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  plot(breaks = 10)


((eval_forecast(Y_forecast, B_estMOD50_list, range = 100, horizon = 10) -
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
    eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) %>% 
  apply(1, mean) %>%
  plot()




## Benefits compared to persistent forecasting.
rfeMOD50 <- ((eval_forecast(Y_forecast, B_estMOD50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
p1 <- rfeMOD50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  ylim(-0.5,1) +
  ggtitle("MOD-VAR")

## Benefits compared to persistent forecasting.
rfeM50 <- ((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
p2 <- rfeM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  ylim(-0.5,1) +
  ggtitle("Multi VAR")

## Benefits compared to persistent forecasting.
rfeVAR <- (eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)
p3 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE)) %>%
  ggplot(aes(x = StepNo, y = error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  ylim(-0.5,1) +
  ggtitle("VAR")


grid.arrange(p1, p2, p3, ncol = 3)



##################################################
##################################################
dMOD50 <- rfeMOD50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: FI")

dM50 <- rfeM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "Multi-VAR")

dV50 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "VAR")

perf_data <- rbind(dMOD50, dM50, dV50) %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE))

head(perf_data)

ggplot(perf_data, aes(x = factor(StepNo), y = error, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    x = "Step Number",
    y = "Error",
    fill = "Type"
  ) +
  theme_minimal()








##################################################
##################################################
##################################################
##################################################
##################################################
##################################################

## ALL OPTIONS OF MOD-VAR MODELS:
Rcpp::sourceCpp(paste0(scriptdir,"/matrix_fista.cpp"))
source(paste0(scriptdir,"/auxfunct.r"))
source(paste0(scriptdir,"/adaweights.r"))
source(paste0(scriptdir,"/bic.modvar.r"))
source(paste0(scriptdir,"/cv.modvar.r"))
source(paste0(scriptdir,"/ada.modvar.r"))

#lambdas1  <- 10^(seq(2, -2, length.out = 10)) # Best so far!
#ratios    <- 10^(seq(2, 0, length.out = 10)) 

lambdas1  <- 10^(seq(0, -5, length.out = 10)) # Best so far!
ratios    <- 10^(seq(2, -2, length.out = 10)) 

lcdf_norm <- lcdf_norm %>% 
  as.matrix()
colnames(lcdf_norm) <- NULL 
MO.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = lcdf_norm,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = FALSE,
    cv.type = "rolling")


IO.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = NULL,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = TRUE,
    cv.type = "rolling")


MI.model <- dfn_list %>% 
  lapply(function(x) {
    x <- as.matrix(x)
    colnames(x) = NULL
    return(x)
  }) %>%
  cv.modvar(
    X = lcdf_norm,
    lambdas1 = lambdas1,
    ratios = ratios,
    multi = TRUE,
    cv.type = "rolling")


save.image("BrainNetworksGeneData/102_Data.RData")
rm(list = ls())
##################################################
##################################################
## Visualizations before comparison:
load("BrainNetworksGeneData/102_Data.RData")

## Verifying correct fit conditions:
str(MO.model$idiographic_coeffs)
str(MO.model$moderator_coeffs)
str(MO.model$joint_coeffs)

str(IO.model$idiographic_coeffs)
str(IO.model$moderator_coeffs)
str(IO.model$joint_coeffs)

str(MI.model$idiographic_coeffs)
str(MI.model$moderator_coeffs)
str(MI.model$joint_coeffs)

## Are the models having the same sparsity pattern?
par(mfcol = c(3,2))
plot(log(abs(MI.model$joint_coeffs) + 1e-10, 10) > -2)
plot(log(abs(MO.model$joint_coeffs) + 1e-10, 10) > -2)
plot(log(abs(IO.model$joint_coeffs) + 1e-10, 10) > -2)

plot(MI.model$joint_coeffs == MO.model$joint_coeffs)
plot(MI.model$joint_coeffs == IO.model$joint_coeffs)


## How does the evaluation matrix look?
par(mfcol = c(3,2))
plot(log(MO.model$eval.mat, base = 10))
plot(log(IO.model$eval.mat, base = 10))
plot(log(MI.model$eval.mat, base = 10))

plot(MO.model$eval.mat == min(MO.model$eval.mat))
plot(IO.model$eval.mat == min(IO.model$eval.mat))
plot(MI.model$eval.mat == min(MI.model$eval.mat))



##################################################
##################################################
## Visualizing results:
B_estMO50_list <- MO.model$bysubject_coeffs
B_estIO50_list <- IO.model$bysubject_coeffs
B_estMI50_list <- MI.model$bysubject_coeffs

rfeMO50 <- ((eval_forecast(Y_forecast, B_estMO50_list, range = 100, horizon = 10) -
               eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeIO50 <- ((eval_forecast(Y_forecast, B_estIO50_list, range = 100, horizon = 10) -
               eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 
rfeMI50 <- ((eval_forecast(Y_forecast, B_estMI50_list, range = 100, horizon = 10) -
               eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 


rfeM50 <- ((eval_forecast(Y_forecast, B_estM50_list, range = 100, horizon = 10) -
              eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) 

rfeVAR <- (eval_forecast(Y_forecast, B_est_list, range = 100, horizon = 10) -
             eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)) /
  eval_forecast(Y_forecast, B_est_list2, range = 100, horizon = 10)

#######################
#######################
## Generating plots:
dMO50 <- rfeMO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: MOD")

dIO50 <- rfeMO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: FI")

dMI50 <- rfeMO50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "MOD-VAR: MOD+FI")

dM50 <- rfeM50 %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "Multi-VAR")

dV50 <- rfeVAR %>% 
  t() %>% as.data.frame() %>%
  tibble() %>% 
  pivot_longer(cols = msfe_step1:msfe_step10,
               names_to = "StepNo",
               names_prefix = "msfe_step",
               values_to = "error") %>%
  mutate(type = "VAR")

perf_data <- rbind(dMO50, dIO50, dMI50,  dM50, dV50) %>%
  mutate(StepNo = factor(StepNo, levels = 1:10, ordered = TRUE))

head(perf_data)

ggplot(perf_data, aes(x = factor(StepNo), y = error, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    x = "Step Number",
    y = "Error",
    fill = "Type"
  ) +
  theme_minimal()


