############################################
############################################
## SETTING PARAMETERS:
wd <- getwd()
.libPaths(paste0(wd,"/req_lib"))

library(readr)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(stringr)
library(lazyeval)
library(pROC)

input <- commandArgs(trailingOnly = TRUE)
version <- as.numeric(input[1])
cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")

###################### Parameter table:
runtype       <- 2 # FOR EXPERIMENT RUNS
#runtype       <- 3 # FOR FULL RUNS
index_old     <- 1 # run index to use
sim_par_table <- expand.grid(
  running_days  = 2,
  range        = 50,
  horizon      = 3,
  entry_min     = 0.3,                            ## entry_min : minimum B0,B1,...,Bp entry magnitude
  entry_max     = 0.9,                            ## entry_max : maximum B0,B1,...,Bp entry magnitude
    
  g_sd          = 0.3,                            ## g_sd  : 
  u_min         = 0.5,                            ## u_min : minimum entry of exogenous X for type = "unif".
  u_max         = 1,                              ## u_max : maximum entry of exogenous X for type = "unif".
  type          = "unif",                         ## type  : distribution of exogenous variables.
  nz_x_prob     = c(0.75),                        ## nz_x_prob : proportion of entries in exogenous
                                                    ##              data matrix X with non-zero values.
  signed        = c(TRUE),                        ## signed    : are entries signed or all positive?


  prob_c        = c(0.25, 0.75),                  ## prob_c     : proportion of common entries.
  prob_tot      = 0.1,                            ## prob_tot   : total proportion of non-zero entries.
  prob_delta    = 0,                              ## prob_delta : total proportion of idiographic entries.

  nsim          = ifelse(
                      runtype == 1, 
                      1, 
                      ifelse(runtype == 2, 2, 5)),  ## nsim     : no of simulation repetitions.
  sigma2        = c(0.1),                         ## sigma2   : variance o VAR error term.
  n             = c(15, 20, 25),                  ## n        : No. of individuals
  T             = c(30, 60, 90),                  ## T        : timepoints per individual.
  p             = c(2, 5),                        ## p        : covariate dimension
  d             = c(10, 20, 30))                  ## d        : Time series dimension
## Add specific parameters to table:
sim_par_table$prob_phi0 <- sim_par_table$prob_tot * sim_par_table$prob_c
sim_par_table$prob_phip <- sim_par_table$prob_tot * (1 - sim_par_table$prob_c)
  
sim_par_table$phi0_min <- sim_par_table$entry_min
sim_par_table$phi0_max <- sim_par_table$entry_max
sim_par_table$phip_min <- sim_par_table$entry_min
sim_par_table$phip_max <- sim_par_table$entry_max
  
sim_par_table$vmax <- sim_par_table$sigma2
sim_par_table$vmin <- 0.5 * sim_par_table$sigma2
attach(sim_par_table)



###################### Creating folders:
subfolder_new        <- paste0("500_AggregatedDataExperiments/")
subfolder_data_new   <- paste0(subfolder_new, "data_all/")
subfolder_plots_new  <- paste0(subfolder_new, "plots_all/")

if (!dir.exists(subfolder_new)) {
       dir.create(subfolder_new)
}
if (!dir.exists(subfolder_data_new)) {
       dir.create(subfolder_data_new)
}
if (!dir.exists(subfolder_plots_new)) {
       dir.create(subfolder_plots_new)
}

##################################################################
##################################################################
###################### Processing Data
##################################################################
##################################################################
method_names <- c(
  "var_standard",
  "mvar_standard", "mvar_adaptive",
  "modvar_cv_bsubj", "modvar_ada.cv")
method_names_clean <- c(
  "VAR",
  "M-VAR: Standard", "M-VAR: Adaptive",
  "MOD-VAR: CV", "Adap. MOD-VAR: CV")
method_lookup <- setNames(method_names_clean, method_names)

highlight_methods <- c("MOD-VAR: CV", "Adap. MOD-VAR: CV")
method_colors <- c(
  "VAR" = cbPalette[1],
  "M-VAR: Standard" = cbPalette[2],
  "M-VAR: Adaptive" = cbPalette[3],
  "MOD-VAR: CV" = "#000000",
  "Adap. MOD-VAR: CV" = "#D55E00"
)

for(sigma2_val in c(0.1)) {
  for (d_val in c(10, 20)) {
  
    ##############################
    ##############################
    ## LOADING ALL DATA WITH T0 = P.
    sim_ind_load    <- which(d == d_val, sigma2 == sigma2_val) #p == p_val)
    type            <- "all"
    results_dir     <- paste0(subfolder_new, "plots_", type, "/")

    for (sim_ind in sim_ind_load) {
      load(paste0(
        subfolder_new, "data_all/",
        "output", sim_ind, ".RData"))
    }

    ##############################
    ##############################
    ## TRUE POSITIVE RATE OF JIC-HD METHODS:

    ## Merge dataset of JIC-HD-derived data.
    output_merged  <- distinct(bind_rows(mget(ls(pattern = '^output\\d+')))) %>%
      filter(method %in% method_names) %>%
      mutate(method = recode(method, !!!setNames(method_names_clean, method_names))) %>%
      filter(sigma2 == sigma2_val, signed == TRUE) %>%
      select(-sigma2) %>%
      group_by(d, p, T, n, prob_c, method) %>%

    summarise(
      meanTPR = mean(TPR),
      sdTPR = sd(TPR),
      meanFPR = mean(FPR),
      sdFPR = sd(FPR),
      meanTime = mean(time),
      sdTime = sd(time),
      meanL1 = mean(l1),
      sdL1 = sd(l1),
      meanFr = mean(Fr),
      sdFr = sd(Fr),
      meanfc1 = mean(fcast1),
      sdfc1 = sd(fcast1),
      meanfc2 = mean(fcast2),
      sdfc2 = sd(fcast2),
      meanfc3 = mean(fcast3),
      sdfc3 = sd(fcast3),
      meanSens = mean(sens),
      sdSens = sd(sens),
      meanSpec = mean(spec),
      sdSpec = sd(spec)
      )

    print(table(output_merged$method))

    highlight_data <- output_merged %>% filter(method %in% highlight_methods)
    other_data <- output_merged %>% filter(!(method %in% highlight_methods))

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_TPR_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p1 <- output_merged %>% ggplot(aes(x = T, y = meanTPR)) +
      geom_ribbon(aes(ymin = meanTPR - sdTPR, ymax = meanTPR + sdTPR, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      geom_hline(yintercept = c(0, 1), linetype = 2) +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p1)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_FPR_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p2 <- output_merged %>% ggplot(aes(x = T, y = meanFPR)) +
      geom_ribbon(aes(ymin = meanFPR - sdFPR, ymax = meanFPR + sdFPR, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      geom_hline(yintercept = c(0), linetype = 2) +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p2)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_L1_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p3 <- output_merged %>% ggplot(aes(x = T, y = meanL1)) +
      geom_ribbon(aes(ymin = meanL1 - sdL1, ymax = meanL1 + sdL1, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      geom_hline(yintercept = c(0), linetype = 2) +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p3)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_Fr_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p4 <- output_merged %>% ggplot(aes(x = T, y = meanFr)) +
      geom_ribbon(aes(ymin = meanFr - sdFr, ymax = meanFr + sdFr, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      geom_hline(yintercept = c(0), linetype = 2) +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p4)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_Time_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p5 <- output_merged %>%
      mutate(logmt = log(meanTime, base = 60)) %>%
      ggplot(aes(x = T, y = logmt)) +
      geom_line(data = highlight_data %>% mutate(logmt = log(meanTime, base = 60)), aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data %>% mutate(logmt = log(meanTime, base = 60)), aes(color = method), linewidth = 0.9, linetype = "dashed") +
      ylab(expression(log[60]("seconds"))) +
      geom_hline(yintercept = c(0, 1, 1.56, 2), linetype = 2, linewidth = 0.25) +
      annotate("text", x = 55, y = 0.15, label = "1 s", size = 2.5) +
      annotate("text", x = 55 * 0.93, y = 1.15, label = "1 m", size = 2.5) +
      annotate("text", x = 55 * 0.93, y = 1.71, label = "10 m", size = 2.5) +
      annotate("text", x = 55 * 0.93, y = 2.15, label = "1 h", size = 2.5) +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p5)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_fc1_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p6 <- output_merged %>% ggplot(aes(x = T, y = meanfc1)) +
      geom_ribbon(aes(ymin = meanfc1 - sdfc1, ymax = meanfc1 + sdfc1, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      geom_hline(yintercept = c(sqrt(sigma2_val)), linetype = 2, col = "red") +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p6)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_fc2_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p7 <- output_merged %>% ggplot(aes(x = T, y = meanfc2)) +
      geom_ribbon(aes(ymin = meanfc2 - sdfc2, ymax = meanfc2 + sdfc2, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      geom_hline(yintercept = c(sqrt(sigma2_val)), linetype = 2, col = "red") +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p7)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_fc3_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p8 <- output_merged %>% ggplot(aes(x = T, y = meanfc3)) +
      geom_ribbon(aes(ymin = meanfc3 - sdfc3, ymax = meanfc3 + sdfc3, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      geom_hline(yintercept = c(sqrt(sigma2_val)), linetype = 2, col = "red") +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p8)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_Sens_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p9 <- output_merged %>% ggplot(aes(x = T, y = meanSens)) +
      geom_ribbon(aes(ymin = meanSens - sdSens, ymax = meanSens + sdSens, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p9)
    dev.off()

    file_name <- paste0(
      subfolder_plots_new,
      "v", version,
      "_542_Spec_ByPnN",
      "_d", d_val,
      "_sigma2_", sigma2_val,
      ".pdf")
    pdf(file_name, width = 8, height = 5)
    p10 <- output_merged %>% ggplot(aes(x = T, y = meanSpec)) +
      geom_ribbon(aes(ymin = meanSpec - sdSpec, ymax = meanSpec + sdSpec, fill = method), alpha = 0.1) +
      geom_line(data = highlight_data, aes(color = method), linewidth = 1.4, linetype = "solid") +
      geom_line(data = other_data, aes(color = method), linewidth = 0.9, linetype = "dashed") +
      geom_hline(yintercept = c(0), linetype = 2) +
      facet_grid(p ~ prob_c + n, scales = "free_x", labeller = label_parsed) +
      scale_colour_manual(values = method_colors) +
      scale_fill_manual(values = method_colors) +
      theme(legend.position = "bottom")
    print(p10)
    dev.off()

    rm(list = grep('^output\\d+', ls(), value = TRUE))
    rm(list = grep("args", ls(), value = TRUE))
    rm(output_merged, file_name, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)

  }
}
