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
library(gridExtra)


###################### Creating folders:
subfolder_new        <- paste0("600_AggregatedDataReal_v3/")
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
  "mvar_adaptive",    "mvar_standard",    
  "modvar_MO",
  "ada.modvar_MO")
method_names_clean <- c(
  "VAR",          
  "MultiVAR: Standard", "MultiVAR: Adaptive",
  "MOD-VAR: Standard",
  "MOD-VAR: Adaptive") 

results_dir     <- paste0(subfolder_new, "plots_all/")

##############################
## LOADING ALL DATA
for(ind in 1:9) {
  load(paste0(
    subfolder_new, "data_all/",
    "output", ind, ".RData"))
}

##############################
##############################

## Create an exact lookup table
method_lookup <- setNames(method_names_clean, method_names)

## Merge dataset of JIC-HD-derived data.
output_merged  <- distinct(bind_rows(mget(ls(pattern = '^output\\d+')))) %>%
  filter(method %in% method_names) %>%
  mutate(method = recode(method, !!!method_lookup)) %>%
  group_by(T, method) %>%

  summarise(
    meanfc1 = mean(rel_msfe1),
    sdfc1 = sd(rel_msfe1),
    meanfc2 = mean(rel_msfe2),
    sdfc2 = sd(rel_msfe2),
    meanfc3 = mean(rel_msfe3),
    sdfc3 = sd(rel_msfe3),
    meanfc4 = mean(rel_msfe4),
    sdfc4 = sd(rel_msfe4),
    meanfc5 = mean(rel_msfe5),
    sdfc5 = sd(rel_msfe5),
    meanfc6 = mean(rel_msfe6),
    sdfc6 = sd(rel_msfe6),
    meanfc7 = mean(rel_msfe7),
    sdfc7 = sd(rel_msfe7),
    meanfc8 = mean(rel_msfe8),
    sdfc8 = sd(rel_msfe8),
    meanfc9 = mean(rel_msfe9),
    sdfc9 = sd(rel_msfe9),
    meanfc10 = mean(rel_msfe10),
    sdfc10 = sd(rel_msfe10)
  )


file_name <- paste0(
  subfolder_plots_new, 
  "622_forecast.pdf")
pdf(file_name, width = 8, height = 5)
p1 <-  output_merged %>% ggplot(aes(x = T, y = meanfc1)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc1 - sdfc1, ymax = meanfc1 + sdfc1, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p1)
    
p2 <-  output_merged %>% ggplot(aes(x = T, y = meanfc2)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc2 - sdfc2, ymax = meanfc2 + sdfc2, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p2)
    
p3 <-  output_merged %>% ggplot(aes(x = T, y = meanfc3)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc3 - sdfc3, ymax = meanfc3 + sdfc3, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p3)
        
p4 <-  output_merged %>% ggplot(aes(x = T, y = meanfc4)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc4 - sdfc4, ymax = meanfc4 + sdfc4, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p4)
    
p5 <-  output_merged %>% ggplot(aes(x = T, y = meanfc5)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc5 - sdfc5, ymax = meanfc5 + sdfc5, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p5)
    
p6 <-  output_merged %>% ggplot(aes(x = T, y = meanfc6)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc6 - sdfc6, ymax = meanfc6 + sdfc6, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p6)
    
p7 <-  output_merged %>% ggplot(aes(x = T, y = meanfc7)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc7 - sdfc7, ymax = meanfc7 + sdfc7, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p7)
        
p8 <-  output_merged %>% ggplot(aes(x = T, y = meanfc8)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc8 - sdfc8, ymax = meanfc8 + sdfc8, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p8)
    
p9 <-  output_merged %>% ggplot(aes(x = T, y = meanfc9)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc9 - sdfc9, ymax = meanfc9 + sdfc9, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p9)
    
p10 <-  output_merged %>% ggplot(aes(x = T, y = meanfc10)) +
  geom_line(aes(col = method, linetype = method), linewidth = 1) +
  geom_ribbon(aes(ymin = meanfc10 - sdfc10, ymax = meanfc10 + sdfc10, fill = method), alpha = 0.1) +
  theme(legend.position="bottom")
print(p10)
        

dev.off()


file_name <- paste0(
  subfolder_plots_new, 
  "623_forecast.pdf")
pdf(file_name, width = 6.5, height = 4)

library(tidyverse)

plot_data <- output_merged %>%
  filter(T <= 100) %>%
  select(T, method, 
         meanfc1, sdfc1,
         meanfc2, sdfc2,
         meanfc3, sdfc3) %>%
  pivot_longer(
    cols = -c(T, method),
    names_to = c(".value", "horizon"),
    names_pattern = "(meanfc|sdfc)([1-3])"
  ) %>%
  mutate(
    horizon = factor(
      horizon,
      levels = c("1", "2", "3"),
      labels = c("1-step", "2-step", "3-step")
    )
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "VAR",
        "MultiVAR: Standard",
        "MultiVAR: Adaptive",
        "MOD-VAR: Standard",
        "MOD-VAR: Adaptive"
      )
    )
  )
  
cbPalette <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)
lw1 <- 1.1
lw2 <- 1.5
ggplot(
  plot_data,
  aes(
    x = T,
    y = meanfc,
    color = method,
    fill = method,
    linetype = method
  )
) +
  geom_ribbon(
    aes(
      ymin = meanfc - sdfc,
      ymax = meanfc + sdfc
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line(aes(linewidth = method)) +
  facet_wrap(~ horizon, nrow = 1) +
  scale_color_manual(
    values = c(
      "VAR" = cbPalette[1],
      "MultiVAR: Standard" = cbPalette[2],
      "MultiVAR: Adaptive" = cbPalette[3],
      "MOD-VAR: Standard" = "red",
      "MOD-VAR: Adaptive" = "black"
    )
  ) +
  scale_fill_manual(
    values = c(
      "VAR" = cbPalette[1],
      "MultiVAR: Standard" = cbPalette[2],
      "M-VAR: Adaptive" = cbPalette[3],
      "MOD-VAR: Standard" = "red",
      "MOD-VAR: Adaptive" = "black"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "VAR" = "dotted",
      "MultiVAR: Standard" = "dashed",
      "MultiVAR: Adaptive" = "dotdash",
      "MOD-VAR: Standard" = "solid",
      "MOD-VAR: Adaptive" = "solid"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      "VAR" = lw1,
      "MultiVAR: Standard" = lw1,
      "MultiVAR: Adaptive" = lw1,
      "MOD-VAR: Standard" = lw2,
      "MOD-VAR: Adaptive" = lw2
    )
  ) +
  labs(
    x = "Time",
    y = "Forecasting Error",
    color = "Method",
    linetype = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  guides(
  linewidth = "none",
  fill = "none",
  color = guide_legend(
    nrow = 2,
    byrow = TRUE,
    override.aes = list(
      linewidth = c(lw1, lw1, lw1, lw2, lw2)
    )
  ),
  linetype = guide_legend(
    nrow = 2,
    byrow = TRUE,
    override.aes = list(
      linewidth = c(lw1, lw1, lw1, lw2, lw2)
    )
  )
)

dev.off()

rm(list = grep('^output\\d+', ls(), value = TRUE))
rm(
  output_merged, file_name, 
  p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)

