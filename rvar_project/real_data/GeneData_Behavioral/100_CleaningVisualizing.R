
library(magrittr)
library(readxl)
library(tidyverse)

df <- read_xlsx("GeneData_Behavioral/bigger_daily_ass.xlsx")

dim(df)
colnames(df)

subj_ids <- unique(df$id)
n <- length(unique(df$id))

sum(is.na(df))
sum(is.na(df)) / (nrow(df) * ncol(df))

df_list <- subj_ids %>%
  lapply(function(x) { df[df$id == x, ] })

sapply(df_list, dim)[1, ] %>% plot(ylim = c(0,150))  


sapply(df_list, function(x) {sum(is.na(x))})


colnames(df)


dates_list <- df_list %>% 
  lapply(function(x) {as.Date(x$true_date)})
str(dates_list)

no_rep_dates <- lapply(
  dates_list, 
  function(x) {length(x) - length(unique(x))}) %>%
  unlist()
no_rep_dates

par(mfrow = c(5,4))
df_list %>% 
  lapply(function(x) {as.Date(x$true_date)}) %>%
  #{.[no_rep_dates > 0]} %>%
  lapply(plot)

df_list %>%
  lapply(
    function(x) {
      t <- length(x)
      
    }
  ) %>% unlist()



## Removing all variables except psychological + date.
df1_sel <- df_list[[1]] %>% 
  select(
    -gender, -gpa, -employed,
    -id, -time, true_date, -DateTime,
    -m.rt, sd.rt,              
    -falsestarts, -cv.rt, -sd.rt, -va_acc,
    -sum_flanker, -total_flanker, -rt_mean_flnkr,
    -date_time, -feedback, -date, -data_set, 
    -(sleep_time:o2_max), -attention_check)

## Imputing missing values
apply(df1_sel, 2, function(x) {sum(is.na(x))})
median_vals <- apply(df1_sel[,-1], 2, function(x) {median(x[!is.na(x)])})
missing_entries <- which(is.na(df1_sel), arr.ind = TRUE)

df1_sel[is.na(df1_sel[,2]),2] <- median_vals[1]
df1_sel[is.na(df1_sel[,8]),8] <- median_vals[7]

sum(is.na(df1_sel))

## Visualizing data:
n <- nrow(df1_sel)
#df1_norm <- scale(df1_sel[, -1])
#df1_norm$dates <- as.Date(df1_sel$true_date)
#df1_pca <- prcomp(df1_norm)


library(ggplot2)

df_long <- df1_sel %>%
  pivot_longer(
    cols = -true_date,
    names_to = "variable",
    values_to = "value")

ggplot(df_long, aes(x = true_date, y = value, group = variable)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw()

##########################################
##########################################
##########################################

for (ind in 1:3) {
  
  ind = 3
  df_ind <- df_list[[ind]]
  df_sel <- df_ind %>% 
    select(
      -gender, -gpa, -employed,
      -id, -time, true_date, -DateTime,
      -m.rt, sd.rt,              
      -falsestarts, -cv.rt, -sd.rt, -va_acc,
      -sum_flanker, -total_flanker, -rt_mean_flnkr,
      -date_time, -feedback, -date, -data_set, 
      -(sleep_time:o2_max), -attention_check)
  
  median_vals <- apply(df_sel[,-1], 2, function(x) {median(x[!is.na(x)])})
  
  for(col_ind in 2:ncol(df_sel)) {
    missing <- is.na(df_sel[ , col_ind]) 
    df_sel[missing, col_ind] <- median_vals[col_ind - 1]
  }
  
  df_sel[, -1] <- scale(df_sel[, -1])
  
  df_long <- df_sel %>%
    pivot_longer(
      cols = -true_date,
      names_to = "variable",
      values_to = "value")
  
  ggplot(df_long, aes(x = true_date, y = value, group = variable)) +
    geom_line() +
    facet_wrap(~ variable, scales = "free_y") +
    theme_bw()

}


df_sel_list <- lapply(
  df_list,
  function(df_ind) {
    df_sel <- df_ind %>% 
      select(
        -gender, -gpa, -employed,
        -id, -time, true_date, -DateTime,
        -m.rt, sd.rt,              
        -falsestarts, -cv.rt, -sd.rt, -va_acc,
        -sum_flanker, -total_flanker, -rt_mean_flnkr,
        -date_time, -feedback, -date, -data_set, 
        -(sleep_time:o2_max), -attention_check)
    
    median_vals <- apply(df_sel[,-1], 2, function(x) {median(x[!is.na(x)])})
    
    for(col_ind in 2:ncol(df_sel)) {
      missing <- is.na(df_sel[ , col_ind]) 
      df_sel[missing, col_ind] <- median_vals[col_ind - 1]
    }
    
    df_sel[, -1] <- scale(df_sel[, -1])
    
    return(df_sel)
  }
)

library(BigVAR)
library(plot.matrix)


par(mfrow = c(3,3),
    mar = c(5.1, 4.1, 4.1, 4.1))
for (ind in 1:18) {
  
  model <- BigVAR::constructModel(
    Y = as.matrix(df_sel_list[[ind]][, -1]), 
    p = 1,
    gran = c(50,10),
    struct = "Basic",
    cv = "Rolling",
    verbose = TRUE,
    ownlambdas = FALSE,
    model.controls=list(intercept=FALSE),
    linear = FALSE)
  
  
  results <- BigVAR::cv.BigVAR(model)
  xlim = max(abs(as.matrix(coef(results))))
  if (xlim == 0) xlim = 0.1
  plot(as.matrix(coef(results)), 
       breaks = seq(-1.1 * xlim, 1.1 * xlim, length.out = 20),
       main = paste0("Subject:", ind))

}

table(df$gender)
summary(df$gpa)
summary(df$employed)

  
    