

install.packages("tseries")
install.packages("MTS")
install.packages("ppcor")
install.packages("TSA")
install.packages("forecast")

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

##############################################
##############################################
## We will analyze the relationship between:
##  Moderators:
##    GPA, Employment, Gender.
##  Time series data:
##    agreeableness, conscientiousness, 
##    emotional_stability, extraversion,
##    openness, depression.

df <- read_xlsx("GeneData_Behavioral/bigger_daily_ass.xlsx")
dim(df)
colnames(df)

subj_ids <- unique(df$id)
n <- length(unique(df$id))

sum(is.na(df))
sum(is.na(df)) / (nrow(df) * ncol(df))

df_list <- subj_ids %>%
  lapply(function(x) { 
    df[df$id == x, ] %>% 
      mutate(date = as.Date(true_date)) %>%
      return()
    })

## Step 1: Lets calculate the average score for each subject:

X <- df %>%
  dplyr::select(id, gender, gpa, employed) %>%
  group_by(id) %>%
  summarize_all(mean)

mean_no_missing <- function(x) {
  mean(x, na.rm = TRUE)
}

Y <- df %>%
  dplyr::select(id, cv.rt, m.rt, va_acc, sleep_time, anxiety, depression)  %>%
  group_by(id) %>% 
  summarize_all(mean_no_missing)




##############################################
##############################################
##############################################
##############################################
## Cognitive:
colnames(df_list[[1]])

missing <- df_list %>% lapply(
  function(x) {
    x %>%
      dplyr::select(date, cv.rt, m.rt, va_acc, sleep_time, anxiety, depression) %>%
      return()
  }
) %>%
  lapply(
    function(x) {
      rename(x,
             cv = cv.rt,
             mrt =  m.rt,
             ac = va_acc,
             sl = sleep_time,
             an = anxiety,
             dp = depression)
    }
  ) %>%
  lapply(
    function(data) {
      date <- as.Date(data$date)
      ts <- data %>% dplyr::select(-date)
    
      apply(ts, 2, function(x){var(x[!is.na(x)]) < 1e-6}) %>%
        {sum(.) > 0}
    }) %>% unlist() %>% which()

s <- 1
p_vals <- df_list[-missing] %>% #[1:120] %>%
  lapply(
    function(x) {
      x %>%
        dplyr::select(date, cv.rt, m.rt, va_acc, sleep_time, anxiety, depression) %>%
        return()
    }
  ) %>%
  lapply(
    function(x) {
      rename(x,
             cv = cv.rt,
             mrt =  m.rt,
             ac = va_acc,
             sl = sleep_time,
             an = anxiety,
             dp = depression)
    }
  ) %>%
  lapply(
    function(data) {
      date <- as.Date(data$date)
      print(class(date))
      ts <- data %>% dplyr::select(-date)
      
      print(apply(ts, 2, function(x){var(x[!is.na(x)])}))
      
      for (i in 1:ncol(ts)) {
        ts.isna <- is.na(ts[ , i]) %>%
          as.vector()
        med_val <- ts[, i] %>%
          pull(1) %>%
          median(na.rm = TRUE)
        ts[ts.isna ,i] <- med_val
      }
      
      #return(lapply(data, function(x) {adf.test(x, k = 2, alternative = "stationary")}))
      print(i)
      return(mq_na(ts, date, lag = s, adj = 0, plot = FALSE))
      #return(ccm(data, lag = s, output = FALSE))
      #output <- mq_na(x, date, lag = 24, adj = 0, plot = TRUE)
      
      #n <- nrow(data)
      #data1 <- data[1:(n - s), ]
      #data2 <- data[(s+1):n, ]
      
      #cor(data1, data2) %>%
      #  abs() %>%
      #  plot(breaks = seq(0,1, length.out = 20))
    }
  )
str(p_vals[[1]][,4])
p_vals[[1]]$ccm

pv <- c()
for(i in 1:length(p_vals)) {
  pv <- c(pv, p_vals[[i]][,4])  
}

par(mfrow = c(1,1))
plot(pv)
abline(h = 0.05)
abline(v = seq(0.5, length(df_list) * s + 0.5, length.out = length(df_list) + 1))

## Conclusion: it seems like approximately 10 patients out of 54 have a p-value below 0.05.
## Conclusion: it seems like approximately 12 patients out of 54 have a p-value below 0.1.

plot(p.adjust(pv, method = "BH"))
abline(h = 0.05)

## After multiple testing adjustment, we have only 3 patients well below 0.05.
## So, not great really.

## Lets see what happens when we look at the neuro-data.

##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


par(mfrow = c(4,4))
df_list[1:16] %>%
  lapply(
    function(x) {
      x %>%
        select(agreeableness, conscientiousness, 
               emotional_stability, extraversion,
               openness) %>%
        return()
    }
  ) %>%
  lapply(
    function(x) {
      rename(x,
        agr = agreeableness, 
        con = conscientiousness,
        est = emotional_stability,
        ext = extraversion,
        opn = openness)
    }
  ) %>%
  lapply(
    function(data) {
      
      for (i in 1:ncol(data)) {
        x.isna <- is.na(data[ , i]) %>%
          as.vector()
        med_val <- data[, i] %>%
          pull(1) %>%
          median(na.rm = TRUE) %>%
          print()
        data[x.isna ,i] <- med_val
      }
      
      n <- nrow(data)
      data1 <- data[1:(n - 1), ]
      data2 <- data[2:n, ]
      
      cor(data1, data2) %>%
        abs() %>%
        plot(breaks = seq(0,1, length.out = 20))
    }
  )


##############################################
##############################################

par(mfrow = c(4,4))
df_list[1:16] %>%
  lapply(
    function(x) {
      x %>%
        select(agreeableness, conscientiousness, 
               emotional_stability, extraversion,
               openness) %>%
        return()
    }
  ) %>%
  lapply(
    function(x) {
      rename(x,
             agr = agreeableness, 
             con = conscientiousness,
             est = emotional_stability,
             ext = extraversion,
             opn = openness)
    }
  ) %>%
  lapply(
    function(data) {
      for (i in 1:ncol(data)) {
        x.isna <- is.na(data[ , i]) %>%
          as.vector()
        med_val <- data[, i] %>%
          pull(1) %>%
          median(na.rm = TRUE) %>%
          print()
        data[x.isna ,i] <- med_val
      }
      
      pca <- prcomp(data) 
      return(pca$x[, 1:2] %*% t(pca$rotation[, 1:2]))
    }
  ) %>%
  lapply(
    function(data) {

      n <- nrow(data)
      data1 <- data[1:(n - 1), ]
      data2 <- data[2:n, ]
      
      cor(data1, data2) %>%
        abs() %>%
        plot()
    }
  )
