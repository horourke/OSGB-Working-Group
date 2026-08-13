library(magrittr)
library(readxl)
library(tidyverse)
library(plot.matrix)

df <- read_xlsx("GeneData_Behavioral/bigger_daily_ass.xlsx")

##############################################
##############################################
## We will analyze the relationship between:
##  Moderators:
##    GPA, Employment, Gender.
##  Time series data:
##    agreeableness, conscientiousness, 
##    emotional_stability, extraversion,
##    openness, depression.

dim(df)
colnames(df)

subj_ids <- unique(df$id)
n <- length(unique(df$id))

sum(is.na(df))
sum(is.na(df)) / (nrow(df) * ncol(df))

df_list <- subj_ids %>%
  lapply(function(x) { df[df$id == x, ] })


## Step 1: Lets calculate the average score for each subject:

X <- df %>%
  select(id, gender, gpa, employed) %>%
  group_by(id) %>%
  summarize_all(mean)

mean_no_missing <- function(x) {
  mean(x, na.rm = TRUE)
}

Y <- df %>%
  select(id, agreeableness, conscientiousness, 
         emotional_stability, extraversion,
         openness, depression)  %>%
  group_by(id) %>% 
  summarize_all(mean_no_missing)

##############################################
##############################################
##############################################
##############################################
## Cognitive:
colnames(df_list[[1]])

par(mfrow = c(4,4))
df_list[1:16] %>%
  lapply(
    function(x) {
      x %>%
        select(cv.rt, m.rt, va_acc, sleep_time, anxiety, depression) %>%
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
             dp = depression,)
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
        plot(breaks = 20)
    }
  )




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
