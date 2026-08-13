
library(magrittr)
library(readxl)
library(tidyverse)

df <- read.csv("BrainNetworksGeneData_Systematic/reduced_network_signals_data.csv")

dim(df)
colnames(df)
table(df$sub)

head(df)


subj_ids <- unique(df$sub)
n <- length(unique(df$sub))

sum(is.na(df))
sum(is.na(df)) / (nrow(df) * ncol(df))

df_list <- subj_ids %>%
  lapply(function(x) { df[df$sub == x, ] }) %>%
  lapply(function(x) x %>% mutate(t = 1:nrow(x)) %>% return()) %T>%
  str()
  
df <- do.call("rbind", df_list)




##########################################
##########################################
##########################################

library(dplyr)
library(tidyr)
library(ggplot2)

## Convert to long format
dat_long <- df %>%
  pivot_longer(
    cols = -c(sub, t),
    names_to = "variable",
    values_to = "value"
  )

## Keep only the first 5 subjects
dat_long <- dat_long %>%
  filter(sub %in% sort(unique(sub))[1:5])

## Variable groups (replace these names with yours if needed)
vars1 <- unique(dat_long$variable)[1:4]
vars2 <- unique(dat_long$variable)[5:8]
vars3 <- unique(dat_long$variable)[9:12]

## Plot 1
p1 <- dat_long %>%
  filter(variable %in% vars1) %>%
  ggplot(aes(t, value, group = 1)) +
  geom_line() +
  geom_point(size = 0.8) +
  facet_grid(sub ~ variable, scales = "free_y") +
  theme_bw() +
  labs(x = "Time", y = NULL)

## Plot 2
p2 <- dat_long %>%
  filter(variable %in% vars2) %>%
  ggplot(aes(t, value, group = 1)) +
  geom_line() +
  geom_point(size = 0.8) +
  facet_grid(sub ~ variable, scales = "free_y") +
  theme_bw() +
  labs(x = "Time", y = NULL)

## Plot 3
p3 <- dat_long %>%
  filter(variable %in% vars3) %>%
  ggplot(aes(t, value, group = 1)) +
  geom_line() +
  geom_point(size = 0.8) +
  facet_grid(sub ~ variable, scales = "free_y") +
  theme_bw() +
  labs(x = "Time", y = NULL)

## Display
p1
p2
p3
rm(p1, p2, p3)

##########################################
##########################################
##########################################

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
##############################################
##############################################
## Cognitive:
par(mfrow = c(3,3))
df_list[1:9] %>%
  lapply(function(x) {
    n <- nrow(x) -1
    X1 <- x[1:n, ] %>% {.[, -c(1, 14)]} %T>% {print(head(.))}
    X2 <- x[2:(n+1), ] %>% {.[, -c(1, 14)]} %T>% {print(head(.))}
    
    plot(cor(X1, X2))
    
  })
## Autocorrelations: non-negligible and heterogeneous. 

par(mfrow = c(3,3))
df_list[1:9] %>%
  lapply(function(x) {
    n <- nrow(x) - 1
    
    x <- x[, -c(1,14)]
    plot(cor(x))
  })
## Covariance matrix: 

par(mfrow = c(5,5))
df_list[1:25] %>%
  lapply(function(x) {
    n <- nrow(x) -1
    x <- x[, -c(1,14)]
    cmat <- cor(x)
    
    
    plot(log(eigen(cmat)$values, base = 10))
  })


par(mfrow = c(1,1))
df_list %>%  
  lapply(function(x) { x[,-c(1,14)] %>% cor() %>% eigen() }) %>%
  lapply(function(x) return(min(x$values))) %>%
  unlist() 
## It seems like most covariances are rank deficient. I don't know why...
par(mfrow = c(1,1))
vectors <- df_list %>%  
  lapply(function(x) { x[,-c(1,14)] %>% cov() %>% eigen() }) %>%
  lapply(function(x) return((x$vectors[, length(x$values)]))) %>%
  {do.call("rbind", .)} 
plot(abs(vectors))
abs(vectors[1,])


colnames(df)
vectors[1, c(2,7,8,9,10)]
(vectors[1,2] * df$DMN_Global + vectors[1,7] * df$DMN_mPFC + vectors[1,8] * df$DMN_PCC_Precuneus + 
  vectors[1,9] * df$DMN_IPL + vectors[1,10] * df$DMN_MTL) %>% abs() %>% log(base = 10) %>% plot()
## Yes: there is collinearity because of the variable global.


## We see that there are 
## A) High degree of autocorrelation
## B) High degree of heterogeneity
## C) It can be removed by eliminating the variable DMN_Global

df <- df %>% tibble() %>% dplyr::select(-DMN_Global)
df_list <- df_list %>% lapply(function(x) return(x%>% dplyr::select(-DMN_Global)))

s <- 20
p_vals <- df_list %>%
  lapply(
    function(data) {
      
      pvals <- data %>% {.[, -c(1,14)]}  %>%
        ccm(lag = s, output = FALSE) %T>%
        ##adf.test(k = 2, alternative = "stationary")
        ##mq(lag = 3, adj = 0) %T>%
        print()
      return(pvals)
      
      #return(lapply(function(x) {adf.test(x, k = 2, alternative = "stationary")}))
      
    }
  )

str(p_vals)

pv <- c()
for(i in 1:length(p_vals)) {
  pv <- c(pv, p_vals[[i]]$pvalue)  
}

par(mfrow = c(1,1))
plot(pv)
abline(h = 0.05)
abline(v = seq(0.5, length(df_list) * s + 0.5, length.out = length(df_list) + 1))

## Conclusion: even at the 10-time step, there is autocorrelation among the variables. 


ls()
rm(dat_long, i,n, p_vals, pv, s, vars1, vars2, vars3, vectors)
save.image("BrainNetworksGeneData_Systematic/100_Data.RData")


rm(list = ls())



