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
## Visualize the correlation between X and Y:

par(mar= c(5.1, 4.1, 4.1, 5.1))

Z <- X %>%
  inner_join(Y, by = "id") %>%
  filter(gender %in% c(1,2)) %T>%
  {print(colnames(.))}
mx <- max(abs(cor(X[X$gender<3, -1], Y[X$gender<3, -1])), na.rm = TRUE)

plot(
  cor(X[X$gender<3, -1], Y[X$gender<3, -1]),
  col = colorRampPalette(c("red", "white", "green"))(11),
  breaks = seq(-mx, mx, length.out = 12)
)

lm(agreeableness ~ gender + gpa + employed, data = Z) %>%
  summary()
lm(agreeableness ~ gender + gpa + employed, data = Z) %>%
  plot()

lm(conscientiousness ~ gender + gpa + employed, data = Z) %>%
  summary()
lm(conscientiousness ~ gender + gpa + employed, data = Z) %>%
  plot()

lm(emotional_stability ~ gender + gpa + employed, data = Z) %>%
  summary()
lm(emotional_stability ~ gender + gpa + employed, data = Z) %>%
  plot()

lm(extraversion ~ gender + gpa + employed, data = Z) %>%
  summary()
lm(extraversion ~ gender + gpa + employed, data = Z) %>%
  plot()

lm(openness ~ gender + gpa + employed, data = Z) %>%
  summary()
lm(openness ~ gender + gpa + employed, data = Z) %>%
  plot()

lm(depression ~ gender + gpa + employed, data = Z) %>%
  summary()
lm(depression ~ gender + gpa + employed, data = Z) %>%
  plot()

## Conclusions:  
## -> It seems that employment is the moderator
##      with the highest correlation to emotional
##      variables. Gender/GPA are not strong.
## -> Employment negatively associated to:
##      aggreableness, concientiousness, extraversion,
##      openness. 
## -> Positively associated with depression.
## 
##############################################
##############################################

Z %>% select(-id, -depression, -gender, -gpa, -employed) %>%
  pairs()
pca <- Z %>% 
  select(-id, -depression, -gender, -gpa, -employed) %>%
  prcomp()


pairs(pca$x[,1:3])

Z %>% 
  select(-id, -depression, -gender, -gpa, -employed) %>%
  {sum(diag(cov(.)))}
plot(pca)


mx <- max(abs(pca$rotation), na.rm = TRUE)
plot(
  pca$rotation,
  col = colorRampPalette(c("red", "white", "green"))(100),
  breaks = seq(-mx, mx, length.out = 101)
)



##############################################
##############################################
## Visualize the correlation between X and Y:

Xf <- df %>%
  filter(gender < 3) %>%
  select(gender, gpa, employed) 

Yf <- df %>%
  filter(gender < 3) %>%
  select(agreeableness, conscientiousness, 
         emotional_stability, extraversion,
         openness, depression) %>%
  mutate(
    across(
      where(is.numeric),
      ~ replace(.x, is.na(.x), median(.x, na.rm = TRUE))
    )
  )
Zf <- cbind(Xf, Yf)

mx <- max(abs(cor(Xf, Yf)), na.rm = TRUE)

plot(
  cor(Xf, Yf),
  col = colorRampPalette(c("red", "white", "green"))(100),
  breaks = seq(-mx, mx, length.out = 101)
)


lm(agreeableness ~ gender + gpa + employed, data = Zf) %>%
  summary()

lm(conscientiousness ~ gender + gpa + employed, data = Zf) %>%
  summary()

lm(emotional_stability ~ gender + gpa + employed, data = Zf) %>%
  summary()

lm(extraversion ~ gender + gpa + employed, data = Zf) %>%
  summary()

lm(openness ~ gender + gpa + employed, data = Zf) %>%
  summary()

lm(depression ~ gender + gpa + employed, data = Zf) %>%
  summary()

dim(Zf)

Zf %>% select(-depression, -gender, -gpa, -employed) %>%
  pairs()
pcaf <- Zf %>% 
  select(-depression, -gender, -gpa, -employed) %>%
  prcomp()

Zf %>% 
  select(-depression, -gender, -gpa, -employed) %>%
  {sum(diag(cov(.)))}

plot(pcaf)

pairs(pcaf$x[,1:3])



##############################################
##############################################
## Compare PCA directions for Individual + Joint
## PCA.

par(mfrow = c(2,1))
mx <- max(abs(pca$rotation), na.rm = TRUE)
PCRM <- -pca$rotation
PCRM[,4] <- -PCRM[,4]
plot(
  PCRM,
  col = colorRampPalette(c("red", "white", "green"))(100),
  breaks = seq(-mx, mx, length.out = 101))

mx <- max(abs(pcaf$rotation), na.rm = TRUE)
sign <- sign(pcaf$rotation[1,])
      
plot(
  pcaf$rotation,
  col = colorRampPalette(c("red", "white", "green"))(100),
  breaks = seq(-mx, mx, length.out = 101))


mx <- max(abs(pcaf$rotation - PCRM), na.rm = TRUE)
plot(
  abs(pcaf$rotation - PCRM),
  col = colorRampPalette(c("white", "green"))(10),
  breaks = seq(0, mx, length.out = 11))


## Conclusions:  
## -> The leading 2 principal components are shared
## across all subjects. 
## -> The third PC is near-common
## -> The 4th and 5th PC are not common at all.
## 
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
## Compare PCA directions for Individual + Joint
## PCA.


X <- df %>%
  select(id, gender, gpa, employed) %>%
  mutate(id = factor(id)) %>%
  group_by(id) %>%
  summarize_all(mean)

mean_no_missing <- function(x) {
  mean(x, na.rm = TRUE)
}
Y <- df %>%
  select(id, agreeableness, conscientiousness, 
         emotional_stability, extraversion,
         openness, depression, subjective_sleep,   
         anxiety, attention_check,distress,
         somatization)  %>%
  mutate(id = factor(id)) %>%
  rename(
    id = id,
    agr = agreeableness, 
    con = conscientiousness,
    est = emotional_stability,
    ext = extraversion,
    opn = openness,
    dep = depression,
    ssl = subjective_sleep,
    anx = anxiety,
    att = attention_check,
    dis = distress,
    som = somatization) %>%
  group_by(id) %>% 
  summarize_all(mean_no_missing) %>%
  mutate(across(where(is.numeric), ~ as.vector(scale(.))))


plot(
  cor(X[X$gender<3, -1], Y[X$gender<3, -1]),
  col = colorRampPalette(c("red", "white", "green"))(11),
  breaks = seq(-mx, mx, length.out = 12)
)

Z <- X %>%
  inner_join(Y, by = "id") %>%
  filter(gender %in% c(1,2)) %T>%
  {print(colnames(.))}

Z %>% select(-id, -dep, -gender, -gpa, -employed) %>%
  pairs()
pca <- Z %>% 
  select(-id, -depression, -gender, -gpa, -employed) %>%
  prcomp()
plot(pca)

pairs(pca$x[,1:3])


## agr = agreeableness, 
## con = conscientiousness,
## est = emotional_stability,
## ext = extraversion,
## opn = openness,
## dep = depression,
## ssl = subjective_sleep,
## anx = anxiety,
## att = attention_check,
## dis = distress,
## som = somatization

t(t(rownames(pca$rotation)))
mx <- max(abs(pca$rotation), na.rm = TRUE)
plot(
  pca$rotation,
  col = colorRampPalette(c("red", "white", "green"))(100),
  breaks = seq(-mx, mx, length.out = 101)
)



library(tidyverse)

vnames <- rownames(pca$rotation)
pca$rotation %>%
  as.data.frame() %>%
  mutate(variable = rownames(pca$rotation)) %>%
  pivot_longer(-variable,
               names_to = "PC",
               values_to = "loading") %>%
  
  mutate(
    PC = factor(
      PC, levels = paste0("PC", 1:10), ordered = TRUE),
    variable = factor(
      variable, levels = vnames, ordered = TRUE)) %>%
  
  ggplot(aes(x = variable,
             y = loading,
             fill = loading > 0)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "green","FALSE" = "red")) +
  labs(title = "PCA Loadings (Rotation Matrix)",
       x = NULL,
       y = "Loading") +
  theme_minimal()


library(corrplot)
cormat <- X %>%
  inner_join(Y, by = "id") %>%
  select(-id) %>%
  filter(gender %in% c(1,2)) %>% 
  cor()
corrplot(cormat)



##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
##############################################
## 


df_sel <- df %>%
  select(id, agreeableness, conscientiousness, 
         emotional_stability, extraversion,
         openness, depression, subjective_sleep,   
         anxiety, distress,
         somatization) %>%
  rename(
    aggre = agreeableness,
    consc = conscientiousness,
    em_st = emotional_stability,
    extrv = extraversion,
    openn = openness,
    deprs = depression,
    subsl = subjective_sleep,
    anxty = anxiety,
    dstrs = distress,
    somtz = somatization)




df_list <- subj_ids[1:16] %>%
  lapply(function(x) { df_sel[df_sel$id == x, -1] }) %>%
  
  lapply(function(x) {
    x %>%
      mutate(
        across(
          where(is.numeric),
          ~ replace(.x, is.na(.x), median(.x, na.rm = TRUE))
        )
      ) %>%
      return()
  })


df_list %>% 
  lapply(is.na) %>%
  sapply(function(x) {apply(x, 2, sum)}) %>%
  t()
  

df_list %>%
  lapply(is.na) %>%
  lapply(sum) %>%
  unlist()


cor_list <- lapply(df_list, cor)
cov_list <- lapply(df_list, cov)


library(corrplot)
par(mfrow = c(4,4),
    mar = c(3.1, 3.1, 3.1, 5.1))
lapply(cor_list, corrplot)
cov_list %>%
  lapply(
    function(x) {
      mx <- max(abs(x))
      plot(
        x,
        col = colorRampPalette(c("red", "white", "green"))(100),
        breaks = seq(-mx, mx, length.out = 101))
    }
  )


group_data <- df %>%
  select(id, gender, gpa, employed) %>%
  group_by(id) %>%
  summarise_all(
    first)
group_data

par(mfrow = c(1,1))
group_data$gender %>% table()
group_data$employed %>% table()
group_data$gpa %>% hist(breaks = 10)





##############################################
##############################################

df_sel <- df %>%
  select(id, agreeableness, conscientiousness, 
         emotional_stability, extraversion,
         openness, depression, subjective_sleep,   
         anxiety, distress,
         somatization) %>%
  
  rename(
    aggre = agreeableness,
    consc = conscientiousness,
    em_st = emotional_stability,
    extrv = extraversion,
    openn = openness,
    deprs = depression,
    subsl = subjective_sleep,
    anxty = anxiety,
    dstrs = distress,
    somtz = somatization)




df_list <- subj_ids[1:16] %>%
  lapply(function(x) { df_sel[df_sel$id == x, -1] }) %>%
  
  lapply(function(x) {
    x %>%
      mutate(
        across(
          where(is.numeric),
          ~ replace(.x, is.na(.x), median(.x, na.rm = TRUE))
        )
      ) %>%
      return()
  })



par(mfrow = c(2,2))
for(i in 1:4) {
  pc_data <- df_list %>% 
    sapply(
      function(x) {
        pca <- prcomp(x) 
        #return(pca$rotation[, 1] * pca$sdev[1])
        return(pca$rotation[, i])
      }) %>% 
    
    t()  
  
  ipi_mat <- matrix(0, nrow(pc_data), nrow(pc_data))
  for(i in 1:nrow(pc_data)) {
    for(j in 1:nrow(pc_data)) {
      ipi_mat[i,j] <- sum(pc_data[i, ] * pc_data[j, ])    
    }
  }
  plot(abs(ipi_mat), breaks = seq(0,1, length.out = 20))
  
}

par(mfrow = c(2,2))
for(i in 1:4) {
  pc_data <- df_list %>% 
    sapply(
      function(x) {
        pca <- prcomp(x, center = TRUE, scale. = TRUE) 
        return(pca$rotation[, i])
      }) %>% 
    
    t()  
  
  ipi_mat <- matrix(0, nrow(pc_data), nrow(pc_data))
  for(i in 1:nrow(pc_data)) {
    for(j in 1:nrow(pc_data)) {
      ipi_mat[i,j] <- sum(pc_data[i, ] * pc_data[j, ])    
    }
  }
  plot(abs(ipi_mat), 
       breaks = seq(0,1, length.out = 20))
}


dim(pc_data)

## What have you learnt about the data?
##
## 1) The data is extremely heterogeneous across subjects.
## 2) It is unclear if there are any auto-regressive effects
##      across subjects. 
## 3) There is clearly ONE common factor across all subjects.
## 4) It seems there aren't more common factors. 
## 5) 
##
##
##
##
##
  

## Idea: examine communities within the inner product matrix:

par(mfrow = c(2,2))
for(i in 1:4) {
  pc_data <- df_list %>% 
    sapply(
      function(x) {
        pca <- prcomp(x)#, center = TRUE, scale. = TRUE) 
        return(pca$rotation[, i])
      }) %>% 
    
    t()  
  
  for (j in 1:16) {
    pc_data[j, ] <- pc_data[j,] * sign(pc_data[j,1])
  }
  rownames(pc_data) <- 1:16

  clustering <- kmeans(
    pc_data, centers = 3, nstart = 50)
  print(paste0("PC", i, ": clusters"))
  print(clustering$cluster)
  
  cluster_order <- lapply(
    1:3,
    function(x) {
      which(clustering$cluster == x)
    }) %>%
    unlist()
  
  pc_data <- pc_data[, ]
  ipi_mat <- matrix(0, nrow(pc_data), nrow(pc_data))
  for(i in 1:nrow(pc_data)) {
    for(j in 1:nrow(pc_data)) {
      ipi_mat[i,j] <- sum(pc_data[cluster_order[i], ] * pc_data[cluster_order[j], ])    
    }
  }
  plot(abs(ipi_mat), 
       breaks = seq(0,1, length.out = 20))

}






