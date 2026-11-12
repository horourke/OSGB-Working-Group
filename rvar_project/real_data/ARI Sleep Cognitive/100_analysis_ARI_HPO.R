

library(magrittr)
library(tidyverse)

source("./002_requirements_lite.R")

source("./041_rvar_supps_vect.R")
source("./043_lm_rvar_vectorized.r")
source("./044_glmnet_rvar_vectorized.r")
source("./045_cv.glmnet_rvar_vectorized.r")
source("./046_bic.glmnet_rvar_vectorized.r")
source("./047_glmnet_rvar_adaptive.r")

library(openxlsx)
ari_data <- read.xlsx("bigger_daily_ass.xlsx")

dim(ari_data)

ari_data <- ari_data %>%
  mutate(across(where(is.character), na_if, "NA"))

# create subset for list with 3 Y variables & 2 X variables

##### EDIT
# X = gender, employed
# Y = depression, anxiety, subjective_sleep

# Code gender as 0 = female, 1 = male
ari_data$f0m1 <- ari_data$gender - 1


# View time points by individual
table(ari_data$id)

# Get variable names
ls(ari_data)

# Rename variables to be more intuitive 
ari_data <- ari_data %>%
  rename(workingmem = va_acc,
         attn_intensity = m.rt,
         attn_consistency = cv.rt,
         exec_fn = sum_flanker)

# compute sleep variables
ari_data <- ari_data %>%
  mutate(
    rem_prop = rem / sleep_time,
    deep_prop = deep / sleep_time
  )

ls(ari_data)

# Aggregate data at individual level
ind <- ari_data %>%
  group_by(id) %>%
  summarize(across(everything(), mean, na.rm = TRUE)
  )

# Examine correlations among mood variables w/ all data
mood <- data.frame(ari_data$depression, ari_data$anxiety, ari_data$distress, ari_data$somatization)
cormood <- round(cor(mood, use = "complete.obs"), 3)
print(cormood)

# Examine correlations among mood variables w/ subject-level data
mood2 <- data.frame(ind$depression, ind$anxiety, ind$distress, ind$somatization)
cormood2 <- round(cor(mood2, use = "complete.obs"), 3)
print(cormood2)

# Examine correlations among sleep variables w/ all data
sleep <- data.frame(ari_data$subjective_sleep, ari_data$awake_time, ari_data$sleep_time, 
                    ari_data$light, ari_data$deep, ari_data$rem)
corsleep <- round(cor(sleep, use = "complete.obs"), 3)
print(corsleep)

# Examine correlations among sleep variables w/ subject-level data
sleep2 <- data.frame(ind$subjective_sleep, ind$awake_time, ind$sleep_time, 
                     ind$light, ind$deep, ind$rem)
corsleep2 <- round(cor(sleep2, use = "complete.obs"), 3)
print(corsleep2)

# Examine correlations among big 5 variables for all data
big5 <- data.frame(ari_data$agreeableness, ari_data$conscientiousness, ari_data$emotional_stability, 
                   ari_data$extraversion, ari_data$openness)
corbig5 <- round(cor(big5, use = "complete.obs"), 3)
print(corbig5)

# Examine correlations among big 5 variables w/ subject-level data
big52 <- data.frame(ind$agreeableness, ind$conscientiousness, ind$emotional_stability, 
                    ind$extraversion, ind$openness)
corbig52 <- round(cor(big52, use = "complete.obs"), 3)
print(corbig52)

pairs(big52)

# Examine correlations among cognitive variables for all data
cog <- data.frame(ari_data$cv.rt, ari_data$m.rt, 
                  ari_data$va_acc, 
                  ari_data$total_flanker)
corcog <- round(cor(cog, use = "complete.obs"), 3)
print(corcog)

# Examine correlations among cognitive variables w/ subject-level data
cog2 <- data.frame(ind$cv.rt, ind$m.rt, 
                   ind$va_acc, 
                   ind$total_flanker)
corcog2 <- round(cor(cog2, use = "complete.obs"), 3)
print(corcog2)
pairs(cog2)

# create sequential observation indices for each participant
ari_data <- ari_data %>%
  arrange(id, true_date) %>%                # Sort chronologically within each participant
  group_by(id) %>%
  mutate(time_seq = row_number()) %>%  # Assigns 1, 2, 3... per participant
  ungroup()

# create lagged variables for each set 
lagged_data <- ari_data %>%
  arrange(id, time_seq) %>%
  group_by(id) %>%
  mutate(
    across(
      where(
        is.numeric
      ), 
      \(x) lag(x, n = 1),
      .names = "{.col}_l"
    )
  ) %>%
  ungroup()

dep <- cbind(lagged_data$depression, lagged_data$depression_l)
cor(lagged_data$m.rt, lagged_data$m.rt_l, use = "complete.obs")

# Aggregate lagged data at individual level
ind_l <- lagged_data %>%
  group_by(id) %>%
  summarize(across(everything(), mean, na.rm = TRUE)
  )

# Examine lagged correlations among mood variables w/ all data
mood_l <- data.frame(lagged_data$depression, lagged_data$anxiety, lagged_data$distress, lagged_data$somatization,
                     lagged_data$depression_l, lagged_data$anxiety_l, lagged_data$distress_l, lagged_data$somatization_l)
print(round(cor(mood_l, use = "complete.obs"), 3))
pairs(mood_l)

# Examine lagged correlations among mood variables w/ subject-level data
mood2_l <- data.frame(ind_l$depression, ind_l$anxiety, ind_l$distress, ind_l$somatization,
                    ind_l$depression_l, ind_l$anxiety_l, ind_l$distress_l, ind_l$somatization_l)
print(round(cor(mood2_l, use = "complete.obs"), 3))
pairs(mood2_l)

# Create data matrix for X's (covariates)

ari_x <- ari_data %>%
  group_by(id) %>%
  summarize(
    f0m1 = mean(f0m1, na.rm = TRUE),
    employed = mean(employed, na.rm = TRUE),
    id = n()
  )

ari_x <- ari_x %>% 
  select(f0m1, employed)


#Set up data for Y
# Check date codes; they do not code time points
table(ari_data$true_date)

# Code time as variable t 
ari_time <- ari_data %>%
  arrange(id, true_date) %>% # It's good practice to arrange by ID and Time first
group_by(id) %>%
  mutate(t = row_number()) %>%
  ungroup()

# subset data to contain only id, t, and Y vars

yvars <- c("depression", "anxiety", "distress", "somatization",
           "agreeableness", "conscientiousness", "emotional_stability", "openness",
           "workingmem", "attn_intensity", "attn_consistency", "exec_fn",
           "hrv", "hrv_deep", "o2_mean", 
           "sleep_time", "rem_prop", "deep_prop", "time_awake", "subjective_sleep")

ari_y <- ari_time %>%
  select(id, t, depression, anxiety, subjective_sleep)
head(ari_y)

# Select complete cases only

ari_y_complete <- ari_y[complete.cases(ari_y), ]
head(ari_y_complete)


#Create list of matrices for Y's (outcomes)

id_vec <- unique(ari_y_complete$id)

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

#MOD-VAR using gender + employment


library(glmnet)
library(plot.matrix)

cd_var_output <- cv.solve_rvar_glmnet_vectorized(
  d = 3, p = 2,
  Y_list = ari_y_list, 
  X = ari_x,
  #### TUNE LAMBDA AND PENALTY FACTOR
  lambda.seq = 10^seq(-12, 0, length.out = 50), 
  penalty.factor = 10^seq(-15, 15, length.out = 50),  
  ####
  nfolds = 3, 
  verbose = T)

library(plot.matrix)
library("psych")

# Overall plot 

layout(mat = matrix(c(1,1,2,3,4,5), nrow = 3, ncol = 2, byrow = TRUE))
par(mar = c(5.1, 4.1, 4.1, 3.1)) 

ymax <- max(abs(cd_var_output$rvar_opt_coeffs[, - (1:3)]))

rwb <- colorRampPalette(colors = c("red", "white", "blue")) 
rwb_vec <- rwb(12)

cd_var_output$rvar_opt_coeffs[, - (1:3)] %>%
  as.matrix() %>%
  #abs() %>%
  plot(breaks = seq(-ymax,ymax, length.out = 12), col = rwb_vec,
       digits = 2, text.cell = list(cex=1))

# Moderation plots

pwg <- colorRampPalette(colors = c("deeppink1", "white", "darkgreen")) 
pwg_vec <- pwg(12)

# According to Jose: rows "var" are y_t, columns are y_t-1
one <- as.matrix(cd_var_output$rvar_opt_coeffs[, -c(1:3)][1:3]) 
two <- as.matrix(cd_var_output$rvar_opt_coeffs[, -c(1:3)][7:9])
three <- as.matrix(cd_var_output$rvar_opt_coeffs[, -c(1:3)][4:6])

ymax <- max(abs(rbind(
    one,
    one + two,
    one + three,
    one + two + three)))

one %>%
  plot(main = "Unemployed female", 
        breaks = seq(-ymax,ymax, length.out = 12), col = pwg_vec, 
       digits = 2, text.cell = list(cex=1))

plot2 <- one + three

plot2 %>%
  plot(main = "Unemployed male", 
        breaks = seq(-ymax,ymax, length.out = 12), col = pwg_vec, 
       digits = 2, text.cell = list(cex=1))

plot3 <- one + two 

plot3  %>%
  plot(main = "Employed female", 
        breaks = seq(-ymax,ymax, length.out = 12), col = pwg_vec, 
       digits = 2, text.cell = list(cex=1))

plot4 <- one + three + two

plot4  %>%
  plot(main = "Employed male", 
        breaks = seq(-ymax,ymax, length.out = 12), col = pwg_vec, 
       digits = 2, text.cell = list(cex=1))

