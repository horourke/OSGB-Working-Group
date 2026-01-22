library(glmnet)

# 2. Prepare Data (Example using built-in 'swiss' dataset)
data(swiss)
# Convert factors to dummy variables if needed (glmnet works with matrices)
x_vars <- model.matrix(Fertility ~ ., data = swiss)[, -1] # Predictors, drop intercept
y_var <- swiss$Fertility # Response

# 3. Split Data (Train/Test)
set.seed(123) # for reproducibility
train_indices <- sample(1:nrow(x_vars), nrow(x_vars) * 0.7)
x_train <- x_vars[train_indices, ]
y_train <- y_var[train_indices]
x_test <- x_vars[-train_indices, ]
y_test <- y_var[-train_indices]

# 4. Fit LASSO model (Cross-validation to find best lambda)
# alpha = 1 specifies LASSO (alpha = 0 is Ridge, 0 < alpha < 1 is Elastic Net)
# family = "gaussian" for continuous outcome (use "binomial" for binary)
lasso_model <- cv.glmnet(x = x_train, y = y_train, alpha = 1, family = "gaussian", nfolds = 10)

# 5. Plot the Cross-Validation Results (shows error vs. lambda)
plot(lasso_model)


# 6. Get the best lambda (lambda.1se is common choice for parsimony)
best_lambda <- lasso_model$lambda.1se
print(paste("Best Lambda:", best_lambda))

# 7. Inspect Coefficients at Best Lambda (shows feature selection)
# glmnet automatically standardizes predictors
coef(lasso_model, s = best_lambda)

# 8. Make Predictions on Test Data
predictions <- predict(lasso_model, s = best_lambda, newx = x_test)

# 9. Evaluate Model Performance (e.g., MSE)
mean((y_test - predictions)^2)