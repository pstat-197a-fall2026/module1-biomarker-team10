## Question 4, part b: 
## Find an alternative panel that achieves improved classification accuracy 
## and benchmark results against the in-class analysis.

## We will be using a random forest model since from previous experience, it 
## provides the most accurate predictions with the trade off being 
## time and space constraints.

## Since our in class analysis used random forest as a conjunction with 
## logistic regression, we ultimately chose the proteins that intersected.

## In other words, part of the model was attributed to logistic regression.
## Our improved panel will only incorporate random forest models and 
## use it iteratively through multiple folds.

library(tidyverse)
library(tidymodels)
library(randomForest)
library(yardstick)

load("data/biomarker-clean.RData")

# Make ASD the positive class
biomarker_clean <- biomarker_clean %>%
  mutate(group = factor(group, levels = c("TD", "ASD")))

# 1. Split BEFORE selecting any proteins
# Reserve 20% for the final test set. Do not use it to choose proteins or tune the model.
set.seed(101422)

biomarker_split <- initial_split(
  biomarker_clean,
  prop = 0.8,
  strata = group
)

train_data <- training(biomarker_split)
test_data  <- testing(biomarker_split)

# Exclude the outcome and ADOS score from predictors
protein_names <- setdiff(
  names(train_data),
  c("group", "ados")
)

# Candidate panel sizes
panel_sizes <- c(5, 10, 15, 20, 30, 50)

# Keep only sizes supported by the available proteins
panel_sizes <- panel_sizes[
  panel_sizes <= length(protein_names)
]

# 2. Cross-validation on the training data
# For each cross-validation fold, train a random forest and rank proteins by MeanDecreaseGini. 
# Repeat with different panel sizes: 5, 10, 15, 20, 30, and 50 proteins.

set.seed(101422)

folds <- vfold_cv(
  train_data,
  v = 5,
  strata = group
)

cv_results <- list()

for (k in panel_sizes) {
  
  fold_accuracy <- numeric(nrow(folds))
  
  for (i in seq_len(nrow(folds))) {
    
    fold_train <- analysis(folds$splits[[i]])
    fold_valid <- assessment(folds$splits[[i]])
    
    # Rank proteins using ONLY this fold's training data
    rf_rank <- randomForest(
      x = fold_train[, protein_names, drop = FALSE],
      y = fold_train$group,
      ntree = 1000,
      importance = TRUE
    )
    
    importance_scores <- importance(rf_rank)[
      , "MeanDecreaseGini"
    ]
    
    ranked_proteins <- names(
      sort(importance_scores, decreasing = TRUE)
    )
    
    selected_proteins <- head(ranked_proteins, k)
    
    # Fit the classifier using the selected panel
    rf_panel <- randomForest(
      x = fold_train[, selected_proteins, drop = FALSE],
      y = fold_train$group,
      ntree = 1000
    )
    
    # Evaluate on this fold's validation data
    predictions <- predict(
      rf_panel,
      newdata = fold_valid[, selected_proteins, drop = FALSE]
    )
    
    fold_accuracy[i] <- mean(
      predictions == fold_valid$group
    )
  }
  
  cv_results[[as.character(k)]] <- tibble(
    n_proteins = k,
    mean_accuracy = mean(fold_accuracy),
    sd_accuracy = sd(fold_accuracy)
  )
}

cv_summary <- bind_rows(cv_results) %>%
  arrange(desc(mean_accuracy))

print(cv_summary)

# Select the panel size with the highest CV accuracy
best_k <- cv_summary$n_proteins[1]

cat("Selected panel size:", best_k, "\n")

# 3. Re-rank proteins using the FULL training set
# Re-rank proteins using the full training set, select the best panel size 
# from cross-validation, and fit the model using that panel.

set.seed(101422)

rf_rank_final <- randomForest(
  x = train_data[, protein_names, drop = FALSE],
  y = train_data$group,
  ntree = 1000,
  importance = TRUE
)

importance_scores <- importance(rf_rank_final)[
  , "MeanDecreaseGini"
]

ranked_proteins <- names(
  sort(importance_scores, decreasing = TRUE)
)

alternative_panel <- head(ranked_proteins, best_k)

cat("Alternative protein panel:\n")
print(alternative_panel)

# 4. Fit the final random forest on training data
# Rank proteins using the full training set, select the chosen number of proteins, 
# fit the final random forest, and evaluate it on the untouched test set.
set.seed(101422)

rf_final <- randomForest(
  x = train_data[, alternative_panel, drop = FALSE],
  y = train_data$group,
  ntree = 1000,
  importance = TRUE
)

# Evaluate on the untouched test data
test_predictions <- predict(
  rf_final,
  newdata = test_data[, alternative_panel, drop = FALSE]
)

test_probabilities <- predict(
  rf_final,
  newdata = test_data[, alternative_panel, drop = FALSE],
  type = "prob"
)[, "ASD"]

# Classification metrics
# Measure accuracy, sensitivity, specificity, and ROC AUC on the untouched test set.

confusion_matrix <- table(
  Truth = test_data$group,
  Predicted = test_predictions
)

accuracy <- mean(test_predictions == test_data$group)

sensitivity <- sens_vec(
  truth = test_data$group,
  estimate = test_predictions,
  event_level = "second"
)

specificity <- spec_vec(
  truth = test_data$group,
  estimate = test_predictions,
  event_level = "second"
)

auc <- roc_auc_vec(
  truth = test_data$group,
  estimate = test_probabilities,
  event_level = "second"
)

print(confusion_matrix)

alternative_results <- tibble(
  Model = "Alternative Random Forest",
  Proteins = length(alternative_panel),
  Accuracy = accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  ROC_AUC = auc
)

print(alternative_results)

# Save the selected panel
write.csv(
tibble(protein = alternative_panel),
  "results/alternative-panel.csv",
  row.names = FALSE
)