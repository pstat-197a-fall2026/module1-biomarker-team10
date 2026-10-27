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

# Alternative biomarker panel: Random Forest
# Fair comparison against the in-class logistic regression approach.
# Feature selection and tuning use training data only; the test set is used once.

library(tidyverse)
library(tidymodels)
library(randomForest)
library(yardstick)

load("data/biomarker-clean.RData")

# Ensure the positive/event class is ASD.
biomarker_clean <- biomarker_clean %>%
  mutate(group = factor(group, levels = c("TD", "ASD")))

# Keep ADOS out of the predictors, matching the in-class analysis.
protein_names <- setdiff(names(biomarker_clean), c("group", "ados"))

# -------------------------------------------------------------------
# 1. Split FIRST so test data do not affect feature selection or tuning.
# -------------------------------------------------------------------
set.seed(101422)
biomarker_split <- initial_split(
  biomarker_clean,
  prop = 0.8,
  strata = group
)

train_data <- training(biomarker_split)
test_data  <- testing(biomarker_split)

# -------------------------------------------------------------------
# 2. In-class baseline, adapted to select proteins using training data.
#    Top 10 adjusted t-test results intersected with top 10 RF proteins.
# -------------------------------------------------------------------
ttest_p <- sapply(protein_names, function(p) {
  t.test(train_data[[p]] ~ train_data$group)$p.value
})

# Use a Benjamini-Hochberg adjustment for the training-only t-tests.
ttest_adj <- p.adjust(ttest_p, method = "BH")
proteins_s1 <- names(sort(ttest_adj))[seq_len(min(10, length(ttest_adj)))]

set.seed(101422)
rf_baseline_rank <- randomForest(
  x = train_data[, protein_names, drop = FALSE],
  y = train_data$group,
  ntree = 1000,
  importance = TRUE
)

baseline_importance <- importance(rf_baseline_rank)[, "MeanDecreaseGini"]
proteins_s2 <- names(sort(baseline_importance, decreasing = TRUE))[
  seq_len(min(10, length(baseline_importance)))
]

baseline_train <- train_data %>%
  select(group, all_of(baseline_panel)) %>%
  mutate(class = group == "ASD") %>%
  select(-group)

baseline_test <- test_data %>%
  select(group, all_of(baseline_panel)) %>%
  mutate(class = group == "ASD") %>%
  select(-group)

baseline_fit <- glm(class ~ ., data = baseline_train, family = binomial())
baseline_prob <- predict(baseline_fit, newdata = baseline_test, type = "response")
baseline_pred <- factor(ifelse(baseline_prob >= 0.5, "ASD", "TD"),
                        levels = c("TD", "ASD"))

baseline_truth <- factor(
  test_data$group,
  levels = c("TD", "ASD")
)

baseline_results <- tibble(
  Model = "In-class Logistic Regression (training-only selection)",
  Proteins = length(baseline_panel),
  Accuracy = accuracy_vec(baseline_truth, baseline_pred),
  Sensitivity = sens_vec(baseline_truth, baseline_pred, event_level = "second"),
  Specificity = spec_vec(baseline_truth, baseline_pred, event_level = "second"),
  ROC_AUC = roc_auc_vec(baseline_truth, baseline_prob, event_level = "second")
)

# -------------------------------------------------------------------
# 3. Alternative Random Forest: cross-validate panel size and mtry.
#    Protein rankings are recalculated within each fold using fold-training
#    data only. Accuracy is used to choose the panel size and mtry.
# -------------------------------------------------------------------
panel_sizes <- c(5, 10, 15, 20, 30, 50)
panel_sizes <- panel_sizes[panel_sizes <= length(protein_names)]

set.seed(101422)
folds <- vfold_cv(train_data, v = 5, strata = group)

cv_results <- list()
result_index <- 1

for (k in panel_sizes) {
  # Try a small set of mtry values for each panel size.
  mtry_values <- unique(pmax(1, pmin(k, c(floor(sqrt(k)), floor(k / 3), k))))
  
  for (fold_i in seq_len(nrow(folds))) {
    fold_train <- analysis(folds$splits[[fold_i]])
    fold_valid <- assessment(folds$splits[[fold_i]])
    
    # Rank proteins using only this fold's training observations.
    set.seed(101422 + fold_i)
    rf_rank <- randomForest(
      x = fold_train[, protein_names, drop = FALSE],
      y = fold_train$group,
      ntree = 500,
      importance = TRUE
    )
    fold_importance <- importance(rf_rank)[, "MeanDecreaseGini"]
    ranked_proteins <- names(sort(fold_importance, decreasing = TRUE))
    selected_proteins <- head(ranked_proteins, k)
    
    for (m in mtry_values) {
      set.seed(20000 + fold_i * 100 + k + m)
      rf_fit <- randomForest(
        x = fold_train[, selected_proteins, drop = FALSE],
        y = fold_train$group,
        ntree = 500,
        mtry = m
      )
      
      fold_pred <- predict(
        rf_fit,
        newdata = fold_valid[, selected_proteins, drop = FALSE]
      )
      
      cv_results[[result_index]] <- tibble(
        n_proteins = k,
        mtry = m,
        fold = fold_i,
        accuracy = mean(fold_pred == fold_valid$group)
      )
      result_index <- result_index + 1
    }
  }
}

cv_detail <- bind_rows(cv_results)
cv_summary <- cv_detail %>%
  group_by(n_proteins, mtry) %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy = sd(accuracy),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_accuracy), n_proteins)

print(cv_summary)
best_setting <- cv_summary %>% slice(1)
best_k <- best_setting$n_proteins[[1]]
best_mtry <- best_setting$mtry[[1]]

cat("Selected panel size:", best_k, "\n")
cat("Selected mtry:", best_mtry, "\n")

# -------------------------------------------------------------------
# 4. Re-rank on all training data, fit the final RF, evaluate once on test.
# -------------------------------------------------------------------
set.seed(101422)
rf_rank_final <- randomForest(
  x = train_data[, protein_names, drop = FALSE],
  y = train_data$group,
  ntree = 1000,
  importance = TRUE
)
final_importance <- importance(rf_rank_final)[, "MeanDecreaseGini"]
ranked_final <- names(sort(final_importance, decreasing = TRUE))
alternative_panel <- head(ranked_final, best_k)

cat("Alternative protein panel:\n")
print(alternative_panel)

set.seed(101422)
rf_final <- randomForest(
  x = train_data[, alternative_panel, drop = FALSE],
  y = train_data$group,
  ntree = 1000,
  mtry = best_mtry,
  importance = TRUE
)

rf_test_pred <- predict(rf_final, newdata = test_data[, alternative_panel, drop = FALSE])
rf_test_prob <- predict(
  rf_final,
  newdata = test_data[, alternative_panel, drop = FALSE],
  type = "prob"
)[, "ASD"]

test_truth <- factor(test_data$group, levels = c("TD", "ASD"))
rf_test_pred <- factor(rf_test_pred, levels = c("TD", "ASD"))

alternative_results <- tibble(
  Model = "Tuned Random Forest",
  Proteins = length(alternative_panel),
  Accuracy = accuracy_vec(test_truth, rf_test_pred),
  Sensitivity = sens_vec(test_truth, rf_test_pred, event_level = "second"),
  Specificity = spec_vec(test_truth, rf_test_pred, event_level = "second"),
  ROC_AUC = roc_auc_vec(test_truth, rf_test_prob, event_level = "second")
)

# -------------------------------------------------------------------
# 5. Compare on the same held-out test subjects and save outputs.
# -------------------------------------------------------------------
comparison <- bind_rows(baseline_results, alternative_results)
print(comparison)

#write.csv(tibble(protein = alternative_panel),
#          "results/alternative-panel.csv", row.names = FALSE)
#write.csv(cv_summary, "results/random-forest-cv-summary.csv", row.names = FALSE)
#write.csv(comparison, "results/model-comparison.csv", row.names = FALSE)