## Q4 draft: alternative panel via lasso on a pooled candidate set
## Suggested location: scripts/drafts/q4-alternative-panel.R
library(tidyverse)
library(randomForest)
library(tidymodels)
library(glmnet)
load('data/biomarker-clean.RData')

## SPLIT FIRST (same seed / same non-stratified split as in-class)
##################################################################
set.seed(101422)
biomarker_all <- biomarker_clean %>%
  select(-ados) %>%
  mutate(class = factor(group == 'ASD', levels = c(FALSE, TRUE))) %>%
  select(-group)

split <- initial_split(biomarker_all, prop = 0.8)
train <- training(split)
test  <- testing(split)

prot_names <- setdiff(names(train), 'class')

## SELECTION ON TRAINING DATA ONLY
###################################
k <- 25  # candidates per method

# method 1: t-test p-values
pvals <- map_dbl(prot_names, ~ t.test(train[[.x]] ~ train$class)$p.value)
top_t <- prot_names[order(pvals)][1:k]

# method 2: random forest importance
set.seed(101422)
rf <- randomForest(x = train[, prot_names], y = train$class,
                   ntree = 1000, importance = TRUE)
top_rf <- names(sort(rf$importance[, 'MeanDecreaseGini'],
                     decreasing = TRUE))[1:k]

# hard-intersection baseline (train-only) and pooled candidate set
panel_baseline <- intersect(top_t[1:10], top_rf[1:10])
candidates <- union(top_t, top_rf)

## ALTERNATIVE PANEL: CROSS-VALIDATED LASSO ON CANDIDATES
##########################################################
x_train <- as.matrix(train[, candidates])
y_train <- train$class

set.seed(101422)
cv_fit <- cv.glmnet(x_train, y_train, family = 'binomial',
                    alpha = 1, type.measure = 'auc', nfolds = 10)

coefs <- coef(cv_fit, s = 'lambda.1se')
panel_lasso <- setdiff(rownames(coefs)[as.vector(coefs != 0)], '(Intercept)')

## EVALUATION (same test set for every panel)
##############################################
class_metrics <- metric_set(sensitivity, specificity, accuracy, roc_auc)

evaluate_panel <- function(panel, label) {
  if (length(panel) == 0) stop('Empty panel: ', label)
  fit <- glm(class ~ ., family = 'binomial',
             data = train %>% select(class, all_of(panel)))
  test %>%
    mutate(pred = predict(fit, newdata = test, type = 'response'),
           est = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
    class_metrics(truth = class, estimate = est, pred,
                  event_level = 'second') %>%
    mutate(panel = label, n_proteins = length(panel))
}

results <- bind_rows(
  evaluate_panel(panel_baseline, 'in-class logic (train-only selection)'),
  evaluate_panel(panel_lasso,    'lasso panel')
)
print(results)
print(panel_lasso)

## TODO: repeat over many random splits (or nested CV) to get a spread of
## each metric, since one ~30-subject test set is noisy.